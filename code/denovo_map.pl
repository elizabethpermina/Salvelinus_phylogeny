#!/usr/bin/env perl
#
# Copyright 2010-2017, Julian Catchen <jcatchen@illinois.edu>
#
# This file is part of Stacks.
#
# Stacks is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Stacks is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Stacks.  If not, see <http://www.gnu.org/licenses/>.
#

#
# Process the data for a genetic map: build stacks in parents and progeny,
# create a catalog from the parents, and match progeny against the catatlog.
# Call genotypes, and load all data into an MySQL database along the way.
#
# For the database interactions to work, the 'mysql' program is expected to be
# on the path and sufficient permissions set to access the specified database.
#

use strict;
use POSIX;
use File::Temp qw/ mktemp /;
use File::Spec;
use constant stacks_version => "2.0Beta7c";

use constant true  => 1;
use constant false => 0;

my $dry_run      = false;
my $sql          = false;
my $create_db    = false;
my $overw_db     = false;
my $gapped_alns  = false;
my $mysql_config = "/home/epermina/soft/stacks-2.0Beta7c/share/stacks/" . "sql/mysql.cnf";
my $mysql_tables = "/home/epermina/soft/stacks-2.0Beta7c/share/stacks/" . "sql/stacks.sql";
my $exe_path     = "/home/epermina/soft/stacks-2.0Beta7c/bin/";
my $out_path     = "";
my $popmap_path  = "";
my $sample_path  = "";
my $db           = "";
my $data_type    = "map";
my $min_cov      = 0;
my $min_rcov     = 0;
my $batch_id     = 1;
my $sample_id    = 1;
my $desc         = ""; # Database description of this dataset
my $date         = ""; # Date relevent to this data, formatted for SQL: 2009-05-31
my $gzip         = false;
my $paired       = false;
my $time         = "";

my @parents;
my @progeny;
my @samples;

my (@_ustacks, @_cstacks, @_sstacks, @_tsv2bam, @_gstacks, @_genotypes, @_populations);

my $cmd_str = $0 . " " . join(" ", @ARGV);

parse_command_line();

my $cnf = (-e $ENV{"HOME"} . "/.my.cnf") ? $ENV{"HOME"} . "/.my.cnf" : $mysql_config;

#
# Check for the existence of the necessary pipeline programs
#
foreach my $prog ("ustacks", "cstacks", "sstacks", "tsv2bam", "gstacks", "populations", "index_radtags.pl") {
    die "Unable to find '" . $exe_path . $prog . "'.\n" if (!-e $exe_path . $prog);
}

my ($log, $log_fh, $sample);

my (@sample_list, %pop_ids, %pops, %grp_ids, %grps, %sample_ids);

parse_population_map(\@sample_list, \%pop_ids, \%pops, \%grp_ids, \%grps);

initialize_samples(\@parents, \@progeny, \@samples, \@sample_list, \%pop_ids, \%grp_ids);

#
# Open the log file
#
$log = "$out_path/denovo_map.log";
open($log_fh, ">$log") or die("Unable to open log file '$log'; $!\n");

print $log_fh
        "denovo_map.pl version ", stacks_version, " started at ", strftime("%Y-%m-%d %H:%M:%S", (localtime(time))), "\n",
        $cmd_str, "\n";

initialize_database($log_fh, \@parents, \@progeny, \@samples, \%sample_ids) if ($sql == true);

execute_stacks($log_fh, $sample_id, \@parents, \@progeny, \@samples, \%sample_ids);

load_sql_data($log_fh, \%pops, \@parents, \@progeny, \@samples) if ($sql == true);

print $log_fh "\ndenovo_map.pl completed at ", strftime("%Y-%m-%d %H:%M:%S", (localtime(time))), "\n";
close($log_fh);

sub check_return_value {
    #
    # $? is a 16 bit int. Exit code is given by `$? & 255` if the process was
    # terminated by a signal, and by `$? >> 8` if it exited normally.
    #
    my ($rv, $log_fh) = @_;
    if ($rv != 0) {
        my $msg = "\ndenovo_map.pl: Aborted because the last command failed. (" . ($rv & 255 ? $rv : 128 + ($rv >> 8) % 128) . ")\n";
        print $log_fh $msg;
        print STDERR $msg;
        exit 1;
    }
}

sub execute_stacks {
    my ($log_fh, $sample_id, $parents, $progeny, $samples, $sample_ids) = @_;

    my (@results, @depths_of_cov);
    my ($pop_cnt, $sample, $num_files, $i, $cmd, $pipe_fh, $cat_file);

    my $minc  = $min_cov  > 0 ? " -m $min_cov"  : "";
    my $minrc = $min_rcov > 0 ? " -m $min_rcov" : $minc;

    $i         = 1;
    $num_files = scalar(@{$parents}) + scalar(@{$progeny}) + scalar(@{$samples});

    #
    # Assemble RAD loci in each individual.
    #
    print STDERR "Indentifying unique stacks...\n";
    print $log_fh "\nustacks\n==========\n";
    foreach $sample (@parents, @progeny, @samples) {

        print $log_fh "\nSample $i of $num_files '$sample->{'file'}'\n----------\n";

        if (scalar(keys %{$sample_ids}) > 0) {
            $sample_id = $sample_ids->{$sample->{'file'}};
        }

        $cmd = $exe_path . "ustacks -t $sample->{'fmt'} -f $sample->{'path'} -o $out_path -i $sample_id";
        if ($sample->{'type'} eq "sample") {
            $cmd .= $minc;
        } elsif ($sample->{'type'} eq "parent") {
            $cmd .= $minc;
        } elsif ($sample->{'type'} eq "progeny") {
            $cmd = $minrc;
        }
        if ($sample->{'path'} !~ /$sample->{'file'}.$sample->{'suffix'}$/) {
            # Guessing the sample name from the input path won't work.
            $cmd .= " --name " . $sample->{'file'};
        }
        foreach (@_ustacks) {
            $cmd .= " " . $_;
        }
        print STDERR  "  $cmd\n";
        print $log_fh "$cmd\n";

        if ($dry_run == false) {
            open($pipe_fh, "$time $cmd 2>&1 |");
            @results = ();
            while (<$pipe_fh>) {
                print $log_fh $_;
                push(@results, $_);
            }
            close($pipe_fh);
            check_return_value($?, $log_fh);

            #
            # Pull the depth of coverage from ustacks.
            #
            my $depthline = (grep(/^Final coverage/, @results))[0];
            my ($depth) = ($depthline =~ /mean=([^;]+)/);
            push(@depths_of_cov, [$sample->{'file'}, $depth]);
        }

        $i++;
        $sample_id++;
    }

    write_depths_of_cov(\@depths_of_cov, $log_fh);
    print STDERR "\n";

    #
    # Generate catalog of RAD loci.
    #
    print STDERR "Generating catalog...\n";
    print $log_fh "\ncstacks\n==========\n";

    $cmd = $exe_path . "cstacks -P $out_path " . join(" ", @_cstacks);
    print STDERR  "  $cmd\n\n";
    print $log_fh "$cmd\n\n";

    if ($dry_run == false) {
        open($pipe_fh, "$time $cmd 2>&1 |");
        while (<$pipe_fh>) {
            print $log_fh $_;
            if ($_ =~ /failed/i) { print STDERR "Catalog construction failed.\n"; exit(1); }
        }
        close($pipe_fh);
        check_return_value($?, $log_fh);
    }

    #
    # Match parents, progeny, or samples to the catalog.
    #
    print STDERR "Matching samples to the catalog...\n";
    print $log_fh "\nsstacks\n==========\n";

    $cat_file = "batch_" . $batch_id;
    $cmd      = $exe_path . "sstacks -P $out_path " . join(" ", @_sstacks);
    print STDERR  "  $cmd\n\n";
    print $log_fh "$cmd\n\n";

    if ($dry_run == false) {
        open($pipe_fh, "$time $cmd 2>&1 |");
        while (<$pipe_fh>) {
            print $log_fh $_;
        }
        close($pipe_fh);
        check_return_value($?, $log_fh);
    }

    #
    # Sort the reads according by catalog locus / run tsv2bam.
    #
    print STDERR "Sorting reads by RAD locus...\n";
    print $log_fh "\ntsv2bam\n==========\n";

    my $file_paths = "";
    if (length($popmap_path) == 0) {
        foreach $sample (@parents, @progeny, @samples) {
            $file_paths .= " -s $sample->{'file'}";
        }
    }

    $cmd = $exe_path . "tsv2bam -P $out_path $file_paths " . join(" ", @_tsv2bam);
    print STDERR  "  $cmd\n";
    print $log_fh "$cmd\n";

    if (!$dry_run) {
        open($pipe_fh, "$time $cmd 2>&1 |");
        while (<$pipe_fh>) {
            print $log_fh $_;
        }
        close($pipe_fh);
        check_return_value($?, $log_fh);
    }

    #
    # Call genotypes / run gstacks.
    #
    print STDERR "Calling variants, genotypes and haplotypes...\n";
    print $log_fh "\ngstacks\n==========\n";

    $cmd = $exe_path . "gstacks -P $out_path " . join(" ", @_gstacks);
    print STDERR  "  $cmd\n\n";
    print $log_fh "$cmd\n\n";

    if (!$dry_run) {
        open($pipe_fh, "$time $cmd 2>&1 |");
        while (<$pipe_fh>) {
            print $log_fh $_;
        }
        close($pipe_fh);
        check_return_value($?, $log_fh);
    }

    if ($data_type eq "map") {
        #
        # Generate a set of observed haplotypes and a set of markers and generic genotypes
        #
        printf(STDERR "Generating genotypes...\n");
        print $log_fh "\ngenotypes\n==========\n";

        $cmd = $exe_path . "genotypes" . " -b $batch_id -P $out_path -r 1 -c -s " . join(" ", @_genotypes);
        print STDERR  "$cmd\n\n";
        print $log_fh "$cmd\n\n";

        if ($dry_run == 0) {
            open($pipe_fh, "$time $cmd 2>&1 |");
            while (<$pipe_fh>) {
                print $log_fh $_;
            }
            close($pipe_fh);
            check_return_value($?, $log_fh);
        }

    } else {
        printf(STDERR "Calculating population-level summary statistics\n");
        print $log_fh "\npopulations\n==========\n";

        $cmd = $exe_path . "populations" . " -P $out_path " . join(" ", @_populations);
        print STDERR  "  $cmd\n\n";
        print $log_fh "$cmd\n\n";

        if ($dry_run == 0) {
            open($pipe_fh, "$time $cmd 2>&1 |");
            while (<$pipe_fh>) {
                print $log_fh $_;
            }
            close($pipe_fh);
            check_return_value($?, $log_fh);
        }
    }
    
    print STDERR  "denovo_map.pl is done.\n";
    print $log_fh "denovo_map.pl is done.\n";
}

sub parse_population_map {
    my ($sample_list, $pop_ids, $pops, $grp_ids, $grps) = @_;

    my ($fh, @parts, $line, $sample);

    return if (length($popmap_path) == 0);

    open($fh, "<$popmap_path") or die("Unable to open population map, '$popmap_path', $!\n");

    while ($line = <$fh>) {
        chomp $line;
        next if ($line =~ /^\s*#/);

        @parts = split(/\t/, $line);
        if (scalar(@parts) > 3) {
            die("Unable to parse population map, '$popmap_path' (map should contain no more than three columns).\n");
        }

        push(@{$sample_list}, $parts[0]);
        $pop_ids->{$parts[0]} = $parts[1];
        $pops->{$parts[1]}++;
        if (scalar(@parts) > 2) {
            $grp_ids->{$parts[0]} = $parts[2];
            $grps->{$parts[2]}++;
        }
    }

    if (scalar(keys %{$grps}) == 0) {
    $grps->{"1"}++;

        foreach $sample (@{$sample_list}) {
            $grp_ids->{$sample} = "1";
        }
    }

    print STDERR "Parsed population map: ", scalar(@{$sample_list}), " files in ", scalar(keys %{$pops});
    scalar(keys %{$pops}) == 1 ?  print STDERR " population" : print STDERR " populations";
    print STDERR " and ", scalar(keys %{$grps});
    scalar(keys %{$grps}) == 1 ? print STDERR " group.\n" : print STDERR " groups.\n";

    close($fh);
}

sub initialize_samples {
    my ($parents, $progeny, $samples, $sample_list, $pop_ids, $grp_ids) = @_;

    if (scalar(@{$sample_list}) > 0 && scalar(@{$samples}) == 0) {
        my @suffixes = ("fq",    "fastq", "fq.gz",   "fastq.gz", "fa",    "fasta", "fa.gz",   "fasta.gz");
        my @fmts     = ("fastq", "fastq", "gzfastq", "gzfastq",  "fasta", "fasta", "gzfasta", "gzfasta");

        #
        # If a population map was specified and no samples were provided on the command line.
        #
        my ($i, $extension, $extension_pe);
        my $first = true;
        foreach $sample (@{$sample_list}) {
            if ($first) {
                $first = false;
                my $found = false;
                for ($i = 0; $i < scalar(@suffixes); $i++) {
                    if (!$paired && -e $sample_path . $sample . "." . $suffixes[$i]) {
                        $found = true;
                        $extension = "." . $suffixes[$i];
                        last;
                    } elsif (-e $sample_path . $sample . ".1." . $suffixes[$i]) {
                        $found = true;
                        $extension = ".1." . $suffixes[$i];
                        $extension_pe = ".2." . $suffixes[$i];
                        last;
                    }
                }
                if (!$found) {
                    print STDERR "Error: Failed to find the first reads file '$sample_path$sample(.1).(fq|fastq|fa|fasta)(.gz)'.\n";
                    exit 1;
                }
                if ($i == 2 || $i == 3 || $i == 6 || $i == 7) {
                    $gzip = true;
                }
            }

            my ($path, $path_pe);

            $path    = $sample_path . $sample . $extension;
            $path_pe = ($paired ? $sample_path . $sample . $extension_pe : "");
            
            die("Error: Failed to open single-end file '$path'.\n") if (! -e $path);
            die("Error: Failed to open paired-end file '$path_pe'.\n") if ($paired && ! -e $path_pe);
            die("Unable to find an entry for '" . $sample . "' in the population map, '$popmap_path'.\n") if (!defined($pop_ids->{$sample}));

            if ($pop_ids->{$sample} eq "parent") {
                push(@{$parents}, {'path'   => $path,
                                       'file'   => $sample,
                                       'suffix' => $suffixes[$i],
                                       'type'   => "parent",
                                       'fmt'    => $fmts[$i],
                                       'path_pe' => $path_pe});
            } elsif ($pop_ids->{$sample} eq "progeny") {
                push(@{$progeny}, {'path'   => $path,
                                       'file'   => $sample,
                                       'suffix' => $suffixes[$i],
                                       'type'   => "progeny",
                                       'fmt'    => $fmts[$i],
                                       'path_pe' => $path_pe});
            } else {
                push(@{$samples}, {'path'   => $path,
                                       'file'   => $sample,
                                       'suffix' => $suffixes[$i],
                                       'type'   => "sample",
                                       'fmt'    => $fmts[$i],
                                       'path_pe' => $path_pe});
            }
        }
    }

    #
    # If a population map was specified, make sure all samples in the list were found (and vice versa) and assign popualtion IDs.
    #
    my %sample_hash;

    foreach $sample (@{$samples}) {
        $sample_hash{$sample->{'file'}}++;

        if (!defined($pop_ids->{$sample->{'file'}})) {
            die("Unable to find an entry for '" . $sample->{'file'} . "' in the population map, '$popmap_path'.\n");
        } else {
            $sample->{'pop_id'} = $pop_ids->{$sample->{'file'}};
        }
        if (!defined($grp_ids->{$sample->{'file'}})) {
            die("Unable to find an entry for '" . $sample->{'file'} . "' in the population map, '$popmap_path'.\n");
        } else {
            $sample->{'grp_id'} = $grp_ids->{$sample->{'file'}};
        }
    }

    foreach $sample (@{$sample_list}) {
        if (!defined($sample_hash{$sample})) {
            die("Unable to find a file corresponding to the population map entry '" . $sample . "' in the population map, '$popmap_path'.\n");
        }
    }

    #
    # Check that no duplicate files were specified.
    #
    my (%files, $file);
    foreach $file (@{$parents}, @{$progeny}, @{$samples}) {
        $files{$file}++;
    }
    foreach $file (keys %files) {
        if ($files{$file} > 1) {
            die("A duplicate file was specified which may create undefined results, '$file'\n");
        }
    }

    print STDERR "Found ", scalar(@{$parents}), " parental file(s).\n\n" if (scalar(@{$parents}) > 0);
    print STDERR "Found ", scalar(@{$progeny}), " progeny file(s).\n\n" if (scalar(@{$progeny}) > 0);
    print STDERR "Found ", scalar(@{$samples}), " sample file(s).\n\n" if (scalar(@{$samples}) > 0);
}

sub initialize_database {
    my ($log_fh, $parents, $progeny, $samples, $sample_ids) = @_;

    my (@results, $sample_id, $sample);

    print $log_fh "Initializing the database...\n";

    #
    # Create the database.
    #
    if ($create_db) {
        #
        # Check that the database doesn't already exist.
        #
        if ($dry_run == false) {
            @results = `mysql --defaults-file=$cnf -N -B -e "SHOW DATABASES LIKE '$db'"`;
            if (scalar(@results) > 0 && $overw_db == false) {
                die("Unable to create database '$db', it already exists.\n");
            }
        }

        if ($overw_db == true) {
            `mysql --defaults-file=$cnf -N -B -e "DROP DATABASE IF EXISTS $db"` if ($dry_run == false);
            print $log_fh "mysql --defaults-file=$cnf -N -B -e \"DROP DATABASE IF EXISTS $db\"\n";
        }

        `mysql --defaults-file=$cnf -e "CREATE DATABASE $db"` if ($dry_run == false);
        print $log_fh "mysql --defaults-file=$cnf $db -e \"CREATE DATABASE $db\"\n";
        `mysql --defaults-file=$cnf $db < $mysql_tables` if ($dry_run == false);
        print $log_fh "mysql --defaults-file=$cnf $db < $mysql_tables\n";
    }

    #
    # Set the SQL Batch ID for this set of loci, along with description and date of
    # sequencing. Insert this batch data into the database.
    #
    `mysql --defaults-file=$cnf $db -e "INSERT INTO batches SET id=$batch_id, description='$desc', date='$date', type='$data_type'"` if ($dry_run == false);

    print $log_fh "mysql --defaults-file=$cnf $db -e \"INSERT INTO batches SET id=$batch_id, description='$desc', date='$date', type='$data_type'\"\n";

    print $log_fh "Loading sample data into the MySQL database...\n";

    my $i = 1;

    foreach $sample (@{$parents}, @{$progeny}, @{$samples}) {

    if ($dry_run == false) {
        `mysql --defaults-file=$cnf $db -e "INSERT INTO samples SET sample_id=$i, batch_id=$batch_id, type='$sample->{'type'}', file='$sample->{'file'}', pop_id='$sample->{'pop_id'}', group_id='$sample->{'grp_id'}'"`;
        @results = `mysql --defaults-file=$cnf $db -N -B -e "SELECT id FROM samples WHERE sample_id=$i AND batch_id=$batch_id AND type='$sample->{'type'}' AND file='$sample->{'file'}'"`;
        chomp $results[0];
        $sample_id = $results[0];

        #
        # Save the sample ID to use when running ustacks.
        #
        $sample_ids->{$sample->{'file'}} = $sample_id;
    }
    print $log_fh "mysql --defaults-file=$cnf $db -e \"INSERT INTO samples SET sample_id=$i, batch_id=$batch_id, type='$sample->{'type'}', file='$sample->{'file'}', pop_id='$sample->{'pop_id'}', group_id='$sample->{'grp_id'}'\"\n";

        $i++;
    }

    print $log_fh "\n";
}

sub load_sql_data {
    my ($log_fh, $pops, $parents, $progeny, $samples) = @_;

    my ($pop_cnt, $sample, $num_files, $i, $file);

    print STDERR "\nComputation is complete, loading results to the database '$db'.\n";

    my $pop_cnt = scalar(keys %{$pops});

    $i         = 1;
    $num_files = scalar(@{$parents}) + scalar(@{$progeny}) + scalar(@{$samples});

    foreach $sample (@{$parents}, @{$progeny}, @{$samples}) {

        printf(STDERR "Loading ustacks output to $db; file % 3s of % 3s [%s]...", $i, $num_files, $sample->{'file'});

        if ($gzip == true) {
            $file = "$out_path/$sample->{'file'}" . ".tags.tsv.gz";
            import_gzsql_file($log_fh, $file, "unique_tags", 1);

            $file = "$out_path/$sample->{'file'}" . ".snps.tsv.gz";
            import_gzsql_file($log_fh, $file, "snps", 1);

            $file = "$out_path/$sample->{'file'}" . ".alleles.tsv.gz";
            import_gzsql_file($log_fh, $file, "alleles", 1);

        } else {
            $file = "$out_path/$sample->{'file'}" . ".tags.tsv";
            import_sql_file($log_fh, $file, "unique_tags", 1);

            $file = "$out_path/$sample->{'file'}" . ".snps.tsv";
            import_sql_file($log_fh, $file, "snps", 1);

            $file = "$out_path/$sample->{'file'}" . ".alleles.tsv";
            import_sql_file($log_fh, $file, "alleles", 1);
        }
        print STDERR "done.\n";

        $i++;
    }

    print STDERR "Importing catalog to $db...";

    my $cat_file = "batch_" . $batch_id;

    if ($gzip == true) {
        $file = "$out_path/$cat_file" . ".catalog.tags.tsv.gz";
        import_gzsql_file($log_fh, $file, "catalog_tags", 1);

        $file = "$out_path/$cat_file" . ".catalog.snps.tsv.gz";
        import_gzsql_file($log_fh, $file, "catalog_snps", 1);

        $file = "$out_path/$cat_file" . ".catalog.alleles.tsv.gz";
        import_gzsql_file($log_fh, $file, "catalog_alleles", 1);

    } else {
        $file = "$out_path/$cat_file" . ".catalog.tags.tsv";
        import_sql_file($log_fh, $file, "catalog_tags", 1);

        $file = "$out_path/$cat_file" . ".catalog.snps.tsv";
        import_sql_file($log_fh, $file, "catalog_snps", 1);

        $file = "$out_path/$cat_file" . ".catalog.alleles.tsv";
        import_sql_file($log_fh, $file, "catalog_alleles", 1);
    }
    print STDERR "done.\n";

    #
    # Load the sstacks results to the database if requested.
    #
    $i = 1;
    foreach $sample (@{$parents}, @{$progeny}, @{$samples}) {
        printf(STDERR "Loading sstacks output to $db; file % 3s of % 3s [%s]...", $i, $num_files, $sample->{'file'});
        if ($gzip == true) {
            $file = "$out_path/" . $sample->{'file'} . ".matches.tsv.gz";
            import_gzsql_file($log_fh, $file, "matches", 1);
        } else {
            $file = "$out_path/" . $sample->{'file'} . ".matches.tsv";
            import_sql_file($log_fh, $file, "matches", 1);
        }
        print STDERR "done.\n";
        $i++;
    }

    if ($data_type eq "map") {
        #
        # Load the outputs generated by genotypes to the database.
        #
        $file = "$out_path/batch_" . $batch_id . ".markers.tsv";
        import_sql_file($log_fh, $file, "markers", 1);

        $file = "$out_path/batch_" . $batch_id . ".genotypes_1.txt";
        import_sql_file($log_fh, $file, "catalog_genotypes", 1);

    } else {
        #
        # Load the outputs generated by populations to the database.
        #
        $file = "$out_path/batch_" . $batch_id . ".markers.tsv";
        import_sql_file($log_fh, $file, "markers", 1);

        $file = "$out_path/batch_" . $batch_id . ".sumstats.tsv";
        import_sql_file($log_fh, $file, "sumstats", $pop_cnt+1);

        $file = "$out_path/batch_" . $batch_id . ".hapstats.tsv";
        import_sql_file($log_fh, $file, "hapstats", $pop_cnt+1);

        #
        # Import the Fst files.
        #
        my $fst_cnt = 0;
        my (@keys, $m, $n);
        @keys = sort keys %pops;
        for ($m = 0; $m < scalar(@keys); $m++) {
            for ($n = 0; $n < scalar(@keys); $n++) {
                $file = "$out_path/batch_" . $batch_id . ".fst_" . $keys[$m] . "-" . $keys[$n] . ".tsv";

                if (-e $file) {
                    import_sql_file($log_fh, $file, "fst", 1);
                    $fst_cnt++;
                }
            }
        }
        print STDERR "Imported $fst_cnt Fst file(s).\n";

        #
        # Import the Phi_st files.
        #
        $fst_cnt = 0;
        for ($m = 0; $m < scalar(@keys); $m++) {
            for ($n = 0; $n < scalar(@keys); $n++) {
                $file = "$out_path/batch_" . $batch_id . ".phistats_" . $keys[$m] . "-" . $keys[$n] . ".tsv";

                if (-e $file) {
                    import_sql_file($log_fh, $file, "phist", 3);
                    $fst_cnt++;
                }
            }
        }
        print STDERR "Imported $fst_cnt Haplotype Fst file(s).\n";
    }

    print $log_fh "\n";

    #
    # Index the radtags database
    #
    my ($cmd, @results);
    print STDERR "Indexing the database...\n";
    $cmd = $exe_path . "index_radtags.pl -D $db -t -c 2>&1";
    print STDERR  "$cmd\n";
    print $log_fh "$cmd\n";
    @results =    `$time $cmd` if ($dry_run == false);
    print $log_fh @results;
}

sub write_results {
    my ($results, $log_fh) = @_;

    my $line;

    foreach $line (@{$results}) {
        if ($line =~ /\r/) {
            $line =~ s/^.+\r(.*\n)$/\1/;
        }
        print $log_fh $line;
    }
}

sub write_depths_of_cov {
    my ($depths, $log_fh) = @_;

    print STDERR "\nDepths of Coverage for Processed Samples:\n";
    print $log_fh "\nDepths of Coverage for Processed Samples:\n";

    foreach $a (@{$depths}) {
        print STDERR  $a->[0], ": ", $a->[1], "x\n";
        print $log_fh $a->[0], ": ", $a->[1], "x\n";
    }
}

sub import_sql_file {
    my ($log_fh, $file, $table, $skip_lines) = @_;

    my (@results, $ignore);

    $ignore = "IGNORE $skip_lines LINES" if ($skip_lines > 0);

    @results = `mysql --defaults-file=$cnf $db -e "LOAD DATA LOCAL INFILE '$file' INTO TABLE $table $ignore"` if ($sql == true && $dry_run == false);

    if ($sql == true) {
        print $log_fh "mysql --defaults-file=$cnf $db -e \"LOAD DATA LOCAL INFILE '$file' INTO TABLE $table $ignore\"\n", @results;
    }
}

sub import_gzsql_file {
    my ($log_fh, $file, $table, $skip_lines) = @_;

    my (@results, $ignore);

    $ignore = "IGNORE $skip_lines LINES" if ($skip_lines > 0);

    #
    # Get a temporary file name and create a named pipe.
    #
    my $tmpdir     = File::Spec->tmpdir();
    my $named_pipe = mktemp($tmpdir . "/denovo_map_XXXXXX");
    if ($sql == true && $dry_run == false) {
        mkfifo($named_pipe, 0700) || die("Unable to create named pipe for loading gzipped data: $named_pipe, $!");
        print $log_fh "Streaming $file into named pipe $named_pipe.\n";
    }

    #
    # Dump our gzipped data onto the named pipe.
    #
    system("gunzip -c $file > $named_pipe &") if ($sql == true && $dry_run == false);

    @results = `mysql --defaults-file=$cnf $db -e "LOAD DATA LOCAL INFILE '$named_pipe' INTO TABLE $table $ignore"` if ($sql == true && $dry_run == false);

    if ($sql == true) {
      print $log_fh "mysql --defaults-file=$cnf $db -e \"LOAD DATA LOCAL INFILE '$named_pipe' INTO TABLE $table $ignore\"\n", @results;
    }

    #
    # Remove the pipe.
    #
    unlink($named_pipe) if ($sql == true && $dry_run == false);
}

sub parse_command_line {
    my ($arg);

    my $ustacks_mismatch = -1;
    my $cstacks_mismatch = -1;

    while (@ARGV) {
        $_ = shift @ARGV;
        if    ($_ =~ /^-v$/) { version(); exit 1; }
        elsif ($_ =~ /^-h$/) { usage(); }
        elsif ($_ =~ /^-d$/ || $_ =~ /^--dry-run$/) { $dry_run   = true; }
        elsif ($_ =~ /^-o$/) { $out_path  = shift @ARGV; }
        elsif ($_ =~ /^-e$/) { $exe_path  = shift @ARGV; }
        # elsif ($_ =~ /^-D$/) { $desc      = shift @ARGV; }
        # elsif ($_ =~ /^-b$/) { $batch_id  = shift @ARGV; }
        # elsif ($_ =~ /^-i$/) { $sample_id = shift @ARGV; }
        # elsif ($_ =~ /^-a$/) { $date      = shift @ARGV; }
        # elsif ($_ =~ /^-S$/) { $sql       = false; }
        # elsif ($_ =~ /^-B$/) { $db        = shift @ARGV; }
        elsif ($_ =~ /^-m$/)        { $min_cov     = shift @ARGV; }
        elsif ($_ =~ /^-P$/)        { $min_rcov    = shift @ARGV; }
        elsif ($_ =~ /^--paired$/)  { $paired      = true; } 
        elsif ($_ =~ /^--samples$/) { $sample_path = shift @ARGV; } 
        elsif ($_ =~ /^-O$/ || $_ =~ /^--popmap$/) {
            $popmap_path = shift @ARGV;
            push(@_cstacks,     "-M " . $popmap_path);
            push(@_sstacks,     "-M " . $popmap_path);
            push(@_tsv2bam,     "-M " . $popmap_path);
            push(@_gstacks,     "-M " . $popmap_path);
            push(@_populations, "-M " . $popmap_path);

        } elsif ($_ =~ /^--gapped$/) {
            $gapped_alns = true;
            push(@_ustacks, "--gapped ");
            push(@_cstacks, "--gapped ");
            push(@_sstacks, "--gapped ");

        # } elsif ($_ =~ /^--create_db$/) {
        #         $create_db = true;

        #     } elsif ($_ =~ /^--overw_db$/) {
        #         $overw_db  = true;
        #         $create_db = true;

        #     } elsif ($_ =~ /^-A$/) {
        #     $arg = shift @ARGV;
        #     push(@_genotypes, "-t " . $arg);

        #     $arg = lc($arg);
        #     if ($arg ne "gen" && $arg ne "cp" && $arg ne "f2" && $arg ne "bc1" && $arg ne "dh") {
        # 	print STDERR "Unknown genetic mapping cross specified: '$arg'\n";
        # 	usage();
        #     }

        } elsif ($_ =~ /^-t$/) {
            push(@_ustacks, "-d ");

        } elsif ($_ =~ /^-T$/) {
            $arg = shift @ARGV;
            push(@_ustacks, "-p " . $arg);
            push(@_cstacks, "-p " . $arg);
            push(@_sstacks, "-p " . $arg);
            push(@_tsv2bam, "-t " . $arg);
            push(@_gstacks, "-t " . $arg);
            push(@_populations, "-t " . $arg);

        } elsif ($_ =~ /^-M$/) {
            $ustacks_mismatch = shift(@ARGV);
            push(@_ustacks,   "-M " . $ustacks_mismatch);

        } elsif ($_ =~ /^-N$/) {
            push(@_ustacks,   "-N " . shift @ARGV);

        } elsif ($_ =~ /^-n$/) {
            $cstacks_mismatch = shift @ARGV;
            push(@_cstacks, "-n " . $cstacks_mismatch);

        } elsif ($_ =~ /^--bound_low$/) {
            push(@_ustacks, "--bound_low " . shift @ARGV);
            push(@_ustacks, "--model_type bounded");

        } elsif ($_ =~ /^--bound_high$/) {
            push(@_ustacks, "--bound_high " . shift @ARGV);
            push(@_ustacks, "--model_type bounded");

        } elsif ($_ =~ /^--alpha$/) {
            push(@_ustacks, "--alpha " . shift @ARGV);

        } elsif ($_ =~ /^-X$/) {
            #
            # Pass an arbitrary command-line option to a pipeline program.
            #
            # Command line option must be of the form '-X "program:option"'
            #
            $arg = shift @ARGV;
            my ($prog, $opt) = ($arg =~ /^(\w+):(.+)$/);
            if ($prog eq "ustacks") {
            	push(@_ustacks, $opt);
            } elsif ($prog eq "cstacks") {
            	push(@_cstacks, $opt);
            } elsif ($prog eq "sstacks") {
            	push(@_sstacks, $opt);
            } elsif ($prog eq "tsv2bam") {
            	push(@_tsv2bam, $opt);
            } elsif ($prog eq "gstacks") {
            	push(@_gstacks, $opt);
            } elsif ($prog eq "genotypes") {
            	push(@_genotypes, $opt);
            } elsif ($prog eq "populations") {
            	push(@_populations, $opt);
            } else {
            	print STDERR "Unknown pipeline program, '$arg'\n";
            	usage();
            }
        } elsif ($_ =~ /^--time-components$/) {
            $time = '/usr/bin/time';
            if (! -e $time) {
                die "Error: '$time': No such file or directory.\n";
            }
        } else {
            print STDERR "Unknown command line option: '$_'\n";
            usage();
        }
    }

    $exe_path = $exe_path . "/"          if (substr($exe_path, -1) ne "/");
    $out_path = substr($out_path, 0, -1) if (substr($out_path, -1) eq "/");

    if ($sql == true && length($date) == 0) {
        $date = strftime("%Y-%m-%d", (localtime(time)));
    }

    if (scalar(@parents) > 0 && scalar(@samples) > 0) {
        print STDERR "You must specify either parent or sample files, but not both.\n";
        usage();
    }

    if (length($popmap_path) == 0) {
        print STDERR "You must specify a population map that lists your sample names (--popmap).\n";
        usage();
    }

    if (length($sample_path) == 0) {
        print STDERR "You must specify the path to the directory containing the samples (--samples).\n";
        usage();
    }

    if (length($sample_path) > 0) {
        $sample_path .= "/" if (substr($sample_path, -1) ne "/");
    }

    if (scalar(@samples) > 0 || length($popmap_path) > 0) {
        $data_type = "population";
    } else {
        $data_type = "map";
    }

    if ($paired == true) {
        if (length($sample_path) == 0) {
            print STDERR "If you want to assemble paired-ends, you must use a population map to specify samples.\n";
            usage();
        } else {
            push(@_tsv2bam, "-R $sample_path");
        }
    }

    #
    # By default, we want ustacks -M to equal cstacks -n.
    #
    if ($cstacks_mismatch == -1 && $ustacks_mismatch > 0) {
        push(@_cstacks, "-n " . $ustacks_mismatch);
    }
}

sub version {
    print STDERR "denovo_map.pl ", stacks_version, "\n";
}

sub usage {
    version();

    print STDERR <<EOQ;
denovo_map.pl --samples dir --popmap path -o dir [--paired] (assembly options) (database options) [-X prog:"opts" ...]

  Input/Output files:
    --samples: path to the directory containing the samples reads files.
    --popmap: path to a population map file (format is "<name> TAB <pop>", one sample per line).
    o: path to an output directory.

  General options:
    A: for a mapping cross, specify the type; one of 'CP', 'F2', 'BC1', 'DH', or 'GEN'.
    X: additional options for specific pipeline components, e.g. -X "populations: -p 3 -r 0.50".
    T: the number of threads/CPUs to use (default: 1).
    d: Dry run. Do not actually execute anything, just print the commands that would be executed.

  Stack assembly options:
    M: number of mismatches allowed between stacks within individuals (for ustacks).
    n: number of mismatches allowed between stacks between individuals (for cstacks).
    --gapped: perform gapped comparisons (for ustacks, cstacks, sstacks; default: off).

  SNP model options:
    --alpha: significance level at which to call genotypes (for ustacks; default: 0.05).
    For the bounded model:
      --bound_low <num>: lower bound for epsilon, the error rate, between 0 and 1.0 (default 0.0).
      --bound_high <num>: upper bound for epsilon, the error rate, between 0 and 1.0 (default 1.0).

  Paired-end options:
    --paired: after assembling RAD loci, assemble mini-contigs with paired-end reads.

  Database options (currently disabled in beta):
    S: disable database interaction.
    B: specify an SQL database to load data into.
    D: a description of this batch to be stored in the database.
    i: starting sample_id, this is determined automatically if database interaction is enabled.
    --create_db: create the database specified by '-B' and populate the tables.
    --overw_db: delete the database before creating a new copy of it (turns on --create_db).

  Miscellaneous:
    --time-components (for benchmarking)
EOQ

    exit 1;
}
