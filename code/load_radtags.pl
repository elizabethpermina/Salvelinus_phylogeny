#!/usr/bin/env perl
#
# Copyright 2011-2014, Julian Catchen <jcatchen@uoregon.edu>
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
# Load a set of output files from the Stacks pipeline into a Stacks MySQL database.
#
# By Julian Catchen <jcatchen@uoregon.edu>
#

use strict;
use POSIX;
use File::Temp qw/ mktemp /;
use File::Spec;
use constant stacks_version => "2.0Beta7c";

my $mysql_config = "/home/epermina/soft/stacks-2.0Beta7c/share/stacks/" . "sql/mysql.cnf";
my $dry_run      = 0;
my $db           = "";
my $in_path      = ".";
my $sample_id    = 0;
my $desc         = "";
my $date         = "";
my $batch_id     = 0;
my $batch        = 0;
my $catalog      = 0;
my $stacks_type  = "";
my $popmap_path  = "";
my $ignore_tags  = 0;
my $white_list   = "";

parse_command_line();

my $cnf = (-e $ENV{"HOME"} . "/.my.cnf") ? $ENV{"HOME"} . "/.my.cnf" : $mysql_config;

if (length($date) == 0) {
    $date = strftime("%Y-%m-%d", (localtime(time)));
}

my (@results, @files, @catalog, @parent_ids, @pop_ids, %pops, %sample_ids);

build_file_list(\@files, \@catalog);

extract_parental_ids(scalar(@files), \@catalog, \@parent_ids);

extract_sample_ids(\@files, \%sample_ids);

parse_population_map(\@files, \%sample_ids, \@pop_ids, \%pops);

print STDERR 
    "Stacks pipeline type: '", $stacks_type, "'\n",
    scalar(@files), " files to process: ", join(", ", @files), "\n",
    scalar(@catalog), " catalog files to process: ", join(", ", @catalog), "\n";
if ($stacks_type eq "map") {
    print STDERR
	scalar(@parent_ids), " parent IDs identified: ", join(", ", @parent_ids), "\n";
}

if ($batch) {
    if (!$dry_run) {
	@results = `mysql --defaults-file=$cnf $db -e "INSERT INTO batches SET id=$batch_id, description='$desc', date='$date', type='$stacks_type'"`;
    }
    print STDERR
	"mysql --defaults-file=$cnf $db ",
	"-e \"INSERT INTO batches SET id=$batch_id, description='$desc', date='$date', type='$stacks_type'\"\n",
	@results;
}

my ($file, $f, $i, $cnt, $type, $pop_id);

#
# Import the catalog
#
if ($catalog) {
    foreach $file (@catalog) {

	$f = $in_path . "/$file" . ".catalog.tags.tsv";
	if (-e $f) {
	    import_sql_file($f, "catalog_tags", 1);
	} elsif (-e $f . ".gz") {
	    $f = $in_path . "/$file" . ".catalog.tags.tsv.gz";
	    import_gzsql_file($f, "catalog_tags", 1);
	}

        $f = $in_path . "/$file" . ".catalog.snps.tsv";
	if (-e $f) {
	    import_sql_file($f, "catalog_snps", 1);
	} elsif (-e $f . ".gz") {
	    $f = $in_path . "/$file" . ".catalog.snps.tsv.gz";
	    import_gzsql_file($f, "catalog_snps", 1);
	}

        $f = $in_path . "/$file" . ".catalog.alleles.tsv";
	if (-e $f) {
	    import_sql_file($f, "catalog_alleles", 1);
	} elsif (-e $f . ".gz") {
	    $f = $in_path . "/$file" . ".catalog.alleles.tsv.gz";
	    import_gzsql_file($f, "catalog_alleles", 1);
	}
    }
}

if ($stacks_type eq "map") {
    $f = "$in_path/batch_" . $batch_id . ".markers.tsv";
    import_sql_file($f, "markers", 1);

    $f = "$in_path/batch_" . $batch_id . ".genotypes_1.txt";
    import_sql_file($f, "catalog_genotypes", 1);

} elsif ($stacks_type eq "population") {
    $f = "$in_path/batch_" . $batch_id . ".markers.tsv";
    import_sql_file($f, "markers", 1);

    $f = "$in_path/batch_" . $batch_id . ".sumstats.tsv";
    import_sql_file($f, "sumstats", scalar(keys %pops) + 1);

    $f = "$in_path/batch_" . $batch_id . ".hapstats.tsv";
    import_sql_file($f, "hapstats", scalar(keys %pops) + 1);

    #
    # Import the Fst files.
    #
    my $fst_cnt = 0;
    my (@keys, $m, $n);
    @keys = sort keys %pops;
    for ($m = 0; $m < scalar(@keys); $m++) {
	for ($n = 0; $n < scalar(@keys); $n++) {
	    $f = "$in_path/batch_" . $batch_id . ".fst_" . $keys[$m] . "-" . $keys[$n] . ".tsv";

	    if (-e $file) {
		import_sql_file($f, "fst", 1);
		$fst_cnt++;
	    }
	}
    }
    print STDERR "Imported $fst_cnt SNP Fst file(s).\n";

    #
    # Import the Phi_st files.
    #
    $fst_cnt = 0;
    for ($m = 0; $m < scalar(@keys); $m++) {
	for ($n = 0; $n < scalar(@keys); $n++) {
	    $f = "$in_path/batch_" . $batch_id . ".phistats_" . $keys[$m] . "-" . $keys[$n] . ".tsv";

	    if (-e $file) {
		import_sql_file($f, "phist", 3);
		$fst_cnt++;
	    }
	}
    }
    print STDERR "Imported $fst_cnt Haplotype Fst file(s).\n";
}

$i = 1;
$cnt = scalar(@files);

foreach $file (sort {$sample_ids{$a} <=> $sample_ids{$b}} @files) {
    print STDERR "Processing sample $i of $cnt\n";

    $f = $in_path . "/$file" . ".matches.tsv";
    if (-e $f) {
	import_sql_file($f, "matches", 1);
    } elsif (-e $f . ".gz") {
	$f = $in_path . "/$file" . ".matches.tsv.gz";
	import_gzsql_file($f, "matches", 1);
    }
    $i++;
}

$i = 1;
foreach $file (sort {$sample_ids{$a} <=> $sample_ids{$b}} @files) {
    print STDERR "Processing sample $i of $cnt\n";

    #
    # Pull out the sample ID and insert it into the database
    #
    $sample_id = $sample_ids{$file};
    if ($stacks_type eq "map") {
	$type = (grep(/^$sample_id$/, @parent_ids) > 0) ? 'parent' : 'progeny';
    } else {
	$type = "sample";
    }

    $pop_id = shift(@pop_ids);

    if (!$dry_run) {
	@results = `mysql --defaults-file=$cnf $db -e "INSERT INTO samples SET id=$sample_id, sample_id=$sample_id, batch_id=$batch_id, type='$type', file='$file', pop_id='$pop_id'"`;
    }
    print STDERR 
	"mysql --defaults-file=$cnf $db ",
	"-e \"INSERT INTO samples SET id=$sample_id, sample_id=$sample_id, batch_id=$batch_id, type='$type', file='$file', pop_id='$pop_id'\"\n", 
	@results;

    $f = $in_path . "/$file" . ".tags.tsv";
    if (-e $f) {
	import_sql_file($f, "unique_tags", 1) if ($ignore_tags == 0);
    } elsif (-e $f . ".gz") {
	$f = $in_path . "/$file" . ".tags.tsv.gz";
	import_gzsql_file($f, "unique_tags", 1) if ($ignore_tags == 0);
    }

    $f = $in_path . "/$file" . ".snps.tsv";
    if (-e $f) {
	import_sql_file($f, "snps", 1);
    } elsif (-e $f . ".gz") {
	$f = $in_path . "/$file" . ".snps.tsv.gz";
	import_gzsql_file($f, "snps", 1);
    }

    $f = $in_path . "/$file" . ".alleles.tsv";
    if (-e $f) {
	import_sql_file($f, "alleles", 1);
    } elsif (-e $f . ".gz") {
	$f = $in_path . "/$file" . ".alleles.tsv.gz";
	import_gzsql_file($f, "alleles", 1);
    }

    $i++;
}

print STDERR "\nDon't forget to index your Stacks database -- run index_radtags.pl\n\n";

sub parse_population_map {
    my ($samples, $sample_ids, $pop_ids, $pops) = @_;

    my ($fh, @parts, $line, %ids, $file, $path);

    if (length($popmap_path) == 0) {
	foreach $path (@{$samples}) {
	    push(@{$pop_ids}, 1);
	    $pops->{1}++;
	}
	return;
    }

    open($fh, "<$popmap_path") or die("Unable to open population map, '$popmap_path', $!\n");

    while ($line = <$fh>) {
	chomp $line;
	@parts = split(/\t/, $line);

	if (scalar(@parts) > 3) {
	    die("Unable to parse population map, '$popmap_path' (map should contain no more than three columns).\n");
	}

	$ids{$parts[0]} = $parts[1];
    }

    foreach $file (sort {$sample_ids->{$a} <=> $sample_ids->{$b}} @{$samples}) {
	if (!defined($ids{$file})) {
	    die("Unable to find '$file' in the population map, '$popmap_path'.\n");
	}

	push(@{$pop_ids}, $ids{$file});
	$pops->{$ids{$file}}++;
    }

    print STDERR "Parsed population map: ", scalar(@{$samples}), " files in ", scalar(keys %pops), " populations.\n";

    close($fh);
}

sub extract_parental_ids {
    my ($sample_cnt, $catalog, $parental_ids) = @_;

    my ($fh, $prefix, $path, $line, @parts, $tag_id, @tag_ids, $id, $tag, %ids);

    print STDERR "Scanning catalog for sample IDs...";

    foreach $prefix (@catalog) {
        $path = $in_path . "/" . $prefix . ".catalog.tags.tsv";

	if (-e $path) {
	    open($fh, "<$path") or die("Unable to open catalog file: '$path', $!\n");
	} elsif (-e $path . ".gz") {
	    open($fh, "gunzip -c " . $path . ".gz |") or die("Unable to open catalog file: '$path', $!\n");
	}

	while ($line = <$fh>) {
	    chomp $line;
	    @parts = split(/\t/, $line);

	    @tag_ids = split(/,/, $parts[8]);

	    foreach $tag_id (@tag_ids) {
		($id, $tag) = split(/_/, $tag_id);
		$ids{$id}++;
	    }
	}

	close($fh);
    }

    @{$parental_ids} = keys %ids;

    #
    # Determine the type of pipeline run: either a 'map' or a 'population' type.
    # If all samples are parental, i.e. in the catalog, then this is a population type
    # otherwise, it is a map type.
    #
    if (length($stacks_type) == 0) {
	$stacks_type = (scalar(@{$parental_ids}) == $sample_cnt) ? "population" : "map";
    }

    print STDERR "done.\n";
}

sub extract_sample_ids {
    my ($files, $sample_ids) = @_;

    my ($file, $f, $line, @results, @parts);

    print STDERR "Collecting sample IDs from Stacks output files...";

    foreach $file (@{$files}) {
	$f = $in_path . "/$file" . ".tags.tsv";

	if (-e $f) {
	    @results = `head -n 2 $f | tail -n 1`;

	} elsif (-e $f . ".gz") {
	    $f = $in_path . "/$file" . ".tags.tsv.gz";
	    @results = `gunzip -c $f | head -n 2 | tail -n 1`;

	} else {
	    die("Unable to find file $f\n");
	}

	chomp $results[0];
	@parts = split(/\t/, $results[0]);

	#
	# Sample ID is expected to be the first column in the *.tags.tsv file.
	#
	$sample_ids->{$file} = $parts[1];
    }

    print STDERR "done.\n";
}

sub import_sql_file {
    my ($file, $table, $skip_lines) = @_;

    my (@results, $ignore);

    if (!-e $file) {
	print STDERR "File '$file' does not exist.\n";
	return;
    }

    $ignore = " IGNORE $skip_lines LINES" if ($skip_lines > 0);

    if (!$dry_run) {
	@results = `mysql --defaults-file=$cnf $db -e "LOAD DATA LOCAL INFILE '$file' INTO TABLE $table$ignore"`;
    }
    print STDERR 
	"mysql --defaults-file=$cnf $db ",
	"-e \"LOAD DATA LOCAL INFILE '$file' INTO TABLE $table$ignore\"\n", 
	@results;
}

sub import_gzsql_file {
    my ($file, $table, $skip_lines) = @_;

    my (@results, $ignore);

    if (!-e $file) {
	print STDERR "File '$file' does not exist.\n";
	return;
    }

    $ignore = "IGNORE $skip_lines LINES" if ($skip_lines > 0);

    #
    # Get a temporary file name and create a named pipe.
    #
    my $tmpdir     = File::Spec->tmpdir();
    my $named_pipe = mktemp($tmpdir . "/denovo_map_XXXXXX");
    if ($dry_run == 0) {
	mkfifo($named_pipe, 0700) || die("Unable to create named pipe for loading gzipped data: $named_pipe, $!");
    }
    print STDERR "Streaming $file into named pipe $named_pipe.\n";

    #
    # Dump our gzipped data onto the named pipe.
    #
    system("gunzip -c $file > $named_pipe &") if ($dry_run == 0);

    @results = `mysql --defaults-file=$cnf $db -e "LOAD DATA LOCAL INFILE '$named_pipe' INTO TABLE $table $ignore"` if ($dry_run == 0);

    print STDERR "mysql --defaults-file=$cnf $db -e \"LOAD DATA LOCAL INFILE '$named_pipe' INTO TABLE $table $ignore\"\n", @results;

    #
    # Remove the pipe.
    #
    unlink($named_pipe) if ($dry_run == 0);
}

sub build_file_list {
    my ($files, $catalog_files) = @_;

    my (@wl, @ls, $line, $prefix);

    # Load a white list of files to process if it is supplied.
    if (length($white_list) > 0) {
        load_white_list(\@wl);
    }

    @ls = `ls -1 $in_path/*.tags.tsv* 2> /dev/null`;

    if (scalar(@ls) == 0) {
	print STDERR "Unable to locate any input files to process within '$in_path'\n";
	usage();
    }

    foreach $line (@ls) {
	chomp $line;

	if ($line =~ /\.tags\.tsv\.gz$/) {
	    ($prefix) = ($line =~ /$in_path\/(.+)\.tags\.tsv\.gz/);
	} else {
	    ($prefix) = ($line =~ /$in_path\/(.+)\.tags\.tsv/);
	}

        next if ($prefix =~ /catalog/);

        if (scalar(@wl) > 0) {
            next if (!grep(/^$prefix$/, @wl));
        }

	push(@{$files}, $prefix);
    }

    if ($catalog > 0) {
        @ls = `ls -1 $in_path/*.catalog.tags.tsv* 2> /dev/null`;

        if (scalar(@ls) == 0) {
            print STDERR "Unable to locate any catalog input files to process within '$in_path'\n";
            usage();
        }

        foreach $line (@ls) {
            chomp $line;

	    if ($line =~ /\.catalog\.tags\.tsv\.gz$/) {
		($prefix) = ($line =~ /$in_path\/(.+)\.catalog\.tags\.tsv\.gz/);
	    } else {
		($prefix) = ($line =~ /$in_path\/(.+)\.catalog\.tags\.tsv/);
	    }

            if (scalar(@wl) > 0) {
                next if (!grep(/^$prefix$/, @wl));
            }

            push(@{$catalog_files}, $prefix);
        }
    }
}

sub load_white_list {
    my ($wl) = @_;

    open(WHITE, "<" . $white_list) 
        or die("Unable to open white list file '$white_list': $!\n");

    my $line   = "";

    while ($line = <WHITE>) {
        chomp $line;

        next if (length($line) == 0);

        push(@{$wl}, $line);
    }

    close(WHITE);
}

sub parse_command_line {
    while (@ARGV) {
	$_ = shift @ARGV;
	if    ($_ =~ /^-p$/) { $in_path     = shift @ARGV; }
	elsif ($_ =~ /^-D$/) { $db          = shift @ARGV; }
	elsif ($_ =~ /^-c$/) { $catalog++; }
        elsif ($_ =~ /^-B$/) { $batch++; }
	elsif ($_ =~ /^-b$/) { $batch_id    = shift @ARGV; }
	elsif ($_ =~ /^-e$/) { $desc        = shift @ARGV; }
	elsif ($_ =~ /^-a$/) { $date        = shift @ARGV; }
        elsif ($_ =~ /^-W$/) { $white_list  = shift @ARGV; }
        elsif ($_ =~ /^-t$/) { $stacks_type = shift @ARGV; }
        elsif ($_ =~ /^-M$/) { $popmap_path = shift @ARGV; }
        elsif ($_ =~ /^-U$/) { $ignore_tags++; }
	elsif ($_ =~ /^-d$/) { $dry_run++; }
	elsif ($_ =~ /^-v$/) { version(); exit(); }
	elsif ($_ =~ /^-h$/) { usage(); }
	else {
	    print STDERR "Unknown command line option: '$_'\n";
	    usage();
	}
    }

    $in_path  = substr($in_path, 0, -1)  if (substr($in_path, -1)  eq "/");

    if (length($db) == 0) {
	print STDERR "You must specify a database.\n";
	usage();
    }
}

sub version {
    print STDERR "load_radtags.pl ", stacks_version, "\n";
}

sub usage {
    version();

    print STDERR <<EOQ; 
load_radtags.pl -D db -p path -b batch_id [-B -e desc] [-c] [-M pop_map] [-d] [-t] [-W path] [-U] [-d] [-h]
    D: Database to load data into.
    p: Path to input files.
    b: Batch ID.
    M: if you have analyzed several populations, specify a population map.
    c: Load the catalog into the database.
    B: Load information into batch table.
    e: batch dEscription.
    d: perform a dry run. Do not actually load any data, just print what would be executed.
    W: only load file found on this white list.
    U: do not load stacks to unique_tags table to save database space.
    t: pipeline type (either 'map' or 'population'), load_radtags.pl will guess based on the number or indiviuduals used to build the catalog.
    h: display this help message.

EOQ

exit(0);
}
