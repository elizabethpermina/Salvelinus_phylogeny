#!/usr/bin/env perl
#
# Copyright 2010, Julian Catchen <jcatchen@uoregon.edu>
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
# Written by Julian Catchen <jcatchen@uoregon.edu>
#

use strict;
use File::Temp qw/tempfile/;
use DBI;

use constant stacks_version => "2.0Beta7c";

my $mysql_config  = "/home/epermina/soft/stacks-2.0Beta7c/share/stacks/" . "sql/mysql.cnf";
my $sql_path      = "/home/epermina/soft/stacks-2.0Beta7c/share/stacks/" . "sql/";
my $sql_tag_table = "";
my $sql_cat_table = "";
my $sql_chr_table = "";
my $db            = "";
my $debug         = 0;
my $catalog_index = 0;
my $tag_index     = 0;
my $radome_index  = 0;

parse_command_line();

my $cnf = (-e $ENV{"HOME"} . "/.my.cnf") ? $ENV{"HOME"} . "/.my.cnf" : $mysql_config;

#
# Make sure the SQL definition files are available
#
if ($catalog_index && !-e $sql_cat_table) {
    print STDERR "Unable to locate catalog SQL definition.\n";
    usage();
}

if ($tag_index && !-e $sql_tag_table) {
    print STDERR "Unable to locate tag_index SQL definition.\n";
    usage();
}

#
# Connect to the database and prepare our queries.
#
my (%sth);

prepare_sql_handles(\%sth);

if ($catalog_index) {
    gen_cat_index(\%sth);
}

if ($tag_index) {
    gen_tag_index(\%sth);
}

close_sql_handles(\%sth);

sub gen_cat_index {
    my ($sth) = @_;

    die ("Unable to find SQL file: '$sql_cat_table'\n") if (!-e $sql_tag_table);

    print STDERR "Generating catalog tag index\n";

    my ($fh, $catalog_file) = tempfile("catalog_index_XXXXXXXX", UNLINK => 1, TMPDIR => 1);

    my ($row, $tag, $count, $par_cnt, $pro_cnt, $allele_cnt, $marker, $uncor_marker, $valid_pro, 
        $chisq_pval, $lnl, $ratio, $ests, $pe_radtags, $blast_hits, $geno_cnt, $ref_type, $ref_id, $bp);

    my (%snps, %markers, %genotypes, %seqs, %hits, %parents, %progeny, %alleles, %chrs, %radome);

    print STDERR "  Fetching catalog SNPs...";
    fetch_catalog_snps(\%sth, \%snps);
    print STDERR "done.\n";

    print STDERR "  Fetching markers...";
    fetch_markers(\%sth, \%markers);
    print STDERR "done.\n";

    print STDERR "  Fetching genotypes...";
    fetch_genotypes(\%sth, \%genotypes);
    print STDERR "done.\n";

    print STDERR "  Fetching associated sequences...";
    sequence_matches(\%sth, \%seqs, \%hits);
    print STDERR "done.\n";

    print STDERR "  Fetching catalog matches...";
    catalog_matches(\%sth, \%parents, \%progeny, \%alleles);
    print STDERR "done.\n";

    $ref_type = "";
    $ref_id   = 0;
    if ($radome_index) {
	print STDERR "  Fetching reference RAD data...";
	radome_ref(\%sth, \%radome);
	print STDERR "done.\n";
    }

    print STDERR "  Assembling catalog tags at the database...";
    $sth->{'cat_tags'}->execute()
	or die("Unable to select results from $db.\n");
    print STDERR "done.\n";

    #my $num_rows = $sth->{'cat_tags'}->rows();
    #my $i = 1;

    print STDERR "  Processing catalog tags\n";

    while ($row = $sth->{'cat_tags'}->fetchrow_hashref()) {
	#
	# Determine the number of SNPs contained within this RAD-Tag
	#
	if (defined($snps{$row->{'batch_id'}}->{$row->{'tag_id'}})) {
	    $count = $snps{$row->{'batch_id'}}->{$row->{'tag_id'}};
	} else {
	    $count = 0;
	}

	#
	# Determine the number of parental samples that match this catalog RAD-Tag
	#
	if (defined($parents{$row->{'batch_id'}}->{$row->{'tag_id'}})) {
	    $par_cnt = scalar(keys %{$parents{$row->{'batch_id'}}->{$row->{'tag_id'}}});
	} else {
	    $par_cnt = 0;
	}

	#
	# Determine the number of progeny samples that match this catalog RAD-Tag
	#
	if (defined($progeny{$row->{'batch_id'}}->{$row->{'tag_id'}})) {
	    $pro_cnt = scalar(keys %{$progeny{$row->{'batch_id'}}->{$row->{'tag_id'}}});
	} else {
	    $pro_cnt = 0;
	}

	#
	# Determine the number of genotypes associated with this RAD-Tag
	#
	if (defined($genotypes{$row->{'batch_id'}}->{$row->{'tag_id'}})) {
	    $geno_cnt = $genotypes{$row->{'batch_id'}}->{$row->{'tag_id'}};
	} else {
	    $geno_cnt = 0;
	}

	#
	# Determine the number of alleles for this catalog RAD-Tag
	#
	if (defined($alleles{$row->{'batch_id'}}->{$row->{'tag_id'}})) {
	    $allele_cnt = scalar(keys %{$alleles{$row->{'batch_id'}}->{$row->{'tag_id'}}});
	} else {
	    $allele_cnt = 0;
	}

	#
	# Annotate the RAD site
	#
	if ($radome_index) {
	    my $key = $row->{'chr'} . "|" . $row->{'bp'} . "|" . $row->{'strand'};
	    
	    if (defined($radome{$key})) {
		$ref_type = $radome{$key}->{'type'};
		$ref_id   = $radome{$key}->{'id'};
	    } else {
		#
		# Check for a read aligned on the other strand.
		#
		$bp  = $row->{'bp'} - 4;
		$key = $row->{'chr'} . "|" . $bp . "|" . $row->{'strand'};

		if (defined($radome{$key})) {
		    $ref_type = $radome{$key}->{'type'};
		    $ref_id   = $radome{$key}->{'id'}
		} else {
		    $ref_type = "genomic";
		    $ref_id   = 0;
		}
	    }
	}

        #
        # Determine if there are any sequences associated with this marker
        #
        $ests       = 0;
        $pe_radtags = 0;
        $blast_hits = 0;
        if (defined($seqs{$row->{'batch_id'}}->{$row->{'tag_id'}}->{'est'})) {
            $ests = $seqs{$row->{'batch_id'}}->{$row->{'tag_id'}}->{'est'};
        }
        if (defined($seqs{$row->{'batch_id'}}->{$row->{'tag_id'}}->{'pe_radtag'})) {
            $pe_radtags = $seqs{$row->{'batch_id'}}->{$row->{'tag_id'}}->{'pe_radtag'};
        }
        if (defined($hits{$row->{'batch_id'}}->{$row->{'tag_id'}})) {
            $blast_hits = $hits{$row->{'batch_id'}}->{$row->{'tag_id'}};
        }

	#
	# Does this RAD-Tag have a mappable marker?
	#
	if (defined($markers{$row->{'batch_id'}}->{$row->{'tag_id'}})) {
	    $tag          = $markers{$row->{'batch_id'}}->{$row->{'tag_id'}};
	    $marker       = $tag->{'marker'};
	    $uncor_marker = $tag->{'uncor_marker'};
	    $chisq_pval   = $tag->{'chisq_pval'};
	    $valid_pro    = $tag->{'valid_pro'};
	    $ratio        = $tag->{'ratio'};
	    $lnl          = $tag->{'lnl'};
	} else {
	    $marker       = "";
	    $uncor_marker = "";
	    $valid_pro    = 0;
	    $chisq_pval   = 1.0;
	    $ratio        = "";
	    $lnl          = 0.0;
	}

	#
	# Record the chromosomes present in this dataset (if any).
	#
	if (length($row->{'chr'}) > 0) {
	    $chrs{$row->{'batch_id'}}->{$row->{'chr'}} = 
		$row->{'bp'} > $chrs{$row->{'batch_id'}}->{$row->{'chr'}} ?
		$row->{'bp'} : $chrs{$row->{'batch_id'}}->{$row->{'chr'}};
	}

	print $fh
	    "0\t",
	    $row->{'batch_id'}, "\t",
	    $row->{'id'}, "\t",
	    $row->{'tag_id'}, "\t",
	    $count, "\t",
	    $par_cnt, "\t",
	    $pro_cnt, "\t",
	    $allele_cnt, "\t",
	    $marker, "\t",
	    $uncor_marker, "\t",
	    $valid_pro, "\t",
	    $chisq_pval, "\t",
	    $lnl, "\t",
	    $ratio, "\t",
            $ests, "\t",
            $pe_radtags, "\t",
            $blast_hits, "\t",
	    $geno_cnt, "\t",
	    $row->{'chr'}, "\t",
	    $row->{'bp'}, "\t",
	    $ref_type, "\t",
	    $ref_id, "\n";
    }

    close($fh);

    `mysql --defaults-file=$cnf $db -e "DROP TABLE IF EXISTS catalog_index"`;
    `mysql --defaults-file=$cnf $db < $sql_cat_table`;

    import_sql_file($catalog_file, 'catalog_index');

    if (scalar(keys %chrs) > 0) {
	my ($batch_id, $chr, $max);
	my ($fh, $chr_file) = tempfile("chr_index_XXXXXXXX", UNLINK => 1, TMPDIR => 1);

	foreach $batch_id (sort keys %chrs) {
	    foreach $chr (sort keys %{$chrs{$batch_id}}) {
		#
		# Round up the maximum chromosome length to the nearest megabase.
		#
		$max  = int($chrs{$batch_id}->{$chr} / 1000000);
		$max += $chrs{$batch_id}->{$chr} % 1000000 > 0 ? 1 : 0;

		print $fh 
		    "0\t",
		    $batch_id, "\t",
		    $chr, "\t",
		    $max, "\n";
	    }
	}

	`mysql --defaults-file=$cnf $db -e "DROP TABLE IF EXISTS chr_index"`;
	`mysql --defaults-file=$cnf $db < $sql_chr_table`;

	import_sql_file($chr_file, 'chr_index');
	close($fh);
    }
}

sub fetch_catalog_snps {
    my ($sth, $snps) = @_;

    my ($row);

    $sth->{'cat_snps'}->execute()
	or die("Unable to select results from $db.\n");

    while ($row = $sth->{'cat_snps'}->fetchrow_hashref()) {
	if (!defined($snps->{$row->{'batch_id'}})) {
	    $snps->{$row->{'batch_id'}} = {};
	}

	$snps->{$row->{'batch_id'}}->{$row->{'tag_id'}}++;
    }
}

sub fetch_genotypes {
    my ($sth, $genotypes) = @_;

    my ($row);

    $sth->{'cat_geno'}->execute()
	or die("Unable to select results from $db.\n");

    while ($row = $sth->{'cat_geno'}->fetchrow_hashref()) {
	if (!defined($genotypes->{$row->{'batch_id'}})) {
	    $genotypes->{$row->{'batch_id'}} = {};
	}

	if ($row->{'genotype'} ne "-" && $row->{'genotype'} ne "--") {
	    $genotypes->{$row->{'batch_id'}}->{$row->{'tag_id'}}++;
	}
    }
}

sub fetch_markers {
    my ($sth, $markers) = @_;

    my ($row, $tag);

    $sth->{'marker'}->execute()
	or die("Unable to select results from $db.\n");

    while ($row = $sth->{'marker'}->fetchrow_hashref()) {
	$tag = {};
	$tag->{'marker'}       = $row->{'type'};
	$tag->{'uncor_marker'} = $row->{'uncor_type'};
	$tag->{'chisq_pval'}   = $row->{'chisq_pval'};
	$tag->{'valid_pro'}    = $row->{'progeny'};
	$tag->{'ratio'}        = $row->{'ratio'};
	$tag->{'lnl'}          = $row->{'lnl'};

	if (!defined($markers->{$row->{'batch_id'}})) {
	    $markers->{$row->{'batch_id'}} = {};
	}
	$markers->{$row->{'batch_id'}}->{$row->{'catalog_id'}} = $tag;
    }
}

sub radome_ref {
    my ($sth, $radome) = @_;

    my ($row, $key, $strand);

    $sth->{'radome'}->execute()
	or die("Unable to select results from $db.\n");

    while ($row = $sth->{'radome'}->fetchrow_hashref()) {
	$strand = $row->{'strand'} == 1 ? "+" : "-";
	$key    = $row->{'chr'} . "|" . $row->{'bp'} . "|" . $strand;

	if (!defined($radome->{$key}) || 
	    $row->{'type'} eq "exon") {
	    $radome->{$key} = {'type' => $row->{'type'},
			       'id'   => $row->{'id'}};
	}
    }
}

sub catalog_matches {
    my ($sth, $parents, $progeny, $alleles) = @_;

    my ($row, $key);

    $sth->{'cat_matches'}->execute()
	or die("Unable to select results from $db.\n");

    while ($row = $sth->{'cat_matches'}->fetchrow_hashref()) {
	$key = $row->{'sample_id'} . "_" . $row->{'tag_id'};

	if ($row->{'type'} eq "parent" || $row->{'type'} eq "sample") {
	    if (!defined($parents->{$row->{'batch_id'}})) {
		$parents->{$row->{'batch_id'}} = {};
	    }
	    if (!defined($parents->{$row->{'batch_id'}}->{$row->{'catalog_id'}})) {
		$parents->{$row->{'batch_id'}}->{$row->{'catalog_id'}} = {};
	    }
	    $parents->{$row->{'batch_id'}}->{$row->{'catalog_id'}}->{$key}++;

	} elsif ($row->{'type'} eq "progeny") {
	    if (!defined($progeny->{$row->{'batch_id'}})) {
		$progeny->{$row->{'batch_id'}} = {};
	    }
	    if (!defined($progeny->{$row->{'batch_id'}}->{$row->{'catalog_id'}})) {
		$progeny->{$row->{'batch_id'}}->{$row->{'catalog_id'}} = {};
	    }
	    $progeny->{$row->{'batch_id'}}->{$row->{'catalog_id'}}->{$key}++;
	}

	if (!defined($alleles->{$row->{'batch_id'}})) {
	    $alleles->{$row->{'batch_id'}} = {};
	}
	if (!defined($alleles->{$row->{'batch_id'}}->{$row->{'catalog_id'}})) {
	    $alleles->{$row->{'batch_id'}}->{$row->{'catalog_id'}} = {};
	}
	$alleles->{$row->{'batch_id'}}->{$row->{'catalog_id'}}->{$row->{'allele'}}++;
    }
}

sub sequence_matches {
    my ($sth, $seqs, $hits) = @_;

    my ($row);

    $sth->{'cat_seqs'}->execute()
	or die("Unable to select results from $db.\n");

    while ($row = $sth->{'cat_seqs'}->fetchrow_hashref()) {
	if (!defined($seqs->{$row->{'batch_id'}})) {
	    $seqs->{$row->{'batch_id'}} = {};
	}
        if (!defined($seqs->{$row->{'batch_id'}}->{$row->{'catalog_id'}})) {
            $seqs->{$row->{'batch_id'}}->{$row->{'catalog_id'}} = {};
        }

        $seqs->{$row->{'batch_id'}}->{$row->{'catalog_id'}}->{$row->{'type'}}++;
    }

    $sth->{'cat_hits'}->execute()
	or die("Unable to select results from $db.\n");

    while ($row = $sth->{'cat_hits'}->fetchrow_hashref()) {
	if (!defined($hits->{$row->{'batch_id'}})) {
	    $hits->{$row->{'batch_id'}} = {};
	}
        $hits->{$row->{'batch_id'}}->{$row->{'catalog_id'}}++;
    }
}

sub gen_tag_index {
    my ($sth) = @_;

    die ("Unable to find SQL file: '$sql_tag_table'\n") if (!-e $sql_tag_table);

    my ($fh, $tag_file) = tempfile("tag_index_XXXXXXXX", UNLINK => 1, TMPDIR => 1);

    print STDERR "Generating unique tag index...\n";

    my ($sample_row, $tags_row, $row, $catalog_id, $i, $num_samples);

    $sth->{'sample'}->execute()
	or die("Unable to select results from $db.\n");
    
    $num_samples = $sth->{'sample'}->rows();
    $i           = 1;

    while ($sample_row = $sth->{'sample'}->fetchrow_hashref()) {

	print STDERR "Processing sample $i of $num_samples       \r";

	my (%depth, %snps, %cats);

	$sth->{'tags'}->execute($sample_row->{'batch_id'}, $sample_row->{'sample_id'})
	    or die("Unable to select results from $db.\n");

	fetch_depth_counts($sth, \%depth, $sample_row->{'sample_id'});
	fetch_snp_counts($sth, \%snps, $sample_row->{'sample_id'});
	fetch_catalog_ids($sth, \%cats, $sample_row->{'batch_id'}, $sample_row->{'sample_id'});

	while ($tags_row = $sth->{'tags'}->fetchrow_hashref()) {

	    print $fh
		"0\t",
		$sample_row->{'batch_id'}, "\t",
		$sample_row->{'sample_id'}, "\t",
		$tags_row->{'tag_id'}, "\t",
		$tags_row->{'id'}, "\t",
		$depth{$tags_row->{'tag_id'}}, "\t",
		defined($snps{$tags_row->{'tag_id'}}) ? $snps{$tags_row->{'tag_id'}} : 0, "\t",
		$cats{$tags_row->{'tag_id'}}, "\t",
		$tags_row->{'deleveraged'}, "\t",
		$tags_row->{'blacklisted'}, "\t",
		$tags_row->{'removed'}, "\n";
	}

	$i++;
    }

    print STDERR "\n";

    close($fh);

    `mysql --defaults-file=$cnf $db -e "DROP TABLE IF EXISTS tag_index"`;
    `mysql --defaults-file=$cnf $db < $sql_tag_table`;

    import_sql_file($tag_file, 'tag_index');
}

sub fetch_depth_counts {
    my ($sth, $depths, $sample_id) = @_;

    my ($row);

    #
    # Determine the depth of coverage for the RAD-Tags in this sample
    #
    $sth->{'depth'}->execute($sample_id)
	or die("Unable to select results from $db.\n");

    while ($row = $sth->{'depth'}->fetchrow_hashref()) {
	$depths->{$row->{'tag_id'}}++;
    }
}

sub fetch_snp_counts {
    my ($sth, $snps, $sample_id) = @_;

    my ($row);

    #
    # Determine the number of SNPs contained within each RAD-Tag in this sample.
    #
    $sth->{'snps'}->execute($sample_id)
	or die("Unable to select results from $db.\n");

    while ($row = $sth->{'snps'}->fetchrow_hashref()) {
	$snps->{$row->{'tag_id'}}++;
    }
}

sub fetch_catalog_ids {
    my ($sth, $cats, $batch_id, $sample_id) = @_;

    my ($row);

    #
    # Determine the catalog ID that corresponds to the RAD-Tags in this sample
    #
    $sth->{'match'}->execute($batch_id, $sample_id)
	or die("Unable to select results from $db.\n");

    while ($row = $sth->{'match'}->fetchrow_hashref()) {
	$cats->{$row->{'tag_id'}} = $row->{'catalog_id'};
    }
}

sub import_sql_file {
    my ($file, $table) = @_;

    my (@results);

    @results = `mysql --defaults-file=$cnf $db -e "LOAD DATA LOCAL INFILE '$file' INTO TABLE $table"`;
    print STDERR "mysql --defaults-file=$cnf $db -e \"LOAD DATA LOCAL INFILE '$file' INTO TABLE $table\"\n", @results;
}

sub prepare_sql_handles {
    my ($sth, $outg) = @_;

    #
    # Connect to the database, check for the existence of a MySQL config file in the home
    # directory first, otherwise use the stacks-distributed one.
    #
    $sth->{'dbh'} = DBI->connect("DBI:mysql:$db:mysql_read_default_file=$cnf")
	or die("Unable to connect to the $db MySQL Database!\n" . $DBI::errstr);

    my $query;

    $query = 
	"SELECT batch_id, id as sample_id, type FROM samples";
    $sth->{'sample'} = $sth->{'dbh'}->prepare($query) or die($sth->{'dbh'}->errstr());

    $query = 
	"SELECT unique_tags.id as id, tag_id, seq, deleveraged, blacklisted, removed FROM unique_tags " . 
	"JOIN samples ON (unique_tags.sample_id=samples.id) " . 
	"WHERE relationship='consensus' AND samples.batch_id=? AND unique_tags.sample_id=?";
    $sth->{'tags'} = $sth->{'dbh'}->prepare($query) or die($sth->{'dbh'}->errstr());

    $query = 
	"SELECT tag_id FROM unique_tags " . 
	"WHERE relationship!='consensus' AND relationship != 'model' AND unique_tags.sample_id=?";
    $sth->{'depth'} = $sth->{'dbh'}->prepare($query) or die($sth->{'dbh'}->errstr());

    $query = 
	"SELECT tag_id FROM snps " . 
	"JOIN samples ON (snps.sample_id=samples.id) " . 
	"WHERE snps.type='E' AND samples.id=?";
    $sth->{'snps'} = $sth->{'dbh'}->prepare($query) or die($sth->{'dbh'}->errstr());

    $query = 
	"SELECT tag_id, catalog_id FROM matches " . 
	"WHERE batch_id=? AND sample_id=?";
    $sth->{'match'} = $sth->{'dbh'}->prepare($query) or die($sth->{'dbh'}->errstr());

    $query = 
	"SELECT catalog_tags.id as id, tag_id, batch_id, chr, bp, strand, seq FROM catalog_tags " . 
	"JOIN batches ON (catalog_tags.batch_id=batches.id) WHERE relationship='consensus'";
    $sth->{'cat_tags'} = $sth->{'dbh'}->prepare($query) or die($sth->{'dbh'}->errstr());

    $query = 
	"SELECT batch_id, tag_id FROM catalog_snps WHERE type='E'";
    $sth->{'cat_snps'} = $sth->{'dbh'}->prepare($query) or die($sth->{'dbh'}->errstr());

    $query = 
	"SELECT batch_id, catalog_id as tag_id, sample_id, genotype FROM catalog_genotypes";
    $sth->{'cat_geno'} = $sth->{'dbh'}->prepare($query) or die($sth->{'dbh'}->errstr());

    $query = 
	"SELECT batch_id, catalog_id, type, uncor_type, progeny, chisq_pval, ratio, lnl FROM markers";
    $sth->{'marker'} = $sth->{'dbh'}->prepare($query) or die($sth->{'dbh'}->errstr());

    $query = 
	"SELECT samples.batch_id, catalog_id, tag_id, matches.sample_id, allele, type FROM matches " . 
	"JOIN samples ON (samples.id=matches.sample_id)";
    $sth->{'cat_matches'} = $sth->{'dbh'}->prepare($query) or die($sth->{'dbh'}->errstr());
#	"JOIN samples ON (samples.sample_id=matches.sample_id AND samples.batch_id=matches.batch_id)";

    $query = 
	"SELECT batch_id, catalog_id, type FROM sequence";
    $sth->{'cat_seqs'} = $sth->{'dbh'}->prepare($query) or die($sth->{'dbh'}->errstr());

    $query = 
	"SELECT batch_id, catalog_id FROM sequence_blast";
    $sth->{'cat_hits'} = $sth->{'dbh'}->prepare($query) or die($sth->{'dbh'}->errstr());

   $query = 
	"SELECT id, chr, bp, strand, type FROM ref_radome";
    $sth->{'radome'} = $sth->{'dbh'}->prepare($query) or die($sth->{'dbh'}->errstr());
}

sub close_sql_handles {
    my ($sth) = @_;

    my $key;

    foreach $key (keys %{$sth}) {
	next if ($key =~ /dbh/);

	$sth->{$key}->finish();
    }

    foreach $key (keys %{$sth}) {
	next if ($key !~ /dbh/);

	$sth->{$key}->disconnect();
    }
}

sub parse_command_line {
    while (@ARGV) {
	$_ = shift @ARGV;
	if    ($_ =~ /^-d$/) { $debug++; }
	elsif ($_ =~ /^-D$/) { $db       = shift @ARGV; }
	elsif ($_ =~ /^-s$/) { $sql_path = shift @ARGV; }
	elsif ($_ =~ /^-c$/) { $catalog_index++; }
	elsif ($_ =~ /^-t$/) { $tag_index++; }
	elsif ($_ =~ /^-r$/) { $radome_index++; }
	elsif ($_ =~ /^-v$/) { version(); exit(); }
	elsif ($_ =~ /^-h$/) { usage(); }
	else {
	    print STDERR "Unknown command line options received: $_\n";
	    usage();
	}
    }

    if (length($db) == 0) {
        print STDERR "You must specify a database to index.\n";
        usage();
    }

    $sql_path .= "/" if (substr($sql_path, -1, 1) ne "/");
    $sql_tag_table = $sql_path . "tag_index.sql";
    $sql_cat_table = $sql_path . "catalog_index.sql";
    $sql_chr_table = $sql_path . "chr_index.sql";
}

sub version {
    print STDERR "index_radtags.pl ", stacks_version, "\n";
}

sub usage {
    version();

    print << "EOQ";
index_radtags.pl -D db [-c] [-t] [-s path] [-d] [-h]
  D: radtag database to examine.
  c: generate a catalog index.
  t: generate a unique tags index.
  s: path to SQL definition files for catalog/tag index tables (if not in default, installed location).
  h: display this help message.
  d: turn on debug output.

EOQ

    exit(0);
}
