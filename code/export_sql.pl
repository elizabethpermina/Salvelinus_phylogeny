#!/usr/bin/env perl
#
# Copyright 2010-2015, Julian Catchen <jcatchen@illinois.edu>
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
# Read in a set of filtering paramters, query a Stacks pipeline database based on
# those filters and write the results into a compact tab-separated values or excel file.
#

use strict;
use DBI;
#use Excel::Writer::XLSX;
use Spreadsheet::WriteExcel;

use constant stacks_version => "2.0Beta7c";

my $mysql_config = "/home/epermina/soft/stacks-2.0Beta7c/share/stacks/" . "sql/mysql.cnf";
my $out_file  = "";
my $type      = "tsv";
my $batch_id  = 0;
my $data_type = "haplo";
my $map_type  = "gen";
my $all_depth = 0;
my $allele_depth_lim = 1;
my $locus_depth_lim  = 0;
my $locus_lnl_lim    = -10000.0;
my $man_cor   = 0;
my $db        = "";

my $translate_genotypes = {'dh'  => \&trans_dh_map,
			   'cp'  => \&trans_cp_map,
			   'bc1' => \&trans_bc1_map,
			   'f2'  => \&trans_f2_map,
			   'gen' => \&trans_gen_map};

my @valid_filters = ("cata", "alle_l", "alle_u", "snps_l", "snps_u", "pare_l", "pare_u",
		     "prog",  "vprog", "mark", "est",  "pe", "blast", "gcnt",
		     "chisq_l", "chisq_u", "lnl_l", "lnl_u", "ref", "loc");

parse_command_line();

my (%sth, %loci, %samples, %depths, %filters);

prepare_sql_handles(\%sth, \%filters);

populate(\%sth, \%loci, \%samples, \%depths, \%filters);

apply_corrected_genotypes(\%sth, \%loci) if ($man_cor > 0);

if ($data_type eq "haplo") {
    write_observed_haplotypes(\%loci, \%samples);

} elsif ($data_type eq "geno") {
    write_genotypes(\%loci, \%samples, \%depths);
}

print "Success\n";

sub populate {
    my ($sth, $loci, $samples, $depths, $filters) = @_;

    my (%delev);
    my ($row, $snp_row, $all_row, $gen_row, $locus);
    my @params;

    #
    # Cache all the stacks that were deleveraged.
    #
    $sth->{'delev'}->execute($batch_id);

    while ($row = $sth->{'delev'}->fetchrow_hashref()) {
        $delev{$row->{'sample_id'} . "_" . $row->{'tag_id'}}++;
    }

    #
    # Pull list of samples for this batch
    #
    $sth->{'samp'}->execute($batch_id);

    while ($row = $sth->{'samp'}->fetchrow_hashref()) {
        $samples->{$row->{'file'}} = $row->{'id'};
    }

    prepare_filter_parameters(\@params, $filters);

    #
    # Fetch the results and populate the array of groups.
    #
    $sth->{'tag'}->execute(@params);

    while ($row = $sth->{'tag'}->fetchrow_hashref()) {
        $locus = {};
        $locus->{'id'}         = $row->{'tag_id'};
        $locus->{'annotation'} = defined($row->{'external_id'}) ? $row->{'external_id'} : "";
        $locus->{'chr'}        = $row->{'chr'};
        $locus->{'bp'}         = $row->{'bp'};
        $locus->{'delev'}      = 0;
        $locus->{'marker'}     = $row->{'marker'};
        $locus->{'seq'}        = $row->{'seq'};
        $locus->{'alleles'}    = "";
        $locus->{'snps'}       = "";
        $locus->{'gtypes'}     = {};

        $locus->{'num_alleles'}   = $row->{'alleles'};
        $locus->{'num_snps'}      = $row->{'snps'};
        $locus->{'num_parents'}   = $row->{'parents'};
        $locus->{'num_progeny'}   = $row->{'progeny'};
        $locus->{'valid_progeny'} = $row->{'valid_progeny'};
        $locus->{'num_ests'}      = $row->{'ests'};
        $locus->{'num_pe_tags'}   = $row->{'pe_radtags'};
        $locus->{'num_blast'}     = $row->{'blast_hits'};
	$locus->{'gcnt'}          = $row->{'geno_cnt'};

        $loci->{$row->{'tag_id'}} = $locus;
    }


    if ($data_type eq "haplo") {
	#
	# Add observed haplotypes
	#
	$sth->{'mat'}->execute($batch_id, $allele_depth_lim);

	while ($gen_row = $sth->{'mat'}->fetchrow_hashref()) {
	    next if (!defined($loci->{$gen_row->{'catalog_id'}}));
	    $locus = $loci->{$gen_row->{'catalog_id'}};

	    if (!defined($locus->{'gtypes'}->{$gen_row->{'file'}})) {
		$locus->{'gtypes'}->{$gen_row->{'file'}} = [];
	    }

	    push(@{$locus->{'gtypes'}->{$gen_row->{'file'}}}, 
		 {'file'   => $gen_row->{'file'},
		  'allele' => $gen_row->{'allele'},
		  'tag_id' => $gen_row->{'tag_id'},
		  'depth'  => $gen_row->{'depth'},
		  'lnl'    => $gen_row->{'lnl'}});

	    #
	    # Check if this particular sample was deleveraged
	    #
	    if (defined($delev{$gen_row->{'id'} . "_" . $gen_row->{'tag_id'}}) &&
		$delev{$gen_row->{'id'} . "_" . $gen_row->{'tag_id'}} >= 1) {
		$locus->{'delev'}++;
	    }
	}
    } elsif ($data_type eq "geno") {
	#
	# Add genotypes
	#
	$sth->{'gtypes'}->execute($batch_id);

	while ($gen_row = $sth->{'gtypes'}->fetchrow_hashref()) {
	    next if (!defined($loci->{$gen_row->{'catalog_id'}}));
	    $locus = $loci->{$gen_row->{'catalog_id'}};

	    if (!defined($locus->{'gtypes'}->{$gen_row->{'file'}})) {
		$locus->{'gtypes'}->{$gen_row->{'file'}} = [];
	    }

	    push(@{$locus->{'gtypes'}->{$gen_row->{'file'}}}, 
		 {'file'   => $gen_row->{'file'},
		  'gtype'  => $gen_row->{'genotype'}});
	}
    }

    #
    # Fetch SNPs and Alleles
    #
    $sth->{'snp'}->execute($batch_id);

    while ($snp_row = $sth->{'snp'}->fetchrow_hashref()) {
	next if (!defined($loci->{$snp_row->{'tag_id'}}));

	$loci->{$snp_row->{'tag_id'}}->{'snps'} .= $snp_row->{'col'} . "," . $snp_row->{'rank_1'} . ">" . $snp_row->{'rank_2'} . ";";
    }

    $sth->{'allele'}->execute($batch_id);

    while ($all_row = $sth->{'allele'}->fetchrow_hashref()) {
	next if (!defined($loci->{$all_row->{'tag_id'}}));

	$loci->{$all_row->{'tag_id'}}->{'alleles'} .= $all_row->{'allele'} . ";";
    }

    #
    # If exporting genotypes and a locus depth limit was specified, fetch locus depths.
    #
    if ($data_type eq "geno" && $locus_depth_lim > 0) {

	$sth->{'depths'}->execute($batch_id);

	while ($row = $sth->{'depths'}->fetchrow_hashref()) {
	    next if (!defined($loci->{$row->{'catalog_id'}}));

	    if (!defined($depths->{$row->{'catalog_id'}})) {
		$depths->{$row->{'catalog_id'}} = {};
	    }
	    $depths->{$row->{'catalog_id'}}->{$row->{'file'}} += $row->{'depth'};
	}
    }
}

sub apply_corrected_genotypes {
    my ($sth, $loci) = @_;

    my (%corrections, $locus, $key, $row, $sample);

    #print STDERR "Applying manually corrected genotypes to export data...\n";

    #
    # Fetch the manual corrections from the database.
    #
    $sth->{'corr'}->execute($batch_id)
	or die("Unable to select results from $db.\n");

    while ($row = $sth->{'corr'}->fetchrow_hashref()) {
        if (!defined($corrections{$row->{'catalog_id'}})) {
            $corrections{$row->{'catalog_id'}} = {};
        }
        $corrections{$row->{'catalog_id'}}->{$row->{'file'}} = $row->{'genotype'};
    }

    foreach $key (keys %{$loci}) {
	$locus = $loci->{$key};

        next if (!defined($corrections{$locus->{'id'}}));

        foreach $sample (keys %{$corrections{$locus->{'id'}}}) {
	    @{$locus->{'gtypes'}->{$sample}} = ();
	    push(@{$locus->{'gtypes'}->{$sample}}, 
		     {'file'   => $sample,
		      'gtype'  => $corrections{$locus->{'id'}}->{$sample}});
        }
    }
}

sub prepare_filter_parameters {
    my ($params, $filters) = @_;

    my ($filter);

    push(@{$params}, $batch_id);

    foreach $filter (keys %{$filters}) {

        if ($filter eq "snps") {
            push(@{$params}, $filters->{'snps_l'});
            push(@{$params}, $filters->{'snps_u'});

        } elsif ($filter eq "alle") {
            push(@{$params}, $filters->{'alle_l'});
            push(@{$params}, $filters->{'alle_u'});

        } elsif ($filter eq "pare") {
            push(@{$params}, $filters->{'pare_l'});
            push(@{$params}, $filters->{'pare_u'});

        } elsif ($filter eq "lnl") {
            push(@{$params}, $filters->{'lnl_l'});
            push(@{$params}, $filters->{'lnl_u'});

        } elsif ($filter eq "chisq") {
            push(@{$params}, $filters->{'chisq_l'});
            push(@{$params}, $filters->{'chisq_u'});

        } elsif ($filter eq "prog") {
            push(@{$params}, $filters->{'prog'});

        } elsif ($filter eq "vprog") {
            push(@{$params}, $filters->{'vprog'});

        } elsif ($filter eq "cata") {
            push(@{$params}, $filters->{'cata'});

        } elsif ($filter eq "gcnt") {
            push(@{$params}, $filters->{'gcnt'});

        } elsif ($filter eq "est") {
            push(@{$params}, 0);

        } elsif ($filter eq "pe") {
            push(@{$params}, 0);

        } elsif ($filter eq "blast") {
            push(@{$params}, 0);

        } elsif ($filter eq "ref") {
            push(@{$params}, $filters->{'ref'});

        } elsif ($filter eq "loc") {
	    push(@{$params}, $filters->{'chr'});
	    push(@{$params}, $filters->{'sbp'} * 1000000);
	    push(@{$params}, $filters->{'ebp'} * 1000000);

	} elsif ($filter eq "mark") {
            if ($filters->{'mark'} eq "Any") {
                push(@{$params}, "%/%");
            } else {
                push(@{$params}, $filters->{'mark'});
            }
        }
    }
}

sub apply_query_filters {
    my ($filters) = @_;

    my ($query, $filter) = "";

    my %sql_filters = 
        ("cata"  => "(catalog_index.tag_id = ?)", 
	 "alle"  => "(alleles >= ? AND alleles <= ?)", 
	 "snps"  => "(snps >= ? AND snps <= ?)",
	 "pare"  => "(parents >= ? AND parents <= ?)",
         "prog"  => "(progeny >= ?)",
         "vprog" => "(valid_progeny >= ?)",
	 "lnl"   => "(lnl >= ? AND lnl <= ?)",
         "mark"  => "(marker LIKE ?)", 
         "est"   => "(ests > ?)",
         "pe"    => "(pe_radtags > ?)",
         "blast" => "(blast_hits > ?)",
	 "gcnt"  => "(geno_cnt >= ?)",
	 "chisq" => "(chisq_pval >= ? AND chisq_pval <= ?)",
	 "ref"   => "(catalog_index.type = ?)",
	 "loc"   => "(catalog_index.chr = ? && catalog_index.bp >= ? && catalog_index.bp <= ?)");

    if (scalar(keys %{$filters}) > 0) {

        foreach $filter (keys %{$filters}) {
	    next if (!defined($sql_filters{$filter}));

            $query .= " AND ";
            $query .= $sql_filters{$filter};
        }
    }

    return $query;
}

sub write_observed_haplotypes {
    my ($loci, $samples, $filters) = @_;

    my ($workbook, $worksheet);

    my ($out_fh, $str, $cat_id, $id, $locus, $gtypes, $types, $tot_depth);

    if ($type eq "xls") {
        $workbook  = Spreadsheet::WriteExcel->new($out_file) or die("Unable to initiate excel spreadsheet.\n");
        $worksheet = $workbook->add_worksheet() or die("Unable to add a worksheet to our excel spreadsheet.\n");
    } else {
        open($out_fh, ">$out_file") or die("Unable to open output file '$out_file'\n");
    }

    #
    # Order the samples by sample ID
    #
    my @ordered_sam = sort {$samples->{$a} <=> $samples->{$b}} keys %{$samples};

    #
    # Print the heading out for the spreadsheet
    #
    my $i = 0;

    $str = "# " if ($type ne "xls");
    
    $str .= 
        "Catalog ID\t" . 
        "Annotation\t" . 
        "Chr\t" .
        "BP\t" .
        "Consensus Sequence\t" .
        "Num Parents\t" .
        "Num Progeny\t" . 
        "Num SNPs\t" .
        "SNPs\t" .
        "Num Alleles\t" .
        "Alleles\t" .
        "Deleveraged\t";

    foreach $id (@ordered_sam) {
        $str .= $id . "\t";
    }
    $str  = substr($str, 0, -1);
    $str .= "\n";

    $type eq "xls" ? write_excel($worksheet, $i, $str) : print $out_fh $str;
    $i++;

    foreach $cat_id (keys %{$loci}) {
        $locus = $loci->{$cat_id};

	$locus->{'snps'}    = substr($locus->{'snps'}, 0, -1) if (length($locus->{'snps'}) > 0);
	$locus->{'alleles'} = substr($locus->{'alleles'}, 0, -1) if (length($locus->{'alleles'}) > 0);

        $str =
            $cat_id . "\t" .
            $locus->{'annotation'} . "\t" .
            $locus->{'chr'} . "\t" .
            $locus->{'bp'} . "\t" .
            $locus->{'seq'} . "\t" .
            $locus->{'num_parents'} . "\t" .
            $locus->{'num_progeny'} . "\t" .
            $locus->{'num_snps'} . "\t" .
            $locus->{'snps'} . "\t" . 
            $locus->{'num_alleles'} . "\t" .
            $locus->{'alleles'} . "\t" .
            $locus->{'delev'} . "\t";

        foreach $id (@ordered_sam) {
            $types = $locus->{'gtypes'}->{$id};

            if (!defined($types)) {
                $str .= "\t";
                next;
            }

	    #
	    # Check total locus depth.
	    #
	    $tot_depth = 0;
            foreach $type (@{$types}) {
		$tot_depth += $type->{'depth'};
            }
	    if ($tot_depth < $locus_depth_lim) {
                $str .= "\t";
                next;
	    }

	    if ($types->[0]->{'lnl'} < $locus_lnl_lim) {
	    	$str .= "\t";
	    	next;
	    }

            foreach $type (@{$types}) {
                $str .= $all_depth ? $type->{'depth'} : $type->{'allele'};
		$str .= "/";
            }

            $str  = substr($str, 0, -1);
            $str .= "\t";
        }
        $str  = substr($str, 0, -1);
        $str .= "\n";

        $type eq "xls" ? write_excel($worksheet, $i, $str) : print $out_fh $str;

        $i++;
    }

    $str = "\n";
    $type eq "xls" ? write_excel($worksheet, $i, $str) : print $out_fh $str;
    $i++;

    foreach $id (@ordered_sam) {
        $str = "\t" . $samples->{$id} . "\t" . $id . "\n";
        $type eq "xls" ? write_excel($worksheet, $i, $str) : print $out_fh $str;
        $i++;
    }

    $type eq "xls" ? $workbook->close() : close($out_fh);
}

sub write_genotypes {
    my ($loci, $samples, $depths) = @_;

    my ($workbook, $worksheet);

    my ($out_fh, $str, $cat_id, $id, $locus, $gtypes, $types);

    if ($type eq "xls") {
        $workbook  = Spreadsheet::WriteExcel->new($out_file) or die("Unable to initiate excel spreadsheet.\n");
        $worksheet = $workbook->add_worksheet() or die("Unable to add a worksheet to our excel spreadsheet.\n");
    } else {
        open($out_fh, ">$out_file") or die("Unable to open output file '$out_file'\n");
    }

    #
    # Order the samples by sample ID
    #
    my @ordered_sam = sort {$samples->{$a} <=> $samples->{$b}} keys %{$samples};

    #
    # Print the heading out for the spreadsheet
    #
    my $i = 0;

    $str = 
        "Catalog ID\t" . 
        "Annotation\t" . 
        "Chr\t" .
        "BP\t" .
	"Marker\t";

    foreach $id (@ordered_sam) {
        $str .= $id . "\t";
    }
    $str  = substr($str, 0, -1);
    $str .= "\n";

    $type eq "xls" ? write_excel($worksheet, $i, $str) : print $out_fh $str;
    $i++;

    my ($trans_marker);

    foreach $cat_id (keys %{$loci}) {
        $locus = $loci->{$cat_id};

	$trans_marker = translate_marker($map_type, $locus->{'marker'});

        $str =
            $cat_id . "\t" .
            $locus->{'annotation'} . "\t" .
            $locus->{'chr'} . "\t" .
            $locus->{'bp'} . "\t" .
	    $trans_marker . "\t";

        foreach $id (@ordered_sam) {
            $types = $locus->{'gtypes'}->{$id};

	    #
	    # Check that there is a genotype for this sample.
	    #
            if (!defined($types)) {
                $str .= "\t";
                next;
            }

	    #
	    # Check the depth of coverage for this locus, if requested.
	    #
	    if ($locus_depth_lim > 0) {
		if (!defined($depths->{$cat_id}->{$id}) || $depths->{$cat_id}->{$id} < $locus_depth_lim) {
		    $str .= "\t";
		    next;
		}
	    }

	    my $trans_gtype = $translate_genotypes->{$map_type}->($trans_marker, $types->[0]->{'gtype'});
            $str .=  $trans_gtype . "\t";
        }
        $str  = substr($str, 0, -1);
        $str .= "\n";

        $type eq "xls" ? write_excel($worksheet, $i, $str) : print $out_fh $str;

        $i++;
    }

    $str = "\n";
    $type eq "xls" ? write_excel($worksheet, $i, $str) : print $out_fh $str;
    $i++;

    foreach $id (@ordered_sam) {
        $str = $samples->{$id} . "\t" . $id . "\n";
        $type eq "xls" ? write_excel($worksheet, $i, $str) : print $out_fh $str;
        $i++;
    }

    $type eq "xls" ? $workbook->close() : close($out_fh);
}

sub write_excel {
    my ($worksheet, $i, $str) = @_;

    chomp $str;
    my @row = split(/\t/, $str);

    my $j = 0;

    foreach my $r (@row) {
        $worksheet->write($i, $j, $r);
        $j++;
    }
}

sub translate_marker {
    my ($map_type, $in_marker) = @_;

    my %dictionary;

    return $in_marker if ($map_type eq "gen");

    $dictionary{"dh"}  = {};
    $dictionary{"cp"}  = {};
    $dictionary{"f2"}  = {};
    $dictionary{"bc1"} = {};

    $dictionary{"dh"}->{"ab/--"}  = "abx--";
    $dictionary{"dh"}->{"--/ab"}  = "--xab";

    $dictionary{"cp"}->{"ab/--"}  = "lmx--";
    $dictionary{"cp"}->{"--/ab"}  = "--xnp";
    $dictionary{"cp"}->{"ab/aa"}  = "lmxll";
    $dictionary{"cp"}->{"aa/ab"}  = "nnxnp";
    $dictionary{"cp"}->{"ab/ab"}  = "hkxhk";
    $dictionary{"cp"}->{"ab/ac"}  = "efxeg";
    $dictionary{"cp"}->{"ab/cd"}  = "abxcd";

    $dictionary{"f2"}->{"aa/bb"}  = "aaxbb";
    $dictionary{"f2"}->{"ab/cd"}  = "abxcd";
    $dictionary{"f2"}->{"ab/aa"}  = "abxaa";
    $dictionary{"f2"}->{"aa/ab"}  = "aaxab";
    $dictionary{"f2"}->{"ab/cc"}  = "abxcc";
    $dictionary{"f2"}->{"cc/ab"}  = "ccxab";

    $dictionary{"bc1"}->{"aa/bb"} = "aaxbb";
    $dictionary{"bc1"}->{"bb/aa"} = "bbxaa";
    $dictionary{"bc1"}->{"ab/cc"} = "abxcc";
    $dictionary{"bc1"}->{"cc/ab"} = "ccxab";

    return defined($dictionary{$map_type}->{$in_marker}) ? 
	$dictionary{$map_type}->{$in_marker} : "";
}

sub trans_bc1_map {
    my ($marker, $in_gtype) = @_;

    my (%types, %dictionary);

    $dictionary{"aaxbb"} = {};
    $dictionary{"bbxaa"} = {};
    $dictionary{"abxcc"} = {};
    $dictionary{"ccxab"} = {};

    $dictionary{"aaxbb"}->{"--"} = "-";
    $dictionary{"aaxbb"}->{"aa"} = "b";
    $dictionary{"aaxbb"}->{"ab"} = "h";
    $dictionary{"aaxbb"}->{"bb"} = "h";

    $dictionary{"bbxaa"}->{"--"} = "-";
    $dictionary{"bbxaa"}->{"aa"} = "h";
    $dictionary{"bbxaa"}->{"ab"} = "h";
    $dictionary{"bbxaa"}->{"bb"} = "a";

    $dictionary{"abxcc"}->{"--"} = "-";
    $dictionary{"abxcc"}->{"ac"} = "h";
    $dictionary{"abxcc"}->{"bc"} = "h";
    $dictionary{"abxcc"}->{"ab"} = "b";
    $dictionary{"abxcc"}->{"aa"} = "b";
    $dictionary{"abxcc"}->{"bb"} = "b";

    $dictionary{"ccxab"}->{"--"} = "-";
    $dictionary{"ccxab"}->{"ac"} = "h";
    $dictionary{"ccxab"}->{"bc"} = "h";
    $dictionary{"ccxab"}->{"ab"} = "a";
    $dictionary{"ccxab"}->{"aa"} = "a";
    $dictionary{"ccxab"}->{"bb"} = "a";

    my $out_gtype = 
	defined($dictionary{$marker}->{lc($in_gtype)}) ? 
	$dictionary{$marker}->{lc($in_gtype)} : 
	"-";

    if (lc($in_gtype) ne $in_gtype) {
	return uc($out_gtype);
    } else {
	return $out_gtype;
    }
}

sub trans_dh_map {
    my ($marker, $in_gtype) = @_;

    my (%types, %dictionary);

    $dictionary{"abx--"} = {};
    $dictionary{"--xab"} = {};

    $dictionary{"abx--"}->{"aa"} = "a";
    $dictionary{"abx--"}->{"bb"} = "b";
    $dictionary{"abx--"}->{"--"} = "-";

    $dictionary{"--xab"}->{"aa"} = "a";
    $dictionary{"--xab"}->{"bb"} = "b";
    $dictionary{"--xab"}->{"--"} = "-";

    my $out_gtype = 
	defined($dictionary{$marker}->{lc($in_gtype)}) ? 
	$dictionary{$marker}->{lc($in_gtype)} : 
	"-";

    if (lc($in_gtype) ne $in_gtype) {
	return uc($out_gtype);
    } else {
	return $out_gtype;
    }
}

sub trans_f2_map {
    my ($marker, $in_gtype) = @_;

    my (%types, %dictionary);

    $dictionary{"aaxbb"} = {};
    $dictionary{"abxcd"} = {};
    $dictionary{"abxaa"} = {};
    $dictionary{"aaxab"} = {};
    $dictionary{"abxcc"} = {};
    $dictionary{"ccxab"} = {};

    $dictionary{"aaxbb"}->{"aa"} = "a";
    $dictionary{"aaxbb"}->{"ab"} = "h";
    $dictionary{"aaxbb"}->{"bb"} = "b";
    $dictionary{"aaxbb"}->{"--"} = "-";

    $dictionary{"abxcd"}->{"aa"} = "a";
    $dictionary{"abxcd"}->{"ab"} = "a";
    $dictionary{"abxcd"}->{"bb"} = "a";
    $dictionary{"abxcd"}->{"cc"} = "b";
    $dictionary{"abxcd"}->{"cd"} = "b";
    $dictionary{"abxcd"}->{"dd"} = "b";
    $dictionary{"abxcd"}->{"ac"} = "h";
    $dictionary{"abxcd"}->{"ad"} = "h";
    $dictionary{"abxcd"}->{"bc"} = "h";
    $dictionary{"abxcd"}->{"bd"} = "h";
    $dictionary{"abxcd"}->{"--"} = "-";

    $dictionary{"abxaa"}->{"aa"} = "-";
    $dictionary{"abxaa"}->{"ab"} = "-";
    $dictionary{"abxaa"}->{"bb"} = "a";
    $dictionary{"abxaa"}->{"--"} = "-";

    $dictionary{"aaxab"}->{"aa"} = "-";
    $dictionary{"aaxab"}->{"ab"} = "-";
    $dictionary{"aaxab"}->{"bb"} = "b";
    $dictionary{"aaxab"}->{"--"} = "-";

    $dictionary{"abxcc"}->{"a"}  = "a";
    $dictionary{"abxcc"}->{"ab"} = "a";
    $dictionary{"abxcc"}->{"bb"} = "a";
    $dictionary{"abxcc"}->{"cc"} = "b";
    $dictionary{"abxcc"}->{"ac"} = "-";
    $dictionary{"abxcc"}->{"bc"} = "-";
    $dictionary{"abxcc"}->{"--"} = "-";

    $dictionary{"ccxab"}->{"aa"} = "b";
    $dictionary{"ccxab"}->{"ab"} = "b";
    $dictionary{"ccxab"}->{"bb"} = "b";
    $dictionary{"ccxab"}->{"cc"} = "a";
    $dictionary{"ccxab"}->{"ac"} = "-";
    $dictionary{"ccxab"}->{"bc"} = "-";
    $dictionary{"ccxab"}->{"--"} = "-";

    my $out_gtype = 
	defined($dictionary{$marker}->{lc($in_gtype)}) ? 
	$dictionary{$marker}->{lc($in_gtype)} : 
	"-";

    if (lc($in_gtype) ne $in_gtype) {
	return uc($out_gtype);
    } else {
	return $out_gtype;
    }
}

sub trans_cp_map {
    my ($marker, $in_gtype) = @_;

    my (%types, %dictionary);

    $dictionary{"lmx--"} = {};
    $dictionary{"--xnp"} = {};
    $dictionary{"lmxll"} = {};
    $dictionary{"nnxnp"} = {};
    $dictionary{"hkxhk"} = {};
    $dictionary{"efxeg"} = {};
    $dictionary{"abxcd"} = {};

    $dictionary{"lmx--"}->{"--"} = "--";
    $dictionary{"lmx--"}->{"aa"} = "ll";
    $dictionary{"lmx--"}->{"bb"} = "lm";
    $dictionary{"lmx--"}->{"ab"} = "lm";

    $dictionary{"--xnp"}->{"--"} = "--";
    $dictionary{"--xnp"}->{"aa"} = "nn";
    $dictionary{"--xnp"}->{"bb"} = "np";
    $dictionary{"--xnp"}->{"ab"} = "np";

    $dictionary{"lmxll"}->{"--"} = "--";
    $dictionary{"lmxll"}->{"aa"} = "ll";
    $dictionary{"lmxll"}->{"ab"} = "lm";

    $dictionary{"nnxnp"}->{"--"} = "--";
    $dictionary{"nnxnp"}->{"aa"} = "nn";
    $dictionary{"nnxnp"}->{"ab"} = "np";

    $dictionary{"hkxhk"}->{"--"} = "--";
    $dictionary{"hkxhk"}->{"ab"} = "hk";
    $dictionary{"hkxhk"}->{"aa"} = "hh";
    $dictionary{"hkxhk"}->{"bb"} = "kk";

    $dictionary{"efxeg"}->{"--"} = "--";
    $dictionary{"efxeg"}->{"ab"} = "ef";
    $dictionary{"efxeg"}->{"ac"} = "eg";
    $dictionary{"efxeg"}->{"bc"} = "fg";
    $dictionary{"efxeg"}->{"aa"} = "ee";

    $dictionary{"abxcd"}->{"--"} = "--";
    $dictionary{"abxcd"}->{"ac"} = "ac";
    $dictionary{"abxcd"}->{"ad"} = "ad";
    $dictionary{"abxcd"}->{"bc"} = "bc";
    $dictionary{"abxcd"}->{"bd"} = "bd";

    my $out_gtype = 
	defined($dictionary{$marker}->{lc($in_gtype)}) ? 
	$dictionary{$marker}->{lc($in_gtype)} : 
	"-";

    if (lc($in_gtype) ne $in_gtype) {
	return uc($out_gtype);
    } else {
	return $out_gtype;
    }
}

sub trans_gen_map {
    my ($marker, $in_gtype) = @_;

    return $in_gtype;
}

sub prepare_sql_handles {
    my ($sth, $filters) = @_;

    #
    # Connect to the database, check for the existence of a MySQL config file in the home
    # directory first, otherwise use the stacks-distributed one.
    #
    my $cnf = (defined($ENV{"HOME"}) && -e $ENV{"HOME"} . "/.my.cnf") ? $ENV{"HOME"} . "/.my.cnf" : $mysql_config;
    $sth->{'dbh'} = DBI->connect("DBI:mysql:$db:mysql_read_default_file=$cnf")
	or die("Unable to connect to the $db MySQL Database!\n" . $DBI::errstr);

    my $query;

    $query = 
        "SELECT sample_id, tag_id FROM tag_index " . 
        "WHERE batch_id=? AND deleveraged=true";
    $sth->{'delev'} = $sth->{'dbh'}->prepare($query) or die($sth->{'dbh'}->errstr());

    $query = 
        "SELECT id, file FROM samples " . 
        "WHERE batch_id=? ORDER BY id";
    $sth->{'samp'} = $sth->{'dbh'}->prepare($query) or die($sth->{'dbh'}->errstr());

    $query = 
        "SELECT tag_id, allele FROM catalog_alleles " . 
        "WHERE batch_id=?";
    $sth->{'allele'} = $sth->{'dbh'}->prepare($query) or die($sth->{'dbh'}->errstr());

    $query = 
        "SELECT tag_id, col, rank_1, rank_2 FROM catalog_snps " . 
        "WHERE batch_id=?";
    $sth->{'snp'} = $sth->{'dbh'}->prepare($query) or die($sth->{'dbh'}->errstr());

    $query = 
        "SELECT catalog_id, samples.id as id, samples.sample_id, samples.type, file, tag_id, allele, depth, lnl " . 
        "FROM matches " . 
        "JOIN samples ON (matches.sample_id=samples.id) " . 
        "WHERE matches.batch_id=? AND matches.depth>?";
    $sth->{'mat'} = $sth->{'dbh'}->prepare($query) or die($sth->{'dbh'}->errstr());

    $query = 
        "SELECT catalog_id, file, depth " . 
        "FROM matches " . 
        "JOIN samples ON (matches.sample_id=samples.id) " . 
        "WHERE matches.batch_id=?";
    $sth->{'depths'} = $sth->{'dbh'}->prepare($query) or die($sth->{'dbh'}->errstr());

    $query = 
        "SELECT catalog_id, samples.id as id, samples.sample_id, samples.type, file, genotype " . 
        "FROM catalog_genotypes " . 
        "JOIN samples ON (catalog_genotypes.sample_id=samples.id) " . 
        "WHERE catalog_genotypes.batch_id=?";
    $sth->{'gtypes'} = $sth->{'dbh'}->prepare($query) or die($sth->{'dbh'}->errstr());

    $query =
	"SELECT gc.catalog_id, gc.sample_id, gc.genotype, file " . 
        "FROM genotype_corrections as gc " . 
        "JOIN samples ON (gc.sample_id=samples.id) " .
        "WHERE gc.batch_id=?";
    $sth->{'corr'} = $sth->{'dbh'}->prepare($query) or die($sth->{'dbh'}->errstr());

    $query = 
        "SELECT catalog_index.tag_id as tag_id, catalog_index.chr, catalog_index.bp, " .
	"snps, alleles, parents, progeny, valid_progeny, " . 
        "seq, marker, chisq_pval, lnl, ratio, ests, pe_radtags, blast_hits, geno_cnt, external_id " .
        "FROM catalog_index " .
        "JOIN catalog_tags ON (catalog_index.cat_id=catalog_tags.id) " . 
        "LEFT JOIN catalog_annotations ON " . 
        "(" . 
        "catalog_index.batch_id=catalog_annotations.batch_id AND " . 
        "catalog_index.tag_id=catalog_annotations.catalog_id" . 
        ") " .
        "WHERE catalog_index.batch_id=?";

    $query .= apply_query_filters($filters);

    $sth->{'tag'} = $sth->{'dbh'}->prepare($query) or die($sth->{'dbh'}->errstr());
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

    my ($filter, $name, $value);

    while (@ARGV) {
        $_ = shift @ARGV;
        if    ($_ =~ /^-f$/) { $out_file  = shift @ARGV; }
        elsif ($_ =~ /^-o$/) { $type      = shift @ARGV; }
        elsif ($_ =~ /^-a$/) { $data_type = lc(shift @ARGV); }
        elsif ($_ =~ /^-b$/) { $batch_id  = shift @ARGV; }
        elsif ($_ =~ /^-D$/) { $db        = shift @ARGV; }
        elsif ($_ =~ /^-m$/) { $map_type  = lc(shift @ARGV); }
        elsif ($_ =~ /^-A$/) { $allele_depth_lim = shift @ARGV; }
        elsif ($_ =~ /^-L$/) { $locus_depth_lim  = shift @ARGV; }
        elsif ($_ =~ /^-I$/) { $locus_lnl_lim    = shift @ARGV; }
        elsif ($_ =~ /^-d$/) { $all_depth++; }
        elsif ($_ =~ /^-c$/) { $man_cor++; }
	elsif ($_ =~ /^-v$/) { version(); exit(); }
        elsif ($_ =~ /^-h$/) { usage(); }
        elsif ($_ =~ /^-F$/) {
            $filter = shift @ARGV;
            ($name, $value) = split(/=/, $filter);

            if (length($name) == 0 || length($value) == 0) {
                print STDERR "Error parsing filter '$filter'\n";
                usage();
            }
            $filters{$name} = $value;

        } else {
            print STDERR "Unknown command line option: '$_'\n";
            usage();
        }
    }

    if (defined($filters{'chr'})) {
	$filters{'loc'} = 1;
	$filters{'sbp'} = 0 if (!defined($filters{'sbp'}));
	$filters{'ebp'} = 500000000 if (!defined($filters{'ebp'}));
    }

    if (defined($filters{'snps_l'}) || defined($filters{'snps_u'})) {
	$filters{'snps'} = 1;
	$filters{'snps_l'} = 1 if (!defined($filters{'snps_l'}));
	$filters{'snps_u'} = 100 if (!defined($filters{'snps_u'}));
    }

    if (defined($filters{'alle_l'}) || defined($filters{'alle_u'})) {
	$filters{'alle'} = 1;
	$filters{'alle_l'} = 1 if (!defined($filters{'alle_l'}));
	$filters{'alle_u'} = 100 if (!defined($filters{'alle_u'}));
    }

    if (defined($filters{'pare_l'}) || defined($filters{'pare_u'})) {
	$filters{'pare'} = 1;
	$filters{'pare_l'} = 1 if (!defined($filters{'pare_l'}));
	$filters{'pare_u'} = 1000 if (!defined($filters{'pare_u'}));
    }

    if (defined($filters{'lnl_l'}) || defined($filters{'lnl_u'})) {
	$filters{'lnl'} = 1;
	$filters{'lnl_l'} = -500 if (!defined($filters{'lnl_l'}));
	$filters{'lnl_u'} =  0   if (!defined($filters{'lnl_u'}));
    }

    if ($out_file eq "") {
	print STDERR "You must specify the file to write data to!\n";
	usage();
    }

    if ($type ne "tsv" && $type ne "xls") {
        print STDERR "Unknown output file type specified '$type'.\n";
        usage();
    }

    if ($data_type ne "haplo" && $data_type ne "geno") {
	print STDERR "Unknown data type specified, 'haplo' and 'geno' are currently accepted.\n";
	usage();
    }

    if ($data_type eq "geno" &&
	$map_type  ne "bc1" && 
	$map_type  ne "dh" && 
	$map_type  ne "f2" && 
	$map_type  ne "cp" && 
	$map_type  ne "gen") {
	print STDERR "Unknown map type specified, 'bc1', 'dh', 'f2', 'cp', and 'gen' are currently accepted.\n";
	usage();
    }

    if ($data_type ne "geno" && $man_cor > 0) {
	print STDERR "You can only specify manual corrections when exporting genotypes.\n";
	usage();
    }

    if ($data_type ne "haplo" && $all_depth > 0) {
	print STDERR "You must use a data type of 'haplo' to export allele depths.\n";
	usage();
    }
}

sub version {
    print STDERR "export_sql.pl ", stacks_version, "\n";
}

sub usage {
    version();

    my $filt;
    my $i = 1;
    foreach my $f (@valid_filters) {
	$filt .= $f . ", ";
	$filt .= "\n        " if ($i % 10 == 0);
	$i++;
    }
    $filt = substr($filt, 0, -11);

    print STDERR 
        "export_sql.pl -D db -b batch_id -a type -f file -o tsv|xls [-m type -c] [-F filter=value ...] [-L lim] [-d] [-h]\n", 
        "    D: database to export from.\n",
        "    b: batch ID of the dataset to export.\n",
	"    a: type of data to export, either 'geno' or 'haplo', for genotypes or observed haplotypes.\n",
        "    f: file to output data.\n",
        "    o: type of data to export: 'tsv' or 'xls'.\n",
	"    d: output depths of alleles instead of the allele values (must use 'haplo' data type).\n",
	"    m: map type. If genotypes are to be exported, specify the map type.\n",
	"    c: include manual corrections if exporting genotypes.\n",
	"  Filters that are applied to select among the catalog loci:\n",
        "      F: one or more filters in the format name=value.\n",
	"         Supported filters: \n",
	"          $filt\n\n",
	"  Filters to be applied to individual sets of haplotype calls (for those selected catalog loci):\n",
	"      A: specify an minimum allele depth limit.\n",
	"      L: specify a minimum locus depth limit.\n",
	"      I: specify a minimum locus log likelihood limit.\n",
        "    h: display this help message.\n\n";
    exit(0);
}
