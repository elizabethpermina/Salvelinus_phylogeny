#!/usr/bin/env perl
#
# Copyright 2011, Julian Catchen <jcatchen@uoregon.edu>
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

use strict;
use constant stacks_version => "2.0Beta7c";

my $debug        = 0;
my $paired       = 0;
my $amos         = 0;
my $ins_len_dist = 0;
my $single_path  = "";
my $paired_path  = "";
my $sanger_path  = "";
my $out_path     = ".";
my $white_list   = "";
my $hash_len     = 27;
my $insert_len   = 0;
my $exp_cov      = 0;
my $cov_cut      = 0.0;
my $min_len      = 0;
my $read_trk     = 0;
my $clean        = 1;
my $collate      = 0;
my $exe_path     = "";
my $velveth      = "velveth";
my $velvetg      = "velvetg";

parse_command_line();

if (length($exe_path) > 0) {
    $velveth = $exe_path . "/" . $velveth;
    $velvetg = $exe_path . "/" . $velvetg;
}

#
# Test that we can execute the velvet programs
#
die ("Unable to find '" . $velveth . "'.\n") if (!-e $velveth || !-x $velveth);
die ("Unable to find '" . $velvetg . "'.\n") if (!-e $velvetg || !-x $velvetg);

my (@locus_files, $num_files, $file, $input_file, $output_file, $hres_file, $gres_file, $collate_fh);

build_file_list(\@locus_files);
$num_files = scalar(@locus_files);

if ($collate) {
    open($collate_fh, ">$out_path/collated.fa") or die("Unable to open collate file, $!\n");
}

my ($sing_data, $pair_data, $sang_data, $ins, $cov, $afg, $cut, $min, $read, $cln);

$ins   = $paired   > 0 ? "-ins_length $insert_len"   : "-ins_length auto";
$cov   = $paired   > 0 ? "-exp_cov $exp_cov"         : "-exp_cov auto";
$cut   = $cov_cut  > 0 ? "-cov_cutoff $cov_cut"      : "-cov_cutoff auto";
$read  = $read_trk > 0 ? "-read_trkg yes"            : "";
$cln   = $clean    > 0 ? "-very_clean yes"           : "";
#$min   = $min_len  > 0 ? "-min_contig_lgth $min_len" : ""

#
# Write out the parameters for this assembly
#
open(PARAM, "> $out_path/velvet_parameters.txt") or die("Unable to open parameter file: $!\n");
print PARAM
    "Single-end Path: ", $single_path, "\n",
    "Paired-end Path: ", $paired_path, "\n",
    "Sanger Path:     ", $sanger_path, "\n",
    "Hash Length:     ", $hash_len, "\n",
    "Insert Length:   ", $ins, "\n",
    "Coverage:        ", $cov, "\n",
    "Coverage Cutoff: ", $cut, "\n",
    "Miniumum contig length: ", $min, "\n",
    "Read tracking:   ", $read, "\n",
    "Very Clean:      ", $cln, "\n"; 
close(PARAM);

my $i = 1;

foreach $file (@locus_files) {

    ($file) = ($file =~ /(.+)\.fa/); 
    $output_file = $out_path . "/" . $file;
    $hres_file   = $out_path . "/" . $file . "-h.output";
    $gres_file   = $out_path . "/" . $file . "-g.output";

    $sing_data = length($single_path) > 0 ? '-short '       . $single_path . "/" . $file . ".fa" : "";
    $pair_data = length($paired_path) > 0 ? '-shortPaired ' . $paired_path . "/" . $file . ".fa" : "";
    $sang_data = length($sanger_path) > 0 ? '-long '        . $sanger_path . "/" . $file . ".fa" : "";

    print STDERR "Assembling locus '$file'; run $i of $num_files        \r";

    # Execute velveth to build hash table, then velvetg to assemble
    print STDERR "$velveth $output_file $hash_len -fasta $sing_data $pair_data $sang_data &> $hres_file\n" if ($debug);
    `$velveth $output_file $hash_len -fasta $sing_data $pair_data $sang_data &> $hres_file`;

    print STDERR "$velvetg $output_file $ins $cov $cut $min $read $cln &> $gres_file\n" if ($debug);
    `$velvetg $output_file $ins $cov $cut $min $read $cln &> $gres_file`;

    collate_and_clean($out_path, $file, $collate_fh) if ($collate);

    $i++;
}

close($collate_fh) if ($collate);

sub collate_and_clean {
    my ($out_path, $file, $collate_fh) = @_;

    my (@seqs, $seq);

    parse_fasta("$out_path/$file/contigs.fa", \@seqs);

    foreach $seq (@seqs) {
	next if (length($seq->{'seq'}) < $min_len);

	$seq->{'id'} = $file . "|" . $seq->{'id'};

	print_fasta($collate_fh, $seq);
    }

    `rm $out_path/$file-g.output`;
    `rm $out_path/$file-h.output`;
    `rm -r $out_path/$file`;
}

sub parse_fasta {
    my ($file, $seqs) = @_;

    my ($fh, $line, $buf, $id, $seq);

    open($fh, "<$file") 
	or die("Unable to open Velvet output file: $file, $!\n");

    while ($line = <$fh>) {
	chomp $line;

	if (substr($line, 0, 1) eq ">") {
	    if (length($buf) > 0) {
		$seq = {};
		$seq->{'id'}  = $id;
		$seq->{'seq'} = $buf;

		push(@{$seqs}, $seq);
		$buf = "";
	    }
	    $id = substr($line, 1);

	} else {
	    $buf .= $line;
	}
    }

    if (length($buf) > 0 && length($id) > 0) {
	$seq = {};
	$seq->{'id'}  = $id;
	$seq->{'seq'} = $buf;
	push(@{$seqs}, $seq);
    }

    close($fh);
}

sub print_fasta {
    my ($fh, $seq) = @_;

    my ($s);

    print $fh ">", $seq->{'id'}, "\n";

    $s = $seq->{'seq'};

    while (length($s) > 60) {
	print $fh substr($s, 0, 60), "\n";
	$s = substr($s, 60);
    }

    print $fh $s, "\n" if (length($s) > 0);
}

sub build_file_list {
    my ($files) = @_;

    my (@ls, $line, $file, $path);

    # Load a white list of files to process if it is supplied.
    my @wl;
    if (length($white_list) > 0) {
	load_white_list(\@wl);
    }

    $path = length($paired_path) > 0 ? $paired_path : $single_path;

    @ls = `ls -1 $path/`;

    foreach $line (@ls) {
	chomp $line;

	next if (length($line) == 0);
        next if ($line !~ /.+\.fa$/ && $line !~ /.+\.fasta$/);	

	($file) = ($line =~ /^(.+\.fas?t?a?)/); 

        if (scalar(@wl) > 0) {
	    next if (!grep(/^$file$/, @wl));
	}

	push(@{$files}, $file);
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
	if    ($_ =~ /^-s$/) { $single_path = shift @ARGV; }
	elsif ($_ =~ /^-p$/) { $paired_path = shift @ARGV; }
	elsif ($_ =~ /^-l$/) { $sanger_path = shift @ARGV; }
	elsif ($_ =~ /^-o$/) { $out_path    = shift @ARGV; }
	elsif ($_ =~ /^-c$/) { $collate++; }
	elsif ($_ =~ /^-W$/) { $white_list  = shift @ARGV; }
	elsif ($_ =~ /^-I$/) { $insert_len  = shift @ARGV; }
	elsif ($_ =~ /^-C$/) { $exp_cov     = shift @ARGV; }
	elsif ($_ =~ /^-T$/) { $cov_cut     = shift @ARGV; }
	elsif ($_ =~ /^-R$/) { $read_trk    = shift @ARGV; }
	elsif ($_ =~ /^-M$/) { $min_len     = shift @ARGV; }
	elsif ($_ =~ /^-H$/) { $hash_len    = shift @ARGV; }
	elsif ($_ =~ /^-P$/) { $paired++; }
	elsif ($_ =~ /^-L$/) { $clean    = 0; }
	elsif ($_ =~ /^-e$/) { $exe_path = shift @ARGV; }
	elsif ($_ =~ /^-v$/) { version(); exit(); }
	elsif ($_ =~ /^-d$/) { $debug++; }
	elsif ($_ =~ /^-h$/) { usage(); }
	else {
	    print STDERR "Unknown command line option '$_'\n";
	    usage();
	}
    }

    $single_path = substr($single_path,  0, -1) if (substr($single_path, -1)  eq "/");
    $paired_path = substr($paired_path,  0, -1) if (substr($paired_path, -1)  eq "/");
    $out_path    = substr($out_path, 0, -1)     if (substr($out_path, -1)     eq "/");
    $exe_path    = substr($exe_path, 0, -1)     if (substr($exe_path, -1)     eq "/");
}

sub version {
    print STDERR "exec_velvet.pl ", stacks_version, "\n";
}

sub usage {
    version();

    print STDERR <<EOQ; 
exec_velvet.pl -p path -s path [-l path] -o path [-c] [-H len] [-P] [-I len] [-C exp_cov] [-T cov_cut]
              [-W file_white_list] [-w marker_white_list] [-L] [-e path] [-d] [-h]
  p: path to the paired-end FASTA files to assemble.
  s: path to the single-end FASTA files to assemble.
  l: path to long, sanger-style reads to assemble.
  o: path to output the assembled files.
  c: collate the resulting velvet runs into a single FASTA file, clean velvet directories.
  W: a white list of files to process in the input path.
  H: length of overlap required for reads (hash length, default 27).
  P: process paired-end reads.
  I: insert length (for paired-end reads; see velvet documentation).
  C: expected coverage (for paired-end reads).
  T: coverage cutoff.
  M: minimum contig length, discard contigs shorter than this value.
  R: turn on velvet's read tracking (uses additional memory).
  L: leave velvet's intermediate files behind.
  e: executable path, location of velvet programs.
  h: display this help message.
  d: turn on debug output.

EOQ

  exit(0);
}
