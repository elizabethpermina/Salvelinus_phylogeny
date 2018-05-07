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
use DBI;

use constant stacks_version => "2.0Beta7c";

my $debug        = 0;
my $mysql_config = "/home/epermina/soft/stacks-2.0Beta7c/share/stacks/" . "sql/mysql.cnf";
my $db           = "";
my $batch_id     = 0;
my $in_path      = "";
my $seq_path     = "";
my $seq_type     = "";

parse_command_line();

my (%sth);

prepare_sql_handles(\%sth);

if ($seq_type eq "est") {
    load_ests(\%sth);
} elsif ($seq_type eq "pe_radtag") {
    load_pe_radtags(\%sth);
} else {
    print STDERR "Unknown sequence type specified '$seq_type'\n";
}

sub load_pe_radtags {
    my ($sth) = @_;

    #
    # Open the file containing sequences and load them into the database. 
    #
    my ($fh, $line, $id, $buf, $catalog_id, $seq_id, $count);

    open($fh, "<$seq_path") 
	or die("Unable to open FASTA input file: '$seq_path', $!\n");

    $count = 0;
    while ($line = <$fh>) {
	chomp $line;

	if (substr($line, 0, 1) eq ">") {
	    if (length($buf) > 0) {

		($catalog_id, $seq_id) = split(/\|/, $id);
		$sth->{'load'}->execute($batch_id, $catalog_id, $seq_type, $seq_id, $buf) 
		    or die("Unable to insert results into $db.\n");
		$count++;

		$buf = "";
	    }
	    $id  = substr($line, 1);

	} else {
	    $buf .= $line;
	}
    }

    if (length($buf) > 0 && length($id) > 0) {

	($catalog_id, $seq_id) = split(/\|/, $id);
	$sth->{'load'}->execute($batch_id, $catalog_id, $seq_type, $seq_id, $buf) 
	    or die("Unable to insert results into $db.\n");
	$count++;
    }

    print STDERR "Loaded $count RAD-Tag paired-end contigs into the database.\n";

    close($fh);
    close_sql_handles(\%sth);
}

sub load_ests {
    my ($sth) = @_;

    #
    # Load a TSV file of BLAST hits linking sequences to markers in the catalog.
    #
    my %hits;

    load_blast_hits(\%hits);

    my ($fh, $line, $id, $buf, $catalog_id, $count);

    #
    # Open the file containing sequences. For each sequence that has a BLAST
    # hit to a marker in the catalog, load the sequence into the database. 
    #
    open($fh, "<$seq_path") 
	or die("Unable to open FASTA input file: '$seq_path', $!\n");

    $count = 0;
    while ($line = <$fh>) {
	chomp $line;

	if (substr($line, 0, 1) eq ">") {
	    if (length($buf) > 0) {

		if (defined($hits{$id})) {
		    foreach $catalog_id (@{$hits{$id}}) {
			$sth->{'load'}->execute($batch_id, $catalog_id, $seq_type, $id, $buf) 
			    or die("Unable to insert results into $db.\n");
			$count++;
		    }
		}
		$buf = "";
	    }
	    $id  = substr($line, 1);

	} else {
	    $buf .= $line;
	}
    }

    if (length($buf) > 0 && length($id) > 0) {
	if (defined($hits{$id})) {
	    foreach $catalog_id (@{$hits{$id}}) {
		$sth->{'load'}->execute($batch_id, $catalog_id, $seq_type, $id, $buf) 
		    or die("Unable to insert results into $db.\n");
		$count++;
	    }
	}
    }

    print STDERR "Loaded $count EST sequences into the database.\n";

    close($fh);
    close_sql_handles(\%sth);
}

sub load_blast_hits {
    my ($hits) = @_;

    my (@parts, $line);

    open(HITS, "<" . $in_path) or die("Unable to open BLAST hits file '$in_path' $!\n");

    #
    # First line is a header.
    #
    $line = <HITS>;

    while ($line = <HITS>) {
        chomp $line;
        @parts = split(/\t/, $line);

        if (!defined($hits->{$parts[7]})) {
            $hits->{$parts[7]} = [];
        }

        if (!grep(/^$parts[2]$/, @{$hits->{$parts[7]}})) {
            push(@{$hits->{$parts[7]}}, $parts[2]);
        }
    }

    print STDERR "Loaded ", scalar(keys %{$hits}), " distinct BLAST hits.\n";

    close(HITS);
}

sub prepare_sql_handles {
    my ($sth) = @_;

    my $cnf = (defined($ENV{"HOME"}) && -e $ENV{"HOME"} . "/.my.cnf") ? $ENV{"HOME"} . "/.my.cnf" : $mysql_config;
    $sth->{'dbh'} = DBI->connect("DBI:mysql:$db:mysql_read_default_file=$cnf")
	or die("Unable to connect to the $db MySQL Database!\n" . $DBI::errstr);

    my $query;

    $query = 
        "INSERT INTO sequence SET batch_id=?, catalog_id=?, type=?, seq_id=?, seq=?";
    $sth->{'load'} = $sth->{'dbh'}->prepare($query) or die($sth->{'dbh'}->errstr());
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
	if    ($_ =~ /^-p$/) { $in_path   = shift @ARGV; }
	elsif ($_ =~ /^-f$/) { $seq_path  = shift @ARGV; }
	elsif ($_ =~ /^-b$/) { $batch_id  = shift @ARGV; }
	elsif ($_ =~ /^-D$/) { $db        = shift @ARGV; }
	elsif ($_ =~ /^-t$/) { $seq_type  = shift @ARGV; }
	elsif ($_ =~ /^-d$/) { $debug++; }
	elsif ($_ =~ /^-h$/) { usage(); }
	else {
	    print STDERR "Unknown command line option: '$_'\n";
	    usage();
	}
    }

    if (length($seq_type) == 0) {
        print STDERR "You must specify the sequence type.\n";
        usage();
    }

    if ($batch_id == 0) {
        print STDERR "You must specify the batch ID.\n";
        usage();
    }
}

sub version {
    print STDERR "load_sequences.pl ", stacks_version, "\n";
}

sub usage {
    version();

    print STDERR <<EOQ; 
load_sequences.pl -b id -D db -f path -t type [-p blast_path] [-d] [-h]
    b: batch ID.
    D: database to examine.
    f: path to FASTA file containing sequences to load.
    t: sequence type to load: either 'pe_radtag' or 'est'.
    p: path to tab-separated file containing BLAST hits.
    h: display this help message.
    d: turn on debug output.

EOQ

    exit(0);
}
