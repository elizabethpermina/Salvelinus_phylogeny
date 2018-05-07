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
# Execute the Stacks Pipeline exporter script, wait for it to finish, send
# a notification email to the person who submitted the export.
#
# Written by Julian Catchen <jcatchen@uoregon.edu>
#

use strict;
use File::Temp qw/tempfile/;
use Net::SMTP;

use constant stacks_version => "2.0Beta7c";

#
# Configuration:
#  exe_path:    Path to the export executable.
#  output_path: Path to web-accessible directory to output the export data
#  url:         URL to reach the directory specified by output_path
#  local_host:  Name of localhost to present to SMTP server
#  smtp_host:   Name of SMTP server through which to send mail
#  from:        email address to use in the 'From' field of the message
#
my $exe_path      = "/home/epermina/soft/stacks-2.0Beta7c/bin/" . "export_sql.pl";
my $output_path   = "/home/epermina/soft/stacks-2.0Beta7c/share/stacks/" . "php/export/"; 
my $url           = "http://stackshost.edu/stacks/export/";
my $local_host    = "localhost";
my $smtp_host     = "localhost";
my $from          = "stacks\@stackshost.edu";


my $debug         = 0;
my $db            = "";
my $batch_id      = 0;
my $out_file      = "";      # Generated file name for output data
my $out_type      = "tsv";   # Output format
my $email         = "";      # Email address to send message to
my $data_type     = "haplo"; # Type of export: observed haplotypes or a genetic map.
my $map_type      = "geno";  # If exporting a genetic map, specify map type.
my $depth_lim     = 1;       # Minimum stack depth limit
my $man_cor       = 0;       # Include manual corrections in the export
my $filter_str    = "";      # Comma separated list of filters

parse_command_line();

my (@filters, @cmd_opts, $filter);

#
# Generate file name 
#
my $template = "stacks_export_XXXXXXXX";
my $suffix   = $out_type eq "xls" ? ".xls" : ".tsv";
my (undef, $out_file) = tempfile($template, OPEN => 0, DIR => $output_path, SUFFIX => $suffix);

#
# Prepare the command line parameters
#
push(@cmd_opts, "-D $db", "-b $batch_id", "-f $out_file", "-o $out_type", "-a $data_type");

if ($data_type eq "geno") {
    push(@cmd_opts, "-m $map_type");
}

if ($data_type eq "geno" && $man_cor > 0) {
    push(@cmd_opts, "-c");
}

if ($data_type eq "haplo" && $depth_lim > 1) {
    push(@cmd_opts, "-L $depth_lim");
}

if (length($filter_str) > 0) {
    @filters = split(/,/, $filter_str);

    foreach $filter (@filters) {
        push(@cmd_opts, "-F $filter");
    }
}

my $cmd = join(" ", @cmd_opts);
$cmd = $exe_path . " " . $cmd;
print STDERR "CMD: $cmd\n" if ($debug);

#
# Execute the exporter program
#
my @results = `$cmd`;
#my @results = `echo Success`;

#
# Check the results, we expext a one line result: either 'Success' or 'Failure'
#
chomp $results[0];

if ($results[0] eq "Success") {
    send_email('success', $out_file);
} else {
    send_email('failure', $out_file);
}


sub send_email {
    my ($result) = @_;

    my $smtp = Net::SMTP->new($smtp_host,
                              'Hello' => $local_host,
                              'Timeout' => 60);
    $smtp->mail($from);
    $smtp->recipient($email);

    my $msg .= 
        "From: $from\r\n" .
        "To: $email\r\n" .
        "Subject: Stacks pipeline export complete\r\n" .
        "\r\n";

    if ($result eq "success") {
        #
        # Trim the path off the output file
        #
        my ($f) = ($out_file =~ /.*\/(stacks_export_\w{8}\.\w{3})$/);

        $msg .= "Your data has been exported and can be downloaded from: " . $url . $f . "\r\n";
    } else {
        
        $msg .= "There has been an error exporting your data, please contact the system administrator.\r\n";
    }

    $smtp->data($msg);

    $smtp->quit();
}

sub parse_command_line {
    while (@ARGV) {
	$_ = shift @ARGV;
	if    ($_ =~ /^-d$/) { $debug++; }
	elsif ($_ =~ /^-D$/) { $db         = shift @ARGV; }
	elsif ($_ =~ /^-b$/) { $batch_id   = shift @ARGV; }
	elsif ($_ =~ /^-e$/) { $email      = shift @ARGV; }
        elsif ($_ =~ /^-a$/) { $data_type  = lc(shift @ARGV); }
        elsif ($_ =~ /^-m$/) { $map_type   = lc(shift @ARGV); }
        elsif ($_ =~ /^-L$/) { $depth_lim  = shift @ARGV; }
        elsif ($_ =~ /^-c$/) { $man_cor++; }
	elsif ($_ =~ /^-t$/) { $out_type   = shift @ARGV; }
	elsif ($_ =~ /^-F$/) { $filter_str = shift @ARGV; }
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

    if (length($email) == 0) {
        print STDERR "You must specify an email for notification.\n";
        usage();
    }

    if ($out_type ne "tsv" && $out_type ne "xls") {
        print STDERR "Output type can only be 'tsv' or 'xls'.\n";
        usage();
    }

    if ($batch_id !~ /^\d{1,4}$/) {
        print STDERR "Batch ID must be a numeric value.\n";
        usage();
    }
}

sub version {
    print STDERR "stacks_export_notify.pl ", stacks_version, "\n";
}

sub usage {
    version();

    print << "EOQ";
stacks_export_notify.pl -e email -D db -b batch_id [-t type] [-F filters] [-d] [-h]
  e: email to use for notification.
  D: radtag database to examine.
  b: batch_id of data set to export.
  a: type of data to export, either 'gen' or 'haplo', for genotypes or observed haplotypes.
  t: output type, either 'tsv' or 'xls'.
  F: comma separated list of filters to apply to the data.
  L: if exporting observed haplotypes, specify a stack depth limit.
  m: map type. If genotypes are to be exported, specify the map type.
  c: include manual corrections if exporting genotypes.
  h: display this help message.
  d: turn on debug output.

EOQ

    exit(0);
}
