#!/usr/bin/perl

# Give this script a directory name where a bunch of pre-parsed alanine
# scanning output lives, and it will perform all the dsasa calculations
# that have not yet been performed.
# This script is not strictly necessary, as this calculation is also 
# performed by the final script that detects secondary structure elements
# and creates the SQL/data table output, but it can be useful in a few
# instances:
# 1. you've changed something about the dSASA calculation and you want to
# redo it
# 2. a few files got screwed up (sometimes this can happen if they are
# too big) and you want to handle it using a specialized procedure
# It'll also tell you if there are missing partner files you need to remake!

use strict;
use warnings;

my ( $partner_dir, $prepped_dir ) = @ARGV;

if ( not defined $partner_dir or not defined $prepped_dir ) {
	print "\n";
	print "Run this script by: perl Naccess_all.pl [partner_dir] [prepped_dir] \n\n";
	print "partner_dir is where the 2 chain partner files live, while\n";
	print "prepped_dir is where the parsed output lives \n\n\n";
	exit;
}


open LOG, ">proc.log";

# Naccess path
my $naccess = "/path/to/Naccess/naccess";

my $max_len = 0;
my %chains;
my @this_chains_keys;
my %hotspots;
my $pdb_code;
my $first_chain;
my $second_chain;

system ("ls $prepped_dir > tmp");
open LIST, "tmp";

while (<LIST>) {
    my $ala_out = $_;
    chomp $ala_out;
   
    $pdb_code = $ala_out;
    $pdb_code =~ s/prepped_alaout\///;
    $pdb_code =~ s/_ala_out_prepped\.txt//;
    my $pdb_file = "$partner_dir/" . $pdb_code . ".pdb";

    if ($pdb_code =~ /.{4}_(.)_(.)/) {
        $first_chain  = $1;
        $second_chain = $2;
    }

    unless (-e "$partner_dir/$pdb_code.rsa") {
        print "$pdb_code: Running 2-chain naccess \n";
        system($naccess . " $partner_dir/$pdb_code.pdb");
        system("mv $pdb_code.asa $partner_dir");
        system("mv $pdb_code.log $partner_dir");
        system("mv $pdb_code.rsa $partner_dir");
    }

    foreach my $chn (sort { $a cmp $b } keys %chains) {
        my $pdb4lett = substr($pdb_code, 0, 4);
       
        # Perform Naccess's dSASA calculations on this file unless done already
        unless (-e "$partner_dir/$pdb4lett\_$chn.rsa") {
            if (-e "$partner_dir/$pdb4lett\_$chn.pdb") {
            	print "$pdb_code: Running naccess on the single-chain file.\n";
            	system ($naccess . " $partner_dir/$pdb4lett\_$chn.pdb");
            	system ("mv $pdb4lett\_$chn.asa $partner_dir");
            	system ("mv $pdb4lett\_$chn.log $partner_dir");
           		system ("mv $pdb4lett\_$chn.rsa $partner_dir");
            } else {
                print LOG "Regenerate $pdb4lett\_$chn\n";
                next;
            }
        }        
    }
}

system( "rm tmp" );
