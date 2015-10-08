#!/usr/bin/perl

use strict;
use warnings;

my ( $pdb_list, $pdb_dir, $disulfide_dir ) = @ARGV;

if ( not defined $pdb_list or not defined $pdb_dir or not defined $disulfide_dir ) {
	print "\n";
	print "Run this script by: perl sequester_interchain_disulfides.pl [pdb_list] [pdb_directory] [disulfide_dir]\n\n";
	print "\tpdb_list: a text file containing a list of PDB files (xxxx.pdb) found in\n";
	print "\tthe directory pdb_dir. (This allows you to only evaluate a sub-set of all\n";
	print "\tfiles in that directory if desired.)\n"
	print "\n";
	print "\tThen, move all such files to a special indicated directory--don't process them\n";
	print "\tfurther, but maybe you'll want to do something special with them (like not move";
	print "\tthe chains apart in alascan or something)\n\n\n";	 
	exit;
}

open CODES, "<$pdb_list";

while (<CODES>) {
	chomp;
	my $c = $_;
	$c =~ s/.pdb//g;
	system ("mv $pdb_dir/$c* $disulfide_dir/");
}
close CODES;
