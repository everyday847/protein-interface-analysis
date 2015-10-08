#!/usr/bin/perl

use strict;
use warnings;

# Some PDBs contain selenomethionine residues (generally not biologically
# significant; just helps crystallization); these should be changed to 
# regular methionine. This script curates a list thereof for further processing.

my ( $pdb_list, $pdb_dir ) = @ARGV;

if ( not defined $pdb_list or not defined $pdb_dir ) {
	print "\n";
	print "Run this script by: perl check_mse.pl [pdb_list] [pdb_directory]\n\n";
	print "\tpdb_list: a text file containing a list of PDB files (xxxx.pdb) found in\n";
	print "\tthe directory pdb_dir. (This allows you to only evaluate a sub-set of all\n";
	print "\tfiles in that directory if desired.)\n\n\n";	 
	exit;
}

open MSE, ">has_mse.txt";
open LIST, "<$pdb_list";

while (<LIST>) {
	my $pdb = $_;
	chomp $pdb;
	print $pdb . "\n";

	# Generate filename from PDB code; open it.
	my $pdb_file = "$pdb_dir/$pdb\.pdb";
	
	open(PDB, "<$pdb_file")
	  or next;#die "ERROR: Could not open PDB file ($pdb).\n";
	my $hadmse = 0;
	foreach (<PDB>) {
		my $line = $_;
		
		if ($line =~ /^HETATM/ and $line =~ /MSE/ ) {
			# rescue MSE
			$hadmse = 1; 
			last;
		}
	        
	}
	if ($hadmse == 1) {
		print MSE "$pdb\n";
	}	
	
}
