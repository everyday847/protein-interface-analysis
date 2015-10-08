#!/usr/bin/perl

# ------------------------------------------------------------------------------
# --- Use curl to obtain PDB file for PDB ID -----------------------------------
# ------------------------------------------------------------------------------

use strict;
use warnings;

if ( not defined $pdb_list or not defined $pdb_dir ) {
	print "\n";
	print "Run this script by: perl curl.pl [pdb_list] [pdb_directory] \n\n";
	print "\tpdb_list: a text file containing a list of PDB codes (xxxx) to be added ton\n";
	print "\tthe directory pdb_dir.";

open IN, "<$pdb_list";
my @pdblines = <IN>;

foreach my $pdb (@pdblines) {
	chomp $pdb;

	# If the PDB is present, then check its first line
	# for a common download error
 	if (-e "$pdb_dir/$pdb.pdb") {
		open PDB, "<$pdb_dir/$pdb.pdb";
		my $l = <PDB>;
		close PDB;
		if ($l =~ /<!DOCTYPE/) {
			system ("rm $pdb_dir/$pdb.pdb");
			print "$pdb\n";
		} else {
			next;
		}
	}

	# Try to download biological assembly (.pdb1 extension on rcsb website)
	system("curl -o $pdb_dir/$pdb.pdb \"http://www.rcsb.org/pdb/files/$pdb.pdb1\" ";

	# Check if that failed--if so, only asymmetric unit is available.
	open PDB, "<$pdb_dir/$pdb.pdb";
	my $line = <PDB>;
	close PDB;
	if ($line =~ /^<html>/) {
		system ("rm $pdb_dir/$pdb.pdb");
		system ("curl -o $pdb_dir/$pdb.pdb \"http://www.rcsb.org/pdb/files/$pdb.pdb\" ");
	}
}
