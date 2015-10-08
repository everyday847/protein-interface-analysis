#!/usr/bin/perl
use strict;
use warnings;

# Some scurrilous PDB files have inter-chain disulfides. This is actually quite
# common in symmetric homo-multimers; the issue really is with the PDB concept
# of a chain vs. our idea of a covalently connected entity that you could 
# meaningfully dock or mimic part of or what have you.
# 
# The issue is that alanine scanning involves sampling the bound and unbound 
# state, and the unbound state's score in these complexes is nonsense--because
# they are covalently connected!--and so we have to omit them.

my ( $pdb_list, $pdb_dir ) = @ARGV;

if ( not defined $pdb_list or not defined $pdb_dir ) {
	print "\n";
	print "Run this script by: perl mark_interchain_disulfides.pl [pdb_list] [pdb_directory]\n\n";
	print "\tpdb_list: a text file containing a list of PDB files (xxxx.pdb) found in\n";
	print "\tthe directory pdb_dir. (This allows you to only evaluate a sub-set of all\n";
	print "\tfiles in that directory if desired.)\n\n\n";	 
	exit;
}

open OUT, ">has_interchain_disulfide.txt";

# examine all relaxed PDBs to check on interchain disulfide
open LSPDB, "<$pdb_list";
while (<LSPDB>) {
chomp;
my $pdb = $_;
open PDB, "$pdb_dir/$pdb.pdb" or next;
print $pdb . "\n";
while (<PDB>) {
	if (/^SSBOND/) {
	my $ch1 = substr ($_, 15, 1);
	my $ch2 = substr ($_, 29, 1);
	if ($ch1 ne $ch2) {
		#inter-chain
		if ($ch1 le $ch2) {
			print OUT "$pdb\_$ch1\_$ch2.pdb"; 
			print "$pdb\_$ch1\_$ch2.pdb\n"; 
		}
		else {
			print OUT "$pdb\_$ch2\_$ch1.pdb";
			print "$pdb\_$ch2\_$ch1.pdb\n";
		}
	}
	}
	if (/ATOM/) {
	last;
	}
}

close PDB;



}
close LSPDB;
