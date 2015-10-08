#!/usr/bin/perl

use strict;
use warnings;

# Some PDB files use A/B/C next to a residue name (ALYS, BLYS)  to indicate
# alternate conformations. These conformations do not need to be retained, 
# because only one can be alanine scanned (and because we are going to be
# refining the sidechain conformations anyway).

# We also can't handle residues that are numbered as e.g. 81A due to insertion.
# We don't choose to renumber them because then the numbers in HIPP would
# fail to match the numbers if users download the PDB themselves. We just omit
# them, since they only have a tenuous relationship with the native protein
# sequence anyway and because they're typically on loops anyway.

my ( $pdb_list, $pdb_dir ) = @ARGV;

if ( not defined $pdb_list or not defined $pdb_dir ) {
	print "\n";
	print "Run this script by: perl check_for_numbering.pl [pdb_list] [pdb_directory]\n\n";
	print "\tpdb_list: a text file containing a list of PDB files (xxxx.pdb) found in\n";
	print "\tthe directory pdb_dir. (This allows you to only evaluate a sub-set of all\n";
	print "\tfiles in that directory if desired.)\n\n\n";	 
	exit;
}

open LIST, "<$pdb_list";

while (<LIST>) {
	chomp;
	my $pdb = $_;
	open PDB, "$pdb_dir/$pdb" or die("died on $pdb\n");
	while (<PDB>) {
		my $line = $_;
		chomp $line;
		if ($line =~ /^ATOM/) {
			my $ins = substr $line, 16, 1;
			my $conf = substr $line, 26, 1;
			my $resn = substr $line, 22, 4;
			if ($ins != " ") {
				print "$pdb insertion at $resn";
			}
			if ($conf != " ") {
				print "$pdb alternate conf at $resn";
			}
		}
	}
	close PDB;
}
