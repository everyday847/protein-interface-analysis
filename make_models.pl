#!/usr/bin/perl

# First, make sure the file 
# This is the first step in processing a PDB: making individual model files that
# you subsequently will relax (prepare for further analysis). The reason this is
# necessary is because NMR structures contain multiple starting points, and we
# want to incorporate all of that information. Ditto crystal structures whose
# asymmetric units contain multiple instances of the biological interface of 
# interest. Take, for example, 1JPW. If you take the asymmetric unit, There is 
# one interface of interest--TCF4-beta catenin--but three sources of data on it.
# Our scripts obtain the biological assembly for the PDB, which represents this
# as multiple separate models.
#
# This script would take 1JPW and create 1JPW_m1, 1JPW_m2, 1JPW_m3.

use strict;
use warnings;

my ( $pdb_list, $pdb_dir, $relax_dir ) = @ARGV;

if ( not defined $pdb_list or not defined $pdb_dir or not defined $relax_dir ) {
	print "\n";
	print "Run this script by: perl check_for_numbering.pl [pdb_list] [pdb_directory] [relax_dir]\n\n";
	print "\tpdb_list: a text file containing a list of PDB files (xxxx.pdb) found in\n";
	print "\tthe directory pdb_dir. (This allows you to only evaluate a sub-set of all\n";
	print "\tfiles in that directory if desired.)\n"
	print "\n";
	print "\trelax_dir: the directory where your refined PDB structures will live, and \n";
	print "\twhere sub-directories for each PDB file's multiple starting models should go.\n\n\n";	 
	exit;
}

open LIST, "<$pdb_list"

while (<LIST>) {
	my $pdb = $_;
	chomp $pdb;
	
	my $model_count = 1;
	print $pdb . "\n";
	
	# Continue on if we have already finished relaxing this PDB (i.e. the NEXT step)
	# and delete the directory that used to hold the models, if necessary
	if (-e "$relax_dir/$pdb.pdb") {
		if (-d "$relax_dir/$pdb\_models") {
			system ("rm -rf pdb_files_new/$pdb\_models");
		}
		next;
	}
	
	# Generate filename from PDB code; open it.
	my $pdb_file = "$pdb_dir/$pdb\.pdb";
	
	open PDB, "<$pdb_file" or next;
	
	if (!-d "$relax_dir/$pdb\_models") { system ("mkdir pdb_files_new/$pdb" . "_models")};
	open MDLOUT, ">$relax_dir/$pdb" . "_models/$pdb" . "_m$model_count.pdb";

	# We can't handle files that only contain CA atoms (these are CA traces from very
	# low resolution structures and require special, challenging structure preparation
	# to use. Totally out of scope! So we keep track of the number of atoms to 
	# eliminate these from consideration.
	my $num_res = 0;
	my $num_atoms = 0;	
	my $last_rn = "XXXX";
	my $remove = 0;
	my $any_but_unk = 0;
	foreach (<PDB>) {
		my $line = $_;
		
		if (/^ATOM/) {
			$num_atoms++; 
			if ((substr ($line, 22, 4)) ne $last_rn ) {
				$last_rn = substr ($line, 22, 4);
				$num_res++;
			}
		}

		# We can't do much with UNK records, but we can at least alascan the other parts
		# of the file. So flag if we ever find a non-UNK residue in this PDB.
		if (substr ($line, 17, 3) eq "UNK") {next;} else {$any_but_unk = 1;}
		
		# Skip alternate sidechain conformations
		if (substr ($line, 16, 1) =~ /[B-Z]/) {next;}

		(substr ($line, 16, 1)) = " ";
	        
		# We are done with most possible substitutions.

		if ($line =~ /^TER/ or $line =~ /^ATOM/) {
			print MDLOUT $line;
		}

		# Replace selenomethionine with methionine.
		if ($line =~/^HETATM/ and $line =~ /MSE/) {
			$line =~ s/HETATM/ATOM  /g;
			$line =~ s/MSE/MET/g;
			$line =~ s/SE/ S/g;
			print MDLOUT $line;
		}
	    
		# Open next model file. Flag this one for removal if it was a CA trace.
		if ($line =~ /^ENDMDL/) {
			close MDLOUT;
			#print "nr $num_res na $num_atoms\n";
			if ($num_atoms <= $num_res + 2 ) {$remove = 1;}

			$num_res = 0;
			$num_atoms = 0;
			$last_rn = "XXXX";
			$model_count++;
			open MDLOUT, ">$relax_dir/$pdb" . "_models/$pdb" . "_m$model_count.pdb";
	    }
	}

	# Remove if CA trace or if all UNK records
	if ($remove == 1 || $any_but_unk == 0) {
		print "REMOVE $pdb\n";
		system ("rm -rf $relax_dir/$pdb\_models/"); 
		system ("rm -rf $pdb_dir/$pdb.pdb");
	}
	
	close MDLOUT;
	close PDB;
	
}
close LIST;
