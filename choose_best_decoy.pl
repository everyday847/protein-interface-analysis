#!/usr/bin/perl

# This scripts presumes that you have done the following:
# 1. Obtained some number of PDB files and parsed them into separate models.
# 2. Used some form of sampling on those model files to generate 1-999
# "decoy" structures per model. (Preferably, for the survey type of study
# we are doing, where we are preparing many different structures in the same
# way, just create 1-5. If you have a system you really care about,
# do more sampling--email me.)
#
# Now, you want go back down to one best structure to start from--which you
# will split into two-chain files and alanine scan. (For example, for 1NVW,
# we would create 1NVW_Q_R, 1NVW_Q_S, 1NVW_R_S.)
#
# Do this by examining the scorefile, finding the tag for the structure
# with the lowest score, and moving it into the desired directory.
# BTW, make sure to auto-detect if PDBs are finished or not so you don't
# just grab a PDB with one finished structure and say it's the best.

use strict;
use warnings;

my ( $pdb_list, $pdb_dir, $relax_dir ) = @ARGV;

if ( not defined $pdb_list or not defined $pdb_dir or not defined $relax_dir ) {
	print "Usage: perl choose_best_decoy.pl [pdb_list] [pdb_dir] [relax_dir]\n\n\n";
}

open LIST, "<$pdb_list";
open REMAINING, ">remaining.txt";
while (<LIST>) {
	chomp;
	if (!-e "$relax_dir/$_.pdb") {print REMAINING "$_\n"; print "REM $_\n"; }
}
close REMAINING;

# Assume that your output files have been accumulating in this directory.
# Find unique associated PDB codes.
system ("ls *m*.pdb|sed s'/_m[0-9]*_....\.pdb//'|uniq > codes");
open CODES, "codes";

open MODELLEN, ">model_lengths";
open REMAINING, "remaining.txt";

foreach my $line (<REMAINING>) {
	print $line;
	chomp $line;
	open PDB, "$pdb_dir/$line.pdb";
	my $max = 0;
	while (<PDB>) {
		if (/^MODEL/) {
			$max++;
		}
	}
	close PDB;
	print MODELLEN "$line\t$max\n";
}
close MODELLEN;

open MODELLEN, "model_lengths";
my %max;
while (<MODELLEN>) {
	my @objs = split;
	$max{$objs[0]} = $objs[1];
	if ($max{$objs[0]} == 0) {
		$max{$objs[0]} = 1;
	}
	#print $max{$objs[0]} . " " . $objs[0] . "\n";
}
close MODELLEN;

while (<CODES>) {
	chomp;
	my $pdbc = $_;
	#print "|$pdbc|\n";
	if (!defined $max{$pdbc}) {
		print "Issue with $pdbc\n";
		next CODELOOP;
	}
	for (my $num = 1; $num <= $max{$pdbc}; $num++) {
		for (my $nn = 1; $nn <= 1; $nn++) {
			if (!-e $pdbc . "_m$num" . "_000$nn.pdb") {print "Incomplete $pdbc" . "_m$num" . "_000$nn.pdb!\n";next CODELOOP;}
		}
	}

	#all files exist
	my $pdb_final;
	my @files = <$pdbc*.pdb>;
	my %scores;

	foreach my $file (@files) {
		#print $file . "\n";
		open FILE, $file;
		my $num = 19;
		while (<FILE>) {
			if (/^label/) {
				my @lineobjs = split;
				for (my $i = 0; $i < (scalar(@lineobjs)); ++$i ) {
					if ($lineobjs[$i] eq "total") {
						$num = $i;
					}
				}
			}
			if (/^pose/) {
				my @lineobjs = split;
				$scores{$file} = $lineobjs[$num];
				last;
			}
		}
	
		close (FILE);
	}
	
	foreach my $file (sort {$scores{$a} <=> $scores{$b}} keys %scores) {
		$pdb_final = $file; #print $pdb_final . "\n";
		last;
	}
	
	print $pdb_final . "\n";
	if (-e "$relax_dir") {
	} else {
		system ("mkdir $relax_dir");
	}
	system ("cp $pdb_final $relax_dir/$pdbc.pdb");
	# Store number of chains, and their IDs in a hash
	my $nchain = 0;
	my %chain_ids;
		   
	print STDERR "Processing PDB file $pdb_final...\n";
	open PDB, $pdb_final;   
	foreach (<PDB>) {
		if (/^ATOM/) {
			my $letter = substr($_,21,1);
			$chain_ids{$letter} = 1;
		}
		$nchain++ if (/^TER/);
	}
	
	close PDB;
	
	
	# Create, per PDB file, each of the possible n(n-1)/2 two-chain files
		   
	# Make and sort @all_chains, so we don't look for B_A instead of A_B
	my @all_chains = sort(keys %chain_ids);
		 
	if (!-e "pdb_partner_files") {
			system ("mkdir ./pdb_partner_files");
	}
	
	# For all possible interchain partners...
	for (my $chain_index_1 = 0; $chain_index_1 < scalar (@all_chains); $chain_index_1++) {
		for (my $chain_index_2 = $chain_index_1 + 1; $chain_index_2 < scalar (@all_chains); $chain_index_2++) {
			my $partner_1 = $all_chains[$chain_index_1];
			my $partner_2 = $all_chains[$chain_index_2];
	
			# Deduce output PDB filename
			my $out_file = "pdb_partner_files/" . $pdbc;
			$out_file .= "_" . $all_chains[$chain_index_1];
			$out_file .= "_" . $all_chains[$chain_index_2];
			$out_file .= ".pdb";
				
			# Skip writing this PDB partner file if it already exists
			if (-e $out_file) {
				next;
			}
				
			# Create file and write to it
			open OUT_PDB, ">$out_file";

			# Print coordinates for first partner to file
			# Look for lines that start ATOM and contain the partner_1 chain ID
			# or partner_2 chain ID.
				
			my @partner_1_lines;
			my @partner_2_lines;
				
			open(PDB, "<$pdb_final")
				or die "ERROR: Could not open PDB file ($pdb_final).\n";
 
			foreach (<PDB>) {
				if (/^ATOM/) {
					#my @info = split;
					#if ($info[0] eq "ATOM") {
					my $cchn = substr $_, 21, 1;
					if ($cchn eq $partner_1) {
						push (@partner_1_lines, $_);
					} elsif ($cchn eq $partner_2) {
						push (@partner_2_lines, $_);
					}
				}
			}

			close PDB;

			# Print these lines to the output file
			# Print TER at end of each chain
			print OUT_PDB @partner_1_lines;
			print OUT_PDB "TER \n";
			print OUT_PDB @partner_2_lines;
			print OUT_PDB "TER \n";
			close OUT_PDB;
		}
	}
		# Finally, print coordinates of chains in individual files, one per chain.
	   
		# So if PDB file 1A0H has chains A, B, D, and E. The files created above will be:
		#	   1A0H_A_B.pdb			1A0H_B_D.pdb
		#	   1A0H_A_D.pdb			1A0H_B_E.pdb
		#	   1A0H_A_E.pdb			1A0H_D_E.pdb
	   
		# This step creates the following individual files:
		#	   1A0H_A.pdb
		#	   1A0H_B.pdb
		#	   1A0H_D.pdb
		#	   1A0H_E.pdb
	   
	foreach my $chain (@all_chains) {
		my $ind_filename = "pdb_partner_files/" . $pdbc . "_" . $chain . ".pdb";
		if (-e $ind_filename) {
			next;
		}

		open IND_FILE, ">$ind_filename"
			or die "ERROR: Could not open individual PDB chain file ($ind_filename).\n";
			   
		open(PDB, "<$pdb_final")
			or die "ERROR: Could not open PDB file ($pdbc).\n";
	 
		# Loop through original PDB, as before, looking for this chain's ATOM lines
		foreach (<PDB>) {
			if (/^ATOM/) {
			my $cchn = substr $_, 21, 1;
				if ($cchn eq $chain) {
					print IND_FILE $_;
				}
			}
		}
		close PDB;

		print IND_FILE "TER \n";
		close IND_FILE;
	}

	# do not move/remove for first use.

	system ("rm $pdbc*.pdb");

	if (-e "$relax_dir/$pdbc" . "_models/") {
		system ("rm -r $relax_dir/$pdbc" . "_models/");
	}
	if (-e "$relax_dir/$pdbc.clean.pdb") {
		system ("rm -r $relax_dir/$pdbc.clean.pdb");
	}
}



# COUNT REMAINING
open UNI, "remaining.txt";
open REM, ">remaining2.txt";
my $c = 0;
while (<UNI>) {
	chomp;
	if (!-e "$relax_dir/$_.pdb") {print REM "$_\n"; print "REM $_\n"; $c++;}
}
print "-----------------\n$c remaining\n------------------\n";
close UNI;
close REM;

system ("mv remaining2.txt remaining.txt");
system ("ls $relax_dir/*_models|grep .pdb > pdb_models_ready.txt");
system ("cat pdb_models_ready.txt|wc");


exit 0;
