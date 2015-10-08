#!/usr/bin/perl

use strict;
use warnings;

# If you have a lot of PDBs and you want to get the easy ones done first,
# one quick way is to get all the PDBs with fewer than n residues and do those.
#
# This can actually be very helpful if, for example, you want to start by
# qsubbing a bunch of jobs with 8 hour walltime, each of which has 20 PDBs to
# process. If most of those PDBs contain 200 residues, those jobs will probably
# finish in 3-4 hours total at most. If one of those PDBs contains 1000 
# residues, that one will itself take over eight hours!

my ( $maximum_residues, $relax_dir ) = @ARGV;

if ( not defined $maximum_residues or not defined $relax_dir ) {
	print "\n";
	print "Run this script by: perl curl.pl [max_res] [relax_dir] \n\n";
	print "\tmax_res: the most residues you want a single PDB model to contain to be relaxed\n";
	print "\trelax_dir: where your model files and relaxed PDBs are.\n\n\n";
	exit;
}

system ("ls $relax_dir/*_models|grep .pdb > pdb_models_ready.txt");
open MODELS, "pdb_models_ready.txt";
open OUT, ">pdb_models_u$maximum_residues.txt";
while (<MODELS>) {
	my $line = $_;
	chomp $line;
	my $pdbc = substr $line, 0, 4;
	my $fn = "$relax_dir/$pdbc\_models/$line";
	if (-z $fn) {next;}
	if (!-e $fn) {next;}
	
	# Sometimes the size of a model is much less than the other models.
	# It might be missing one or more chains. This is unusual but can happen.
	# It isn't your fault--it's a PDB thing. Omit these.
	if (-s $fn < ((-s "$relax_dir/$pdbc\_models/$pdbc\_m1.pdb") - 20000 )) {system ("rm $fn"); next;}

	open PDB, $fn;
	my $res_number = "XXXX";
	my $n_residues = 0;
	while (<PDB>) {
		if (/^ATOM/) {
			if (substr $_, 22, 4 ne $res_number) {
				$n_residues++;
			}
		}
	}
	if ($n_residues < $maximum_residues) {
		print "$line\n";
		print OUT "$line\n";
	}
}

system( "rm pdb_models_ready.txt" );
