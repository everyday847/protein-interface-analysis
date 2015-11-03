#!/usr/bin/perl

use strict;
use warnings;

sub trim {
	my $x = $_[0];
	$x =~ s/^\s*//;
	$x =~ s/\s*$//;
	return $x;
}

# your list of PDB codes of interest
open UNIQUE, "unique.txt";
my %uni;
while (<UNIQUE>) {
	chomp;
	$uni{$_} = 1;
}

system ("ls pdb_files_new/*_models|grep .pdb > pdb_models_ready.txt");
open LIST, "pdb_models_ready.txt";

my @files = <LIST>;

# Feel free to implement a batching system like in make_alascan_pbs.pl.
foreach my $file (@files) {
	chomp $file;
	my $pdbcode = substr $file, 0, 4;
	if (!defined $uni{$pdbcode} ) {next;} 
	my $firstpart = $file;
	$firstpart =~ s/\.pdb//g;

	# Skip empty files (e.g. the extra file make_models produces)
	if (-z "pdb_files_new/$pdbcode" . "_models/$file"){next;}
	
	# If no script has been made for this one before
	# and it hasn't been relaxed before
	# and ALL its models haven't been relaxed before
	if (!-e "relax_pbs/relax_$firstpart.pbs" and !-e "$firstpart" . "_0001.pdb" and !-e "pdb_files_new/$pdbcode.pdb") { #this model hasn't been relaxed yet
		print "$file\n";
		open OUT, ">relax_pbs/relax_$file.pbs";
		my $time = "8";
		my $mem = "4";
		print OUT "#PBS -l nodes=1:ppn=1,mem=$mem" . "GB,walltime=$time:00:00
#PBS -N relax_$firstpart
#PBS -M amw579\@nyu.edu
#PBS -m abe
#PBS -e /scratch/amw579/new_HIPP/\${PBS_JOBNAME}.e\${PBS_JOBID}
#PBS -o /scratch/amw579/new_HIPP/\${PBS_JOBNAME}.o\${PBS_JOBID}
		
cd /scratch/amw579/new_HIPP/\n\n"; 
		
		# These flags haven't been optimized--it's possible that -ex3 and -ex4 are 
		# critical or that Cartesian minimization is a good idea. Experiment!
		print OUT "if [[ ! -f pdb_files_new/$pdbcode.pdb ]] ; then\n\t/work/amw579/Rosetta/main/source/bin/relax.linuxiccrelease -database /work/amw579/Rosetta/main/database/ -relax:coord_constrain_sidechains -relax:ramp_constraints false -ex1 -ex2 -use_input_sc -flip_HNQ -no_optH false -relax:constrain_relax_to_start_coords -in:file:fullatom -ignore_zero_occupancy false -talaris2014 -linmem_ig 10 -in:file:s pdb_files_new/$pdbcode" . "_models/$file\nfi\n";
		close OUT;
	}
}
