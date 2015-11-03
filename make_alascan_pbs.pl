#!/usr/bin/perl

use strict;
use warnings;

# This script is provided as-is. I suggest modifying all the names of 
# the directories yourself, or imitating the syntax that I use in the
# other scripts to process command-line arguments.
#
# I do provide ample comments!

# We need two lists. First, we need a list of all the PDB partner files
# that have been successfully prepped and procced. These are done! We
# also need a list of all the jobs that stopped mid-way through.
# This is an optimization mostly for the situation where you think some
# jobs will take much longer than others. Statistically, the jobs where
# the PBS script exits are probably longer than average. So this way,
# you have the option of treating them separately.
system ("ls prepped_alaout/|grep -G ^...._._._ala_out_prepped\\\.txt|sed s'/_ala_out_prepped.txt//' >pdbpartnerlist_done.txt");
system ("ls procced_alaout/|grep -G ^...._._._ala_out_prepped\\\.txt|sed s'/_ala_out_prepped.txt//' >>pdbpartnerlist_done.txt");
system ("ls rosetta_output/|grep -G ^...._._._ala_out\\\.txt|sed s'/_ala_out.txt//' >too_long.txt");

my %done;
my %too_long;
my %unique;
my %mse;

# Just populate dicts with everything we know. unique_pdbs.txt is a list
# of the PDB 4-letter codes you're interested in.  
open DONE, "pdbpartnerlist_done.txt";
while (<DONE>) {
	chomp;
	$done{$_} = 1;
}

open TOO_LONG, "too_long.txt";
while (<TOO_LONG>) {
	chomp;
	$too_long{$_} = 1;
}

open UNIQUE, "unique_pdbs.txt";
while (<UNIQUE>) {
	chomp;
	$unique{$_} = 1;
} 

# List all partner files.
system ("ls /scratch/amw579/new_HIPP/pdb_partner_files| grep -G '^...._._.\\\.pdb' |sed s'/\.pdb//g' >partners.txt");
open LIST, "partners.txt";
my @files = <LIST>;

# $i is a file number increment. $c is the number of jobs printed to 
# that particular PBS script so far.
my $i = 1;
my $c = 0;
my $cmax = 8;

open OUT, ">alascan_pbs/alascan_$i.pbs";
print OUT "#PBS -l nodes=1:ppn=1,mem=4GB,walltime=48:00:00
#PBS -N ala_scan_$i
#PBS -M amw579\@nyu.edu
#PBS -m abe
#PBS -e /scratch/amw579/new_HIPP/\${PBS_JOBNAME}.e\${PBS_JOBID}
#PBS -o /scratch/amw579/new_HIPP/\${PBS_JOBNAME}.o\${PBS_JOBID}

cd /scratch/amw579/new_HIPP/\n\n";

# rosetta_output is basically a "scratch directory" for the Rosetta logs.
# After the alanine scanning job completes, it is moved to 
# rosetta_output_new (so that incomplete files are never found there and
# the prep_ala_out script can assume its inputs are complete logfiles).
# Also, the alanine scanning script outputs a new file identical to
# the starting file like XXXX_Y_Z_0001.pdb that can be removed.
foreach my $file (@files) {
	chomp $file;
	$file =~ s/\.pdb//g;

	if (!defined $done{$file}) {

		# Separate output for files suspected of being too long!
		# This way, you can mostly have batches of $cmax jobs
		# but for big files, do them alone.
		if (defined $too_long{$file}) {
				open LOUT, ">long_alascan_pbs/$file.pbs";
				print LOUT "#PBS -l nodes=1:ppn=1,mem=4GB,walltime=96:00:00
#PBS -N long_ala_$file
#PBS -M amw579\@nyu.edu
#PBS -m abe
#PBS -e /scratch/amw579/new_HIPP/\${PBS_JOBNAME}.e\${PBS_JOBID}
#PBS -o /scratch/amw579/new_HIPP/\${PBS_JOBNAME}.o\${PBS_JOBID}

cd /scratch/amw579/new_HIPP/\n";
			print LOUT "if [[ ! -f rosetta_output_new/$file\_ala_out.txt && ! -f prepped_alaout/$file\_ala_out_prepped.txt ]]; then\n\t/work/amw579/Rosetta/main/source/bin/rosetta_scripts.linuxiccrelease -database /work/amw579/Rosetta/main/database/ -s /scratch/amw579/new_HIPP/pdb_partner_files/$file.pdb -use_input_sc -jd2:ntrials 2 -ex1 -ex2 -parser:protocol /scratch/amw579/new_HIPP/alascan.xml -talaris2014 > rosetta_output/$file" . "_ala_out.txt\n
\tmv rosetta_output/$file" . "_ala_out.txt rosetta_output_new/ 
\trm $file" . "_0001.pdb \nfi\n";
			#$i++;
				close LOUT;
		} else {

			print OUT "if [[ ! -f rosetta_output_new/$file\_ala_out.txt && ! -f prepped_alaout/$file\_ala_out_prepped.txt ]]; then\n\t/work/amw579/Rosetta/main/source/bin/rosetta_scripts.linuxiccrelease -database /work/amw579/Rosetta/main/database/ -s /scratch/amw579/new_HIPP/pdb_partner_files/$file.pdb -use_input_sc -jd2:ntrials 2 -ex1 -ex2 -parser:protocol /scratch/amw579/new_HIPP/alascan.xml -talaris2014 > rosetta_output/$file" . "_ala_out.txt\n
\tmv rosetta_output/$file" . "_ala_out.txt rosetta_output_new/
\trm $file" . "_0001.pdb \nfi\n";
			$c++;

			if ($c > $cmax) {
				$c = 0;
				$i++;
				close OUT;

				open OUT, ">alascan_pbs/alascan_$i.pbs";
				print OUT "#PBS -l nodes=1:ppn=1,mem=8GB,walltime=48:00:00
#PBS -N ala_scan_$i
#PBS -M amw579\@nyu.edu
#PBS -m abe
#PBS -e /scratch/amw579/new_HIPP/\${PBS_JOBNAME}.e\${PBS_JOBID}
#PBS -o /scratch/amw579/new_HIPP/\${PBS_JOBNAME}.o\${PBS_JOBID}

cd /scratch/amw579/new_HIPP/\n";
			}
		}
	}
}
