#!/usr/bin/perl

# Once you have created your models, you will want to look through those models
# and remove any residues that are missing backbone atoms. This is most common 
# for terminal residues (N-term with an amine, in particular) and simply means
# you should truncate.
#
# Rosetta won't discard these residues automatically--instead, it'll quit with
# an error asking you what you want to do. So do this first to pre-empt that
# behavior. Run this script AFTER make_models.pl.

use strict;
use warnings;

my ( $model_dir ) = @ARGV;
 
if ( not defined $model_dir ) {
	print "\n";
	print "This script identifies all the potentially problematic chains\n";
	print "that may be missing backbone atoms. Check these manually, though!\n";
	print "Run this script by: perl remove_incomplete_backbone.pl [model_dir] \n\n";
	exit;
}

sub trim {
	my $x = $_[0];
	$x =~ s/^\s*//;
	$x =~ s/\s*$//;
	return $x;
}

system ("ls $model_dir/*_models|grep .pdb > pdb_models_ready.txt");
open LIST, "pdb_models_ready.txt";

my @files = <LIST>;
my $c = 0;
my $i = 1;

foreach my $file (@files) {
	chomp $file;
	my $pdbcode = substr $file, 0, 4;
	open PDB, "$model_dir/$pdbcode\_models/$file" or next;
	
	my $last_rn = "XXXX";
	my $bb_num = 0;
	while (<PDB>) {
		my $line = $_;
		chomp $line;
		if (/^ATOM/) {
			if ((substr ($line, 22, 4)) ne $last_rn) {
				if ($bb_num < 4 and $last_rn ne "XXXX") {
					print "$file $last_rn\n";
				}
				$last_rn = substr ($line, 22, 4);
				$bb_num = 0;
			}
			if ((substr ($line, 12, 4)) eq " CA " ) {$bb_num++;}
			if ((substr ($line, 12, 4)) eq " N  " ) {$bb_num++;}
			if ((substr ($line, 12, 4)) eq " C  " ) {$bb_num++;}
			if ((substr ($line, 12, 4)) eq " O  " ) {$bb_num++;}
		}
										
	}
	close PDB;
}

system( "rm pdb_models_ready.txt");
