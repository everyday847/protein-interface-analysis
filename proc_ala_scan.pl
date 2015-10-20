#!/usr/bin/perl
# This script processes Rosetta alanine scan output (our custom version)
# and triggers any necessary Naccess jobs.
# A. Watkins

use strict;
use warnings;

# Naccess path
# Note! You probably already have this in Naccess_all.pl!
my $naccess = "/path/to/Naccess/naccess";


# There are so many potential command line arguments (that are basically facts of your directory structure)
# that we'll just take care of all of them here. This way, you don't have to input a bunch of paths every time.
# If you want to use them as command line arguments, just use the same kind of syntax as you already saw in the
# other scripts.

# where the parsed alanine scanning output lives
my $prepped_dir = "prepped_alaout";
# where the two-chain PDB files live
my $partner_dir = "pdb_partner_files";
# where the files that have been processed by this script should be moved (so we don't repeat)
# another way would be to maintain a list of file names that have been processed, and skip files
# on that list, but I thought this way might be simplest to understand.
my $procced_dir = "procced_alaout";
# the directory where PDB headers are stored
my $header_dir = "pdb_headers";

sub strip {
	my $x = $_[0];
	$x =~ s/^\s*//g;
	$x =~ s/\s*$//g;
	return $x;
}

my $max_len = 0;

use Class::Struct;
struct(residue => [chain => '$', pdbresnum => '$', aa => '$',  ddg => '$', dssp => '$', dsasa => '$']);

my %aa_3_to_1 = ("ala", 'A', "cys", 'C', "asp", 'D', "glu", 'E', "phe", 'F',
				 "gly", 'G', "his", 'H', "ile", 'I', "lys", 'K', "leu", 'L',
				 "met", 'M', "asn", 'N', "pro", 'P', "gln", 'Q', "arg", 'R',
				 "ser", 'S', "thr", 'T', "val", 'V', "trp", 'W', "tyr", 'Y' );

sub yr {
	my $in = $_[0];
	my $year;
	if (length($in) == 2) {
		# do the best we can; assume that 60 means 60s
		if ($in > 60) {
			$year = "19" . $in;
		} else {
			$year = "20" . $in;
		}
	} elsif (length == 4) {
		$year = $in;
	}
	return $year;
}

sub mon {
	my $month;
	# convert to caps
	$_[0] = uc($_[0]);

	if (     $_[0] eq 'JAN') { $month = '01';
	} elsif ($_[0] eq 'FEB') { $month = '02';
	} elsif ($_[0] eq 'MAR') { $month = '03';
	} elsif ($_[0] eq 'APR') { $month = '04';
	} elsif ($_[0] eq 'MAY') { $month = '05';
	} elsif ($_[0] eq 'JUN') { $month = '06';
	} elsif ($_[0] eq 'JUL') { $month = '07';
	} elsif ($_[0] eq 'AUG') { $month = '08';
	} elsif ($_[0] eq 'SEP') { $month = '09';
	} elsif ($_[0] eq 'OCT') { $month = '10';
	} elsif ($_[0] eq 'NOV') { $month = '11';
	} elsif ($_[0] eq 'DEC') { $month = '12';
	} elsif (length($_[0]) == 2) { $month = $_[0]; }

	return $month;
}

# Globals
my %chains;
my %chain1;
my %chain2;
my @this_chains_keys;
my %hotspots;
my $pdb_code;
my $first_chain;
my $second_chain;

system ("ls $prepped_dir > lspreppedalaout");

#open LIST, "$filter.txt";
open LIST, "lspreppedalaout";
print "opened list\n";
open(SQLOUT, ">>old_method.sql");
open (TBLOUT, ">>old_method.tbl");
#print TBLOUT "pdbcode\ttitle\t\tresolution\torganism\tproteins\tkeywords\tinterfacechains\tchain\tsecstruct\tavg_SS_ddg\ttotal_SS_ddg\ttotal_chain_ddg\tpercent_SS_ddg\tavg_SS_dsasa\ttotal_SS_dsasa\ttotal_chain_dsasa\tpercent_SS_dsasa\tnum_SS_hs\ths_SS_ids\ths_SS_relpos\ths_SS_len\tSS_length\tfirst_SS_pos\tlast_SS_pos\tSS_sequence\tarorascore\tnum_faces\n";
print "about to start codeloop\n";
# Loop through each PDB code that needs to be processed today.
CODELOOP: while (<LIST>) {#foreach my $ala_out (@routs) {
	my $ala_out = $_;
	chomp $ala_out;
	#print $ala_out;
#chomp $pdb_file;
   
	# Isolate PDB code from line in interface_chains_pdb.txt
	# Go from "pdb_partner_files/1A0H_A_B.pdb" to "1A0H_A_B"
	$pdb_code = $ala_out;
	$pdb_code =~ s/$prepped_dir\///;
	$pdb_code =~ s/_ala_out_prepped\.txt//;
	my $pdb_file = "$partner_dir/" . $pdb_code . ".pdb";
	
	# Grab chain letters
	if ($pdb_code =~ /.{4}_(.)_(.)/) {
		$first_chain  = $1;
		$second_chain = $2;
	}
   
	##-----------------------------------------------------------------------##
	## Reading in residue data from PDB file --------------------------------##
	##-----------------------------------------------------------------------##
   
	# Hashes for both chains
	undef(%chain1);
	undef(%chain2);
   
	unless (open(PDB, "<$pdb_file")) {
		print "ERROR: Could not open PDB file $pdb_file.\n";
		next;
	}
   
   
	foreach (<PDB>) {
		my $line = $_;
		# Skip if not an ATOM line
		#if (substr $line, 0, 4 ne 'ATOM') {
		if (!($line =~ /ATOM/)) {
			next;
		}
		#print $line;
   
		my $pdb_res   = substr $line, 22, 5;
		$pdb_res =~ s/^\s+//;
		$pdb_res =~ s/\s+$//;
			   
		my $pdb_chain = substr $line, 21, 1;
			   
		my $pdb_aa = substr $line, 17, 3;
		if ($pdb_chain eq $first_chain) {
			$chain1{$pdb_res} = new residue;
			$chain1{$pdb_res}->chain	 ($pdb_chain);
			$chain1{$pdb_res}->pdbresnum ($pdb_res);
			$chain1{$pdb_res}->aa		($aa_3_to_1{lc $pdb_aa});
		} elsif ($pdb_chain eq $second_chain) {
		   
			$chain2{$pdb_res} = new residue;
			$chain2{$pdb_res}->chain	 ($pdb_chain);
			$chain2{$pdb_res}->pdbresnum ($pdb_res);
			$chain2{$pdb_res}->aa		($aa_3_to_1{lc $pdb_aa});
		   
		}
	   
	}
   

	##-----------------------------------------------------------------------##
	## Processing the 2-chain Naccess file ----------------------------------##
	##-----------------------------------------------------------------------##
   
	# produces "3SL9_A_C.rsa"
	unless (-e "$partner_dir/$pdb_code.rsa") {
		print "$pdb_code: Running 2-chain naccess \n";
		system($naccess . " $partner_dir/$pdb_code.pdb");
		system("mv $pdb_code.rsa $partner_dir");
		system("rm $pdb_code.log");
		system("rm $pdb_code.asa");
	}

	open(COMPLEXIN, "$partner_dir/$pdb_code.rsa")
		or die ("Couldn't open the complex file from naccess");
   
	my %complexsasas;
   
	foreach (<COMPLEXIN>) {
		if (/^RES/) {
			my @line = split;
			# %complexsasas is a hash of hashes
			# first index: chain (e.g. "A")
			# second index: resid (e.g. 42)
			my $chn = $line[2];
			my $resid = $line[3];
			my $protocomplexsasas = $line[4];
			$complexsasas{$chn}{$resid} = $protocomplexsasas;
		}
	}
   
	close COMPLEXIN;
   
	##-----------------------------------------------------------------------##
	## Reading in the AlaScan output file -----------------------------------##
	##-----------------------------------------------------------------------##
   
	# Alanine scanning output in files like rosetta_output/1A0H_A_B_ala_out.txt
	unless (open(FILE, "<$prepped_dir/$ala_out")) {
		print STDERR "Can't open ala out file [$ala_out]. Skipping for now\n";
		next;
	}
	print "Considering $pdb_code...\n";
	while (<FILE>) {
		# If there's no interface then we don't need to do anything more.
		if (/^ /) {close FILE; system ("mv $prepped_dir/$ala_out $procced_dir/$ala_out"); next CODELOOP;}
		my @lineobjs = split;
		my $i = $lineobjs[1];
		if ($lineobjs[2] eq $first_chain) {
			if (defined $chain1{$i}) {
				$chain1{$i}->ddg	   ($lineobjs[4]);
				$chain1{$i}->dssp	  ($lineobjs[3]);
			} else {
				next;
			}
		} else {
			if (defined $chain2{$i}) {
				$chain2{$i}->ddg	   ($lineobjs[4]);
				$chain2{$i}->dssp	  ($lineobjs[3]);
			} else {
				next;
			}
		}
	}
   
	##-----------------------------------------------------------------------##
	## Initialize per-chain structures --------------------------------------##
	##-----------------------------------------------------------------------##
   
	undef (%chains);
 
	$chains{$first_chain}  = \%chain1; # $chains{"A"} is reference to %chain1	
	my $chain1_len = scalar (keys %chain1);
   
	$chains{$second_chain} = \%chain2;
	my $chain2_len = scalar (keys %chain2);
	
	undef %hotspots;
	foreach my $ch (keys %chains) {
		foreach my $id (keys %{$chains{$ch}}) {
			$hotspots{$ch}{$id} = 0;
		}
	}
   
	##-----------------------------------------------------------------------##
	## Per-chain processing -------------------------------------------------##
	##-----------------------------------------------------------------------##
   
	foreach my $chn (sort { $a cmp $b } keys %chains) {
		#print $chn . "\n";
		my $pdb4lett = substr($pdb_code, 0, 4);
	   
		# Perform Naccess's dSASA calculations on this file unless done already
		unless (-e "$partner_dir/$pdb4lett\_$chn.rsa") {
			if (-e "$partner_dir/$pdb4lett\_$chn.pdb") {
			print "$pdb_code: Running naccess on the single-chain file.\n";
			system ($naccess . " $partner_dir/$pdb4lett\_$chn.pdb");
			system ("mv $pdb4lett\_$chn.rsa $partner_dir");
			system ("rm $pdb4lett\_$chn.log ");
			system ("rm $pdb4lett\_$chn.rsa ");
			} else {
				print STDERR "had to move on :-( \n";
				next;
			}
		}		
	   
		# Process the result.
		open(CHAININ, "<$partner_dir/$pdb4lett\_$chn.rsa") or die;
		foreach (<CHAININ>) {
			if (/^RES/) {
				# this is a residue line
				my @line = split;
			   
				my $resid = $line[3];
				my $sasa = $line[4];
			   
				if (!defined $chains{$chn}{$resid}) { next; }
				if (!defined $complexsasas{$chn}{$resid} ) {next; } 
				$chains{$chn}{$resid}->dsasa ($complexsasas{$chn}{$resid} - $sasa);
			}
		}
		close CHAININ;
		
	##-------------------------------------------------------------------##
		## Secondary structure and hotspot analysis--------------------------##
		##-------------------------------------------------------------------##
	   
		my $num_H_res = 0; # num of residues in a one-turn-plus helix
		my $num_E_res = 0; # num of residues in a strand
		my $secstruct = ""; # string to store the secondary structures present
		my $last_SS_pos = 0;
   
		undef(@this_chains_keys);
		foreach (keys %{$chains{$chn}}) {
			if ($_ !~ /[A-z]/) {
				push @this_chains_keys, $_;
			}
		}
		
		@this_chains_keys = sort { $a <=> $b } @this_chains_keys;
		
	my $current_SS = "none	  ";
		for (my $n_chain_index = 0; $n_chain_index < scalar @this_chains_keys; $n_chain_index++) {
			my $chain_i = $this_chains_keys[$n_chain_index];
			#print "at $chain_i\n";
			# Record this as a hotspot if it meets the criterion
			if (defined ($chains{$chn}{$chain_i}->ddg) &&
						 $chains{$chn}{$chain_i}->ddg ne 'n/a' &&
						 #$chains{$chn}{$chain_i}->ddg <= 8 &&
						 $chains{$chn}{$chain_i}->aa ne 'G' &&
						 $chains{$chn}{$chain_i}->aa ne 'P' &&
						 $chains{$chn}{$chain_i}->aa ne 'A' &&
						 $chains{$chn}{$chain_i}->ddg >= 1) {
		   
				$hotspots{$chn}{$chain_i} = 1;
			}
			   
			if (!exists ($chains{$chn}{$chain_i}) || !defined ($chains{$chn}{$chain_i}->dssp)) {
				#print "blank at $chain_i\n";
				if ($num_H_res != 0) {
						$current_SS  = "none	 ";
						if ($num_H_res >= 4) {
							put_secstruct_to_DB ("alpha helix", $chn, $last_SS_pos-$num_H_res+1, $last_SS_pos);
						}
						$num_H_res = 0;
				}
				if ($num_E_res != 0) {
						$current_SS  = "none	 ";
						if ($num_E_res >= 3) {
							put_secstruct_to_DB ("beta strand", $chn, $last_SS_pos-$num_E_res+1, $last_SS_pos);
						}
						$num_E_res = 0;
				}
			} else {
				if ($chains{$chn}{$chain_i}->dssp eq 'H') {
					if ($current_SS eq "alpha helix") {
						$last_SS_pos = $chain_i;
					} else {
						$num_H_res = 0;
						$current_SS = "alpha helix";
					}
					$num_E_res = 0;
					$num_H_res++;
				} elsif ($chains{$chn}{$chain_i}->dssp eq 'E') {
					if ($current_SS eq "beta strand") {
						$last_SS_pos = $chain_i;
					} else {
						$num_E_res = 0;
						$current_SS = "beta strand";
					}
					$num_H_res = 0;
					$num_E_res++;
				}
				else {
						$current_SS  = "none	 ";
						if ($num_H_res >= 4) {
							put_secstruct_to_DB ("alpha helix", $chn, $last_SS_pos-$num_H_res+1, $last_SS_pos);
						}
						$num_H_res = 0;
						$current_SS  = "none	 ";
						if ($num_E_res >= 3) {
							put_secstruct_to_DB ("beta strand", $chn, $last_SS_pos-$num_E_res+1, $last_SS_pos);
						}
						$num_E_res = 0;
						$current_SS = "none	  ";
				}
			}
		}
	} # end foreach chain
	   
	close(FILE);
  
	system ("mv $prepped_dir/$ala_out $procced_dir/$ala_out");
   
}

#close (LIST);
exit 0;



sub put_secstruct_to_DB {
   
	my $secstruct = $_[0];
	my $chn = $_[1];
	my $first = $_[2];
	my $last = $_[3];
	my $total_chain_ddg = 0;			# total ddg
	my $total_chain_dsasa = 0;		  # total dsasa
	my $total_hs_ddg = 0;			   # total hotspot ddg, formerly $hsddg
	my $total_hs_dsasa = 0;			 # total hotspot dsasa, formerly $hsdsasa
	my $total_SS_ddg = 0;
	my $total_SS_ddg_no_neg = 0;
	my $total_SS_dsasa = 0;
	my $SS_sequence = "";	
	my $total_chain_hs = 0; # number hotspots
	my $total_SS_hs = 0;
	my $hs_1 = 0;
	my $hs_2 = 0;
	my $hs_3 = 0;
   
	my $last_hs_pos = 0;
	my $first_hs_pos = 0;
	my $abs_SS_hs_pos_list = "";
	my $rel_SS_hs_pos_list = "i; ";
   
	my $arorascore;
	my $length = 0;
	my @n_faces = (0,0,0,0,0,0,0);
	   
	for (my $n_chain_index = 1; $n_chain_index < scalar @this_chains_keys; $n_chain_index++) {
   
		my $chain_i = $this_chains_keys[$n_chain_index];
	   
		if (defined ($chains{$chn}{$chain_i}->ddg) && $chains{$chn}{$chain_i}->ddg =~ /^[\d\.\-]+$/ && $chains{$chn}{$chain_i}->ddg <= 8) {
			#print "totally valid ddg" . $chains{$chn}{$chain_i}->ddg . "\n";
			$total_chain_ddg   += $chains{$chn}{$chain_i}->ddg;
		}
		if (defined ($chains{$chn}{$chain_i}->dsasa) && $chains{$chn}{$chain_i}->dsasa =~ /^[\d\.\-]+$/) {
			$total_chain_dsasa += $chains{$chn}{$chain_i}->dsasa;		
		}
		if ($chain_i >= $first && $chain_i <= $last) {
			$length = $length + 1;
			
			$SS_sequence = $SS_sequence . $chains{$chn}{$chain_i}->aa;
			if (defined ($chains{$chn}{$chain_i}->ddg) && $chains{$chn}{$chain_i}->ddg =~ /^[\d\.\-]+$/ && $chains{$chn}{$chain_i}->ddg <= 8) {	
				$total_SS_ddg += $chains{$chn}{$chain_i}->ddg;
				if ($chains{$chn}{$chain_i}->ddg > 0) {
					$total_SS_ddg_no_neg += $chains{$chn}{$chain_i}->ddg;
				}

				if ($secstruct eq "alpha helix" ) {
					for (my $facej = 0; $facej < 7; $facej++) {
						my $factor;
						$factor = ((($chain_i-$first+$facej)*10)%36)/10;
						#print STDOUT ($chain_i - $first+$facej) . "\t$factor\n";
						if ($chains{$chn}{$chain_i}->ddg > 0) {$n_faces[$facej] += $chains{$chn}{$chain_i}->ddg * $factor;}
					} 
				} elsif ($secstruct eq "beta strand") {
					if ($chains{$chn}{$chain_i}->ddg > 0) {
						if ($chain_i%2) {$n_faces[0] += $chains{$chn}{$chain_i}->ddg;$n_faces[1] -= $chains{$chn}{$chain_i}->ddg;}
						else			{$n_faces[0] -= $chains{$chn}{$chain_i}->ddg;$n_faces[1] += $chains{$chn}{$chain_i}->ddg;}
					}
				}
			}
			if (defined ($chains{$chn}{$chain_i}->dsasa) && $chains{$chn}{$chain_i}->dsasa =~ /^[\d\.\-]+$/) {
				$total_SS_dsasa += $chains{$chn}{$chain_i}->dsasa;
			}
	   
			if ($hotspots{$chn}{$chain_i} == 1) {
				$total_chain_hs++;
				if ($chain_i >= $first && $chain_i <= $last) {
					 if ($chains{$chn}{$chain_i}->ddg > $hs_1) {
							$hs_3 = $hs_2;
							$hs_2 = $hs_1;
							$hs_1 = $chains{$chn}{$chain_i}->ddg;
					} elsif ($chains{$chn}{$chain_i}->ddg > $hs_2) {
							$hs_3 = $hs_2;
							$hs_2 = $hs_1 = $chains{$chn}{$chain_i}->ddg;
					} elsif ($chains{$chn}{$chain_i}->ddg > $hs_3) {
							$hs_3 = $chains{$chn}{$chain_i}->ddg;
					}

					$total_hs_ddg   += $chains{$chn}{$chain_i}->ddg;
					if (defined ($chains{$chn}{$chain_i}->dsasa) && $chains{$chn}{$chain_i}->dsasa =~ /^[\d.\-]+$/) {
						$total_hs_dsasa += $chains{$chn}{$chain_i}->dsasa;
					}
					$total_SS_hs++;
					$last_hs_pos = $chain_i;
					if ($first_hs_pos == 0) {
						$first_hs_pos = $chain_i; # save the location of the first hs
					}
					$abs_SS_hs_pos_list .= $chains{$chn}{$chain_i}->aa;
					$abs_SS_hs_pos_list .= $chains{$chn}{$chain_i}->pdbresnum;
					$abs_SS_hs_pos_list .= ", ";
					$abs_SS_hs_pos_list .= $chains{$chn}{$chain_i}->ddg;
					$abs_SS_hs_pos_list .= "; ";
					my $diff = $chain_i - $first_hs_pos;
					if ($diff > 0) {
						$rel_SS_hs_pos_list .= "i + $diff; ";
					}
				}
			}
		}
	}

	##-------------------------------------------------------------------##
	## Add to database --------------------------------------------------##
	##-------------------------------------------------------------------##
	if ($total_SS_ddg <= 0) {return;}   
	my $tablename;
	if ($secstruct eq "alpha helix" ) {$tablename  = "helidb_tbl";}
	else							  {$tablename  = "sheetdb_tbl";}


	# Some calculations for fields
	my $pdbc = substr $pdb_code, 0, 4;

	# use this to fetch organism etc information
	while (!-e "$header_dir/$pdbc.pdb" ){
		system ("curl -o $header_dir/$pdbc.pdb http://www.rcsb.org/pdb/files/$pdbc.pdb?headerOnly=YES"); 
	}
	open ORIG, "$header_dir/$pdbc.pdb" or die ($!);
	my $title;
	my $resolution;
	my $date;
	my $organism;
	my %proteins;
	my $keywords = "";
   
	while (<ORIG>) {
	my $lll = $_;
	#print "Irrelevant line: $lll";
	if (defined $title and defined $resolution and defined $date and defined $organism) {
		#print "last $lll";
		last;
	}
	if ($lll =~ /COMPND\s*\d*\s*MOLECULE:/) {
		chomp $lll;
		$lll =~ s/COMPND\s*\d*\s*MOLECULE: //g;
		$lll =~ s/;*\s*$//g;
		$lll =~ s/\"/\\\"/g;
		$lll =~ s/\'/\\\'/g;
		$lll =~ s/\*/\\\*/g;	
		$lll =~ s/\\\\\"/\\\"/g;
		$lll =~ s/\\\\\'/\\\'/g;
		#print "Molecule line: $lll\n";
		OVERALL: while (<ORIG>) {
			my $nnn = $_;
			chomp $nnn;
			#print "Chain line: $nnn\n";
			if ($nnn =~ /COMPND\s*\d*\s*CHAIN: /) {
				$nnn =~ s/COMPND\s*\d*\s*CHAIN: //g;
				$nnn =~ s/;*\s+$//g;
				#print "Parsed chain line: $nnn\n";
				my @thesechains = split(",", $nnn);
				foreach my $thischain (@thesechains) {
					$thischain = strip($thischain);
					#print "Protein $thischain: $lll\n";
					$proteins{$thischain} = $lll;
				}
				while ( <ORIG> ){
					$nnn = $_;
					chomp $nnn;
					if ($nnn =~ /COMPND\s*\d*\s*SYNONYM/  || $nnn =~ /COMPND\s*\d*\s*MOL/ || $nnn =~ /SOURCE/ ) {
						last OVERALL;
					}
					$nnn =~ s/COMPND\s*\d*\s*//g;
			        $nnn =~ s/;*\s+$//g;
	                #print "Line: $nnn\n";
		            my @thesechains = split(",", $nnn);
		            foreach my $thischain (@thesechains) {
	                    $thischain = strip($thischain);
						#print "Protein $thischain: $lll\n";
						$proteins{$thischain} = $lll;
					}
				}
			}
		}
	}

	if ($lll =~ /ORGANISM_SCIENTIFIC/) {
		#print $lll;
			$organism = $lll;
			$organism =~ s/SOURCE\s+[0-9]+\s+ORGANISM_SCIENTIFIC:*\s*//;
			$organism =~ s/;//;
			$organism =~ s/\s+$//;
			$organism =~ s/\"/\\\"/g;
			$organism =~ s/\'/\\\'/g;
			$organism =~ s/\*/\\\*/g;
		$organism =~ s/\\\\\"/\\\"/g;
		$organism =~ s/\\\\\'/\\\'/g;
			chomp $organism;
		}
		if ($lll =~/ORGANISM_COMMON/) {
			my $organism2 = $lll;
			$organism2 =~ s/SOURCE\s+[0-9]+\s+ORGANISM_COMMON:*\s*//;
			$organism2 =~ s/;//;
			$organism2 =~ s/\s+$//;
			$organism2 =~ s/\"/\\\"/g;
			$organism2 =~ s/\'/\\\'/g;
			$organism2 =~ s/\*/\\\*/g;
		$organism2 =~ s/\\\\\"/\\\"/g;
		$organism2 =~ s/\\\\\'/\\\'/g;
			chomp $organism2;
		$organism = "$organism ($organism2)";
		}
			   
		if ($lll =~ /^TITLE/) {
			if (!defined $title) {
				$title = substr $lll, 10;
			} else {
				$title .= substr $lll, 10;
			}
			$title =~ s/\s*$//;
			$title =~ s/\"/\\\"/g;
		$title =~ s/\'/\\\'/g;
		$title =~ s/\*/\\\*/g;
		$title =~ s/\\\\\"/\\\"/g;
		$title =~ s/\\\\\'/\\\'/g;
		}

	if ($lll =~ /^KEYWDS/) {
		chomp $lll;
		$lll =~ s/KEYWDS\s+\d*\s+//g;
		$lll =~ s/\s*$//g;
			$lll =~ s/\"/\\\"/g;
		$lll =~ s/\'/\\\'/g;
		$lll =~ s/\*/\\\*/g;
		$lll =~ s/\\\\\"/\\\"/g;
		$lll =~ s/\\\\\'/\\\'/g;
		$keywords .= " " . $lll;
	}
				
		if ($lll =~ /^REMARK   2 RESOLUTION./) {
		#print $lll;
			if (!defined $resolution || $resolution ne 'NMR') {
				$resolution = substr $lll, 26;
				$resolution =~ s/\s*$//;
			}
		}
			   
		if ($lll =~ /^EXPDTA/) {
			if ($lll =~ /NMR/) {
				$resolution = "NMR";
			}
		}
			   
		if ($lll =~/^HEADER/) {
			$date = substr $lll, 50, 9;
			$date =~ s/\s*$//;
		}
	}
	   
	my @keywds = split (/,\s+/, $keywords);
	$keywords = join (', ', @keywds);
	$keywords =~ s/^\s+//g;
	$keywords =~ s/\s+$//g;
	#if (length($keywords) > $max_len) {$max_len = length($keywords);}
	# $ probably looks like dd-MMM-yy
	if (!defined $date ){ $date = "01-JAN-00";}
	my @datefields = split '-', $date;
	$date = yr($datefields[2]) . '-' . mon($datefields[1]) . '-' . $datefields[0];
	if (!defined($organism) or $organism eq ' ') {
	#	print "organism unknown";
	#	exit;
			$organism = 'unknown';
	}
	if (!defined $resolution) {$resolution = "n/a";}
	#print $resolution . "\n";
	$resolution =~ s/NGSTROMS\.//;
	close ORIG;

	# The orientation of the two IDs in interfacechains field of the database
	#  changes. The first should be this chain and the second is its partner.
	my $interfacechains;
	if ($chn eq $first_chain) {
		$interfacechains =  substr ($pdb_code, 5, 1) . " " . substr ($pdb_code, 7, 1);
	} else {
		$interfacechains =  substr ($pdb_code, 7, 1) . " " . substr ($pdb_code, 5, 1);
	}
	 
	my $avg_SS_ddg;

	if ($length > 0) {
		$avg_SS_ddg = $total_SS_ddg / $length;
	} else {
		return;
	}
	   
	my $pct_ddg_SS;
	if ($total_chain_ddg != 0) {
		$pct_ddg_SS = $total_SS_ddg / $total_chain_ddg;
	} else {
		$pct_ddg_SS = 'n/a';
	}
	   
	my $avg_SS_dsasa;
		$avg_SS_dsasa = $total_SS_dsasa / $length;
	   
	my $pct_dsasa_SS;
	if ($total_chain_dsasa != 0) {
		$pct_dsasa_SS = $total_SS_dsasa / $total_chain_dsasa;
	} else {
		$pct_dsasa_SS = 'n/a';
	}
   
	if (substr $rel_SS_hs_pos_list, -2 eq "; ") {
		$rel_SS_hs_pos_list = substr $rel_SS_hs_pos_list, 0, -2;
	}
   
	if (substr $abs_SS_hs_pos_list, -2 eq "; ") {
		$abs_SS_hs_pos_list = substr $abs_SS_hs_pos_list, 0, -2;
	}
   
	my $SS_hs_distance = $last_hs_pos - $first_hs_pos + 1;
   
	# more helix ddg is better, scaling linearly
	#if ($avg_SS_ddg ne "n/a") {
	#	$arorascore = $total_SS_hs * $avg_SS_ddg;
	#i}
	$arorascore = $hs_1 + $hs_2 + $hs_3;
	my $onechain = substr $interfacechains, 0, 1;
	my $otherchain = substr $interfacechains, 2, 1;
	if ( !defined $proteins{$onechain}) {
		print "Undefined protein for |$onechain|";
		#exit;
	}
	if ( !defined $proteins{$otherchain} ) {
		print "Undefined protein for |$otherchain|";
		#exit;
	}
	my $protstring = $proteins{$onechain} . " / " . $proteins{$otherchain};	
  
	my $final_n_faces = 100;
	
	if ($secstruct eq "alpha helix") {
	for (my $ffi = 0; $ffi < 7; $ffi++) {
		$n_faces[$ffi] = $n_faces[$ffi]/$total_SS_ddg_no_neg+1;
		if ($n_faces[$ffi]**2 < $final_n_faces**2) {
		$final_n_faces = $n_faces[$ffi];
		}
	}
	} else {
	$final_n_faces = 2-$n_faces[0]/$total_SS_ddg_no_neg;
	if (2-$n_faces[1]/$total_SS_ddg_no_neg < $final_n_faces) {$final_n_faces = 2-$n_faces[1]/$total_SS_ddg_no_neg;}
	} 
 
	# Round DDG averages, percent DDG helix, Arora score, percent helix dsasa, and		
	# average helix dsasa to 4 decimal places		
	$avg_SS_ddg   = sprintf ("%.4f", $avg_SS_ddg)   if ($avg_SS_ddg =~ /\d/);		
	$total_SS_ddg   = sprintf ("%.4f", $total_SS_ddg)   if ($total_SS_ddg =~ /\d/);		
	$total_chain_ddg   = sprintf ("%.4f", $total_chain_ddg)   if ($total_chain_ddg =~ /\d/);		
	$pct_ddg_SS   = sprintf ("%.4f", $pct_ddg_SS)   if ($pct_ddg_SS =~ /\d/);		
	$pct_dsasa_SS = sprintf ("%.4f", $pct_dsasa_SS) if ($pct_dsasa_SS =~ /\d/);		
	$avg_SS_dsasa = sprintf ("%.4f", $avg_SS_dsasa) if ($avg_SS_dsasa =~ /\d/);		
	$arorascore   = sprintf ("%.4f", $arorascore)   if ($arorascore =~ /\d/);
	$final_n_faces   = sprintf ("%.4f", $final_n_faces)   if ($final_n_faces =~ /\d/);
   
   
	my $SQL = "INSERT INTO $tablename (pdbcode, title, date_added, resolution, ";
	$SQL   .=  "organism, protein_names, keywords, ";
	$SQL   .=  "interfacechains, chain, secstruct, ";
	$SQL   .=  "avg_SS_ddg, total_SS_ddg, total_chain_ddg, percent_SS_ddg, ";
	$SQL   .=  "avg_SS_dsasa, total_SS_dsasa, total_chain_dsasa, percent_SS_dsasa, ";
	$SQL   .=  "num_SS_hs, hs_SS_ids, hs_SS_relpos, hs_SS_len, ";
	$SQL   .=  "SS_length, first_SS_pos, last_SS_pos, SS_sequence, arorascore, num_faces) ";
   
	$SQL   .= "VALUES (\"$pdbc\", \"$title\", \"$date\", \"$resolution\", ";
	$SQL   .=  "\"$organism\", \"$protstring\", \"$keywords\", ";
	$SQL   .=  "\"$interfacechains\", \"$onechain\", \"$secstruct\", ";
	$SQL   .=  "\"$avg_SS_ddg\", \"$total_SS_ddg\", \"$total_chain_ddg\", \"$pct_ddg_SS\", ";
	$SQL   .=  "\"$avg_SS_dsasa\", \"$total_SS_dsasa\", \"$total_chain_dsasa\", ";
	$SQL   .=  "\"$pct_dsasa_SS\", \"$total_SS_hs\", \"$abs_SS_hs_pos_list\", ";
	$SQL   .=  "\"$rel_SS_hs_pos_list\", \"$SS_hs_distance\", ";
	$SQL   .=  "\"$length\", \"$first\", \"$last\", \"$SS_sequence\", ";
	$SQL   .=  "\"$arorascore\", \"$final_n_faces\");";
	   
	
	my $TBL =  "$pdbc\t$title\t$date\t$resolution\t$organism\t$protstring\t$keywords\t$interfacechains\t$onechain\t$secstruct\t$avg_SS_ddg\t";
	$TBL   .=  "$total_SS_ddg\t$total_chain_ddg\t$pct_ddg_SS\t";
	$TBL   .=  "$avg_SS_dsasa\t$total_SS_dsasa\t$total_chain_dsasa\t";
	$TBL   .=  "$pct_dsasa_SS\t$total_SS_hs\t$abs_SS_hs_pos_list\t";
	$TBL   .=  "$rel_SS_hs_pos_list\t$SS_hs_distance\t";
	$TBL   .=  "$length\t$first\t$last\t$SS_sequence\t";
	$TBL   .=  "$arorascore\t$final_n_faces";
	   
	# Only put on DB if helical hotspotS		
	if ($total_SS_hs >= 2) {
		if (($secstruct eq "alpha helix" and $length >= 4) or ($secstruct eq "beta strand" and $length >= 3)) {		
		#print SQLOUT "GOOD: ";
		print SQLOUT $SQL . "\n"; # we be debuggini
	print TBLOUT $TBL . "\n";
	#print $max_len . "\n";
		#print "GOOD: ";
		}
	}

	#print $SQL . "\n";
}

close(SQLOUT);

