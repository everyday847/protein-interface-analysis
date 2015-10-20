#!/usr/bin/perl
use bignum;
use strict;
use warnings;

my $dssp_executable = "/path/to/dssp";

sub trim {
	my $x = $_[0];
	$x =~ s/^\s+//g;
	$x =~ s/\s+$//g;
	return $x;
}

my ( $output_dir, $prep_dir, $dssp_dir, $relax_dir ) = @ARGV;

if ( not defined $output_dir or not defined $prep_dir or not defined $dssp_dir or not defined $relax_dir) {
    print "\n";
    print "Run this script by: \n";
	print "perl prep_ala_out.pl [output_directory] [prepped_dir] [dssp_dir] [relax_dir]\n\n";
	print "This will parse the alanine scanning output files in output_directory\n";
	print "and put the results in prepped_dir.\n\n";
	print "dssp_dir contains the dssp files that'll get incorporated into the \n";
	print "final prepped output.\n\n";
	print "relax_dir contains the relaxed PDB files that you might need to make dssp output from\n";
	print "final prepped output.\n\n\n";
    exit;
}

system ("ls $output_dir >lsrout.txt");
#system ("perl dl_headers.pl");
open ROUTS, "lsrout.txt";

CODELOOP: while (<ROUTS>) {
	my $file = $_;
	chomp $file;
	$file =~ s/\.txt//g;
	print "Prepping $file\n";
		
	my $pdbc = substr $file, 0, 4;
	my $dssp = lc($pdbc);
	unless (open DSSP, "$dssp_dir/$dssp.dssp") {
            if (-e "$relax_dir/$pdbc.pdb") {
                system ("$dssp_executable $relax_dir/$pdbc.pdb $dssp_dir/$dssp.dssp");
	        open DSSP, "$dssp_dir/$dssp.dssp" or next;
            } else { system ("rm $output_dir/$file.txt"); next; } 
        }
	my %dssp;
	while (<DSSP>) {
		if (/^  #  RESID/) {last;}
	}
	while (<DSSP>) {
		my $aaaa = $_;
		my $rn = trim (substr $aaaa, 5, 6);
		my $chn = substr $aaaa, 11, 1;
		my $dssptrial = trim (substr $aaaa, 14, 3);
		if ($dssptrial eq "") {$dssptrial = "L";}
		if ($rn ne "") { $dssp{"$chn$rn"} = $dssptrial;} #int "dssp for $rn appears to be $dssptrial\n";}	
	}
	
	close DSSP;

	open ROUT, "$output_dir/$file.txt";
	open PREP, ">$prep_dir/$file" . "_prepped.txt";
	my $trigger = 0;
	my $fail = 1;


	while (<ROUT>) {
		my $line = $_;
		if ($line =~ /protocols.geometry.RB_geometry: centroids_by_jump_int called but no interface detected!!/
		 or $line =~ /protocols.rosetta_scripts.ParsedProtocol: Mover Docking reports failure!/) {
			close ROUT;
			# At this point, if you chose you could also archive the output in a directory of your choice instead of removing it entirely. 
			#system ("mv $output_dir/$file.txt $some_other_directory/");
			system ("rm $output_dir/$file.txt");
			print PREP " ";
			close PREP;
			next CODELOOP;
		}
		if ($line =~ /^protocols.jd2.JobDistributor: ...._._._0001 reported success/) {
			$fail = 0;
		}
		if ($trigger == 1 and $line =~ /^protocols.rosetta_scripts.ParsedProtocol.REPORT: /) { 
			print $line; 
			$line =~ s/protocols.rosetta_scripts.ParsedProtocol.REPORT: //g;
			my @lineobjs = split(" ", $line);
			my $rn = $lineobjs[1];
			my $chn = $lineobjs[2];
			if (!defined $dssp{"$chn$rn"}) {$dssp{"$chn$rn"} = "L";}
			if ($rn > 2 ** 60) {print "caught!";$rn = -2 ** 64 + $rn; } 
			print PREP $lineobjs[0] . "\t$rn\t$chn\t" . $dssp{"$chn$rn"} . "\t" . $lineobjs[4] . "\n";
			last; 
		} else {$trigger = 0;}
		if ($line =~ /^protocols.simple_filters.AlaScan: Energy /) { $trigger = 1; }
	}
	
	while (<ROUT>) {
		my $line = $_;
		if ($line =~ /^$/ or $line =~/^protocols.rosetta_scripts.ParsedProtocol.REPORT: ============End report for scan==================/) {last;}
		if ($line =~ /^protocols.jd2.JobDistributor: ...._._._0001 reported success/) {
			$fail = 0;
		}
		my @lineobjs = split (" ", $line);
		my $rn = $lineobjs[1];
		my $chn = $lineobjs[2];
		if (!defined $dssp{"$chn$rn"}) {$dssp{"$chn$rn"} = "L";}
		if ($rn > 2 ** 60) {print "caught!";$rn = -2 ** 64 + $rn; } 
		print PREP $lineobjs[0] . "\t$rn\t$chn\t" . $dssp{"$chn$rn"} . "\t" . $lineobjs[4] . "\n";
	}

	while (<ROUT>) {
		if (/^protocols.jd2.JobDistributor: ...._._._0001 reported success/) {
			$fail = 0;
		}
	}
	close ROUT;
	close PREP;
	if ($fail == 1) {
		print "$file failed \n";
		system ("rm $prep_dir/$file" . "_prepped.txt");
		system ("rm $output_dir/$file" . ".txt");
	} else {
		# At this point, if you chose you could also archive the output in a directory of your choice instead of removing it entirely. 
		#system ("mv $output_dir/$file.txt $some_other_directory/");
		system ("rm $output_dir/$file.txt");
	}
}

system( "rm lsrout.txt" );
