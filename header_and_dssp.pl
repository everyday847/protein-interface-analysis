#!/usr/bin/perl

# Downloads the special header containing biological and physical information
# (used in proc_ala_scan; it's considerably more convenient if you do this
# step separately so that PAS doesn't have to slow down its work to download files)
#
# Also uses dssp to perform DSSP secondary structure prediction on said files.

$DSSP = "/path/to/dssp/executable";

my ( $relax_dir, $pdb_header_dir, $dssp_dir ) = @ARGV;

if ( not defined $relax_dir or not defined $pdb_header_dir or not defined $dssp_dir ) {
	print "\n";
	print "Run this script by: perl header_and_dssp.pl [relax_directory] [header_dir] [dssp_dir]\n\n";
	print "\trelax_dir: where the refined PDB files are stored\n";
	print "\theader_dir: where the headers will be stored\n";
	print "\tdssp_dir: where the dssp files will be stored\n\n\n";
	exit;
}

# Detect all the files that have been successfully relaxed and download their headers.
system ("ls $relax_dir|grep -G \"^....\.pdb\" > to_be_processed.txt");

open LIST, "to_be_processed.txt";

while (<LIST>) {
    my $pdb_code = $_;
    chomp $pdb_code;
    $pdb_code =~ s/\.pdb//;
    
	my $pdbc = substr $pdb_code, 0, 4;
    # use this to fetch organism etc information
    while (!-e "$pdb_header_dir/$pdbc.pdb" || -z "$pdb_header_dir/$pdbc.pdb" ){
    	print "$pdbc\n";
        system ("curl -o $pdb_header_dir/$pdbc.pdb http://www.rcsb.org/pdb/files/$pdbc.pdb?headerOnly=YES");
    }

    my $dssp = lc $pdbc;
    if ( !-e "$dssp_dir/$dssp.dssp") {
		if (-e "$relax_dir/$pdbc.pdb") {
			system ("$DSSP $relax_dir/$pdbc.pdb $dssp_dir/$dssp.dssp");
		}
    }
}

system( "rm to_be_processed.txt" );
