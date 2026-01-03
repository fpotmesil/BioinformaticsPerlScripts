#!/usr/local/bin/perl

use POSIX;
use Bio::SeqIO;
use Getopt::Std;
use Getopt::Long;


#parse_blastout.pl
#G. Lu @ UNO

&init;
&process;

sub usage {
	my $program = `basename $0`;
	chop($program);
	print STDERR "
		$program [ -i ] [ -c ] [-h] seqfile

		Parse the blast result based on %identity and %coverage given. The default values for identity is 80 and 50 for coverage.

		-i  	: %identity
		-c  	: %coverage
		-h 	: this message
		seqfile   : entry file (or stdin)

		"; 
}

sub init {

	getopts(':ich');

	if ($opt_i) {
		$identity=$opt_i;
	} else { $identity = 80;
	}

	if ($opt_c) {
		$coverage = $opt_c;
	} else {
		$coverage = 50;
	}

	if ($opt_h) {
		&usage;
		exit;
	}  

	$seqfile = ($ARGV[0]) ? $ARGV[0] : "-";

}


sub process
{

	open (INPUT_FILE,"$seqfile") || die "cannot open for reading";
	open (OUTPUT_FILE,">$seqfile.blastParseOut") || die "cannot create out file";

	while (<INPUT_FILE>)
	{
		@fields = split;

# find homologe genes, in case of the same genome, paralog genes
		if ($fields[0] ne $fields[1])
		{
# Get the record with Identities >= 50%
			if ($fields[2] >= $identity)
			{
# Get the Query sequence length; Lu modified
				($q_start) = (split (/\|/,$fields[0]))[1];
				($q_end) = (split (/\|/,$fields[0]))[2];
				$X = $q_end - $q_start;
# Get the Subject match positions
				($s_start) = (split (/\|/,$fields[1]))[1];
				($s_end) = (split (/\|/,$fields[1]))[2];
				$Y = $s_end - $s_start; 
				$Coverage = ($fields[3]/MIN($X,$Y))*100;

				if ($Coverage > $coverage)
				{
					print OUTPUT_FILE "$fields[0]\t$fields[1]\t$Coverage\t$fields[2]\n";
				}
			} # End of if

		} # End of if
	} # End of while loop
	close (INPUT_FILE) || die "can't close reading file";
	close (OUTPUT_FILE) || die "can't close out file";

#} # End of for loop


} #End process
#=====================================================================================

sub MIN
{
	if ($_[0] < $_[1])
	{
		return $_[0];
	}
	else
	{
		return $_[1];
	}
}
