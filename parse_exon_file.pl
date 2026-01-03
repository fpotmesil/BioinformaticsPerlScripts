#!/usr/bin/perl

use warnings;
use strict;

# find_exon.pl
# This script will filter exons with length < min_length
# and count the number of exons selected and filtered out
# G. Lu, University of Nebraska at Omaha

if( (!$ARGV[0]) || (!$ARGV[1]) ) 
{
	print STDERR "
		./find_exo.pl exon_size seqfile;

	Extracts exon sequences with length > s
		exon_size 	: exon size
		seqfile 	: entry file (or stdin)
		\n\n";
}
else
{
	my $exon_size=$ARGV[0];
	my $infile = $ARGV[1];
	my $output = $infile.".out";
	my @exons_input;
	my %exon_hash = ();
	my $exon_key;
#
#
open( EXONS, "<$infile") or die "Cannot open $infile for reading: $!\n\n";
@exons_input = <EXONS>;
close(EXONS);

#
#
open( OUT, ">$output") or die "Cannot open $output for writing: $!\n\n";


	foreach my $line (@exons_input)
	{
		#chomp $line; # get rid of newline at end of string.
		$line =~ s/\r?\n$//;

		if( $line =~ /^>/ )
		{
			#$line =~ s/>//;
			$exon_key = $line;
			$exon_hash{$exon_key} = "";
		}
		else
		{
			$exon_hash{$exon_key} = $exon_hash{$exon_key}.$line;
		}
	}

my $counter = 0;
my @temp_array;

foreach my $temp_key (keys(%exon_hash))
{
	#print "Key: $temp_key\n";
	#print "Value: $exon_hash{$temp_key}\n\n";

	
	if( $temp_key =~ m/(.)(length\()(\d+)(\).)/ )
	{
		#print "MATCH FOR Length is: $3\n\n";
		
		if( $3 > $exon_size )
		{
			print OUT "$temp_key\n";
			print OUT "$exon_hash{$temp_key}\n\n";
		}
	}
	else
	{
		@temp_array = split(/\|/,$temp_key);

		if( ($temp_array[4]-$temp_array[3]) > $exon_size )
		{
			print OUT "$temp_key\n";
			print OUT "$exon_hash{$temp_key}\n\n";
		}

		#print "Key lenghts are: $temp_array[3] and $temp_array[4]\n";
		#print "Exon length: " . ($temp_array[4]-$temp_array[3]) . "\n\n";
		#print "Value: " . "$exon_hash{$temp_key}\n\n";
	}


	#$counter++;
	#if( $counter >= 10 )
	#{
	#	last;
	#}
}


	close(OUT);
}


#
# EOF
#
