#!/usr/bin/perl

#
# 
#
#
use strict;
use warnings;


print "starting the first trial of Single Copy Gene Exon finder...\n";

#
# these are the two input files.
#
my @homologs;
my @paralogs;

#
# this is the array that will store the entries until,
# we are ready to write everything to a file.
# we use this to only store unique entries,
# by means of a very inefficient search through the
# array before pushing a new entry in.
#
my @temp_output;

my $homologs_file = "fugu-zebra.blastParse";
my $paralogs_file = "fugu-fugu-zebra-zebra-blastparse";
my $output_file = "compare_output.txt";

#
# these are the exon nucleotide sequence files.
# this somehow needs to be made more dynamic,
# to allow for different input files.
#
my $fugu_exon_file = "fugu";
my $zebra_exon_file = "zebra";

#
# this is the hash that will store the 
# exon name as the key and the nucleotide string as the value.
# we can then use this to get the Single Copy Exons.
#
my %exon_hash = ();
my $exon_key;
my $exon_value;


#
# this is the first file to open...
#
open( PARALOGS, "<$paralogs_file") or die "Cannot open file: $!\n\n";
@paralogs = <PARALOGS>;
close(PARALOGS);

#
# this is the second file that we open...
#
open( HOMOLOGS, "<$homologs_file") or die "Cannot open file: $!\n\n";
@homologs = <HOMOLOGS>;
close(HOMOLOGS);

#
# this is the output file, we will write to stdout also for now.
#
open( OUTFILE, ">$output_file") or die "Cannot open file: $!\n\n";

foreach (@homologs)
{	
	my @homo_strings = split(/\s+/);

#
# if we have an exon that is a homolog, but
# not a paralog, then we have a Single Copy Gene.
# I think that is the algorithm...
#
	
	my $found_paralog = 0;

	foreach (@paralogs)
	{	
		my @para_strings = split(/\s+/);

		if( ($homo_strings[0] eq $para_strings[0])||
			($homo_strings[1] eq $para_strings[1]))
		#if( $homo_strings[0] eq $para_strings[0] )
		{
			$found_paralog = 1;
			last;
		}

	}

	if( ! $found_paralog )
	{
		store_unique_entry( $homo_strings[0] );
	}
}


foreach( @temp_output )
{
	print ("ENTRY: " . $_ . "\n");
	print OUTFILE $_ . "\n";
}


#
# this is the fugu exon file.
# SINFRUE is what the lines will start with.
#
open( FUGU, "<$fugu_exon_file") or die "Cannot open file: $!\n\n";
my @fugu_exons = <FUGU>;
close(FUGU);

#
# this is the zebra exon file.
# ENSDARE is what the lines will start with.
#
open( ZEBRA, "<$zebra_exon_file") or die "Cannot open file: $!\n\n";
my @zebra_exons = <ZEBRA>;
close(ZEBRA);

#
# now, we want to go through the files of exons, 
# and find what we have in our temporary output.
# then we will have the exons and the corresponding
# nucleotide strings that we are looking for.
#
	foreach my $line (@fugu_exons)
	{
		chomp $line; # get rid of newline at end of string.

		if( $line =~ /^>/ )
		{
			$line =~ s/>//;
			$exon_key = $line;
			$exon_hash{$exon_key} = "";
		}
		else
		{
			$exon_hash{$exon_key} = $exon_hash{$exon_key}.$line;
		}
	}

	foreach my $line (@zebra_exons)
	{
		chomp $line; # get rid of newline at end of string.

		if( $line =~ /^>/ )
		{
			$line =~ s/>//;
			$exon_key = $line;
			$exon_hash{$exon_key} = "";
		}
		else
		{
			$exon_hash{$exon_key} = $exon_hash{$exon_key}.$line;
		}
	}


#my $temp_key;

#foreach $temp_key (keys(%exon_hash))
#{
#	print "Key: $temp_key, Value: $exon_hash{$temp_key}\n\n";
#}

#
# now we have all the exons and the values in our hash, 
# we need to go through the output file, get the exon names,
# and do a lookup from the hash to get the nucleotide strings.
#
print OUTFILE "\n\n\nSingle Copy Exons:\n\n";

foreach( @temp_output )
{
	print ("EXON: " . $exon_hash{$_} . "\n");
	print OUTFILE ("Exon: $_\n$exon_hash{$_}\n\n");
}





close(OUTFILE);

print "That is all now! Results in $output_file.  Have a nice day!\n";


#
# a simple search through an array to check for 
# identical entries before storing the entry in 
# the array.  Not the best solution...
#
sub store_unique_entry
{
	my $entry = shift;
	my $found_entry = 0;

	foreach( @temp_output )
	{
		if( $_ eq $entry )
		{
			$found_entry = 1;
			last;
		}
	}

	if( ! $found_entry )
	{
		push( @temp_output, $entry );
	}
}



#
# End Of File.
#
