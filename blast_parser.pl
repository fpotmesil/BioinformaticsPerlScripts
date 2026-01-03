#!/usr/bin/perl

use strict;
use warnings;


use File::Basename;
use File::Path;
use Getopt::Long;

#
# to help me see what is going on, 
# this produces a lot of output!!
#
my $debug = 0;

#
# the file to open, must be genbank format.
#
my $infile = "";

#
# this is the output file name.
#
my $outfile = "";

#
# the file with the exon sequences 
#
my $exon_file = "";


#
# this is the minimum bit score we need to 
# consider the gene to have a copy.
#
my $min_bit_score = 100;

GetOptions( "outfile=s"=>\$outfile,
        "infile=s"=>\$infile,
        "exon_file=s"=>\$exon_file,
        "min_bit_score=i"=>\$min_bit_score,
        "debug"=>\$debug );

if( (!$infile) || (!$outfile) || (!$exon_file) )
{
    usage();
}



if( $debug )
{
    print "Running in debug mode!\n";
}

print "The Minimum Bit Score we are considering is: $min_bit_score\n";


#
# this gets the path from the input file...
# we add this path to the output filenames for now.
#
my $path;

if( $infile )
{
    (undef,$path,undef) = fileparse($infile);

    my $temp_outfile_name = $outfile;

    $outfile = ">". $path . substr($temp_outfile_name,0,40) . ".fas";
    open(OUTPUT,$outfile) || die "Can't open the output file $outfile: $! \n";
}

#
# this is the hash that will store the 
# exon name as the key and the nucleotide string as the value.
# we can then use this to get the Single Copy Exons.
#
my %exon_hash = ();
my $exon_key;


#
# this is the array that will store the entries until,
# we are ready to write everything to a file.
# we use this to only store unique entries,
# by means of a very inefficient search through the
# array before pushing a new entry in.
#
my @temp_output;


#
# the tabbed blast output headers from "blastall -m 9" 
# look just like so...
# gives us all the field id's and such...
#
# BLASTN 2.2.11 [Jun-05-2005]
# Query: lcl|ENSDARG00000024951:ENSDARE00000085456:203387:204325:1
# Database: zebra_exons.fas
# Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score
#
#

my $query_id;
my $subject_id;
my $identity;
my $align_length;
my $mismatches;
my $gap_openings;
my $query_start;
my $query_end;
my $seq_start;
my $seq_end;
my $e_value;
my $bit_score;

my $record;
my @record_lines;

my $line;
my @line_fields;

my @query_fields;
my $query_name;


my $old_record_separator = $/;

#
# make it so we read in each sequence search at a time,
# we call that a 'paragraph' and parse accordingly...
#
local $/ = "# BLASTN 2.2.11 [Jun-05-2005]\n";	

# BLASTN 2.2.11 [Jun-05-2005]
open( INFILE, "<$infile") or die "Cannot open file: $!\n\n";

#
# count the number of exons that we process.
#
my $loop_counter = 0;

#
# count of the number of exons that we find as Single Copy,
# and hopefully write out to the file along with the exon's
# corresponding sequence...
#
my $number_single_copy = 0;


my $found_high_score = 0;
#
# I need to do this to get rid of the first 
# line of the input file, it is considered a 
# 'record', but is only the record separator 
# line, "# BLASTN 2.2.11 [Jun-05-2005]"...
#
my $garbage = <INFILE>;

#if( <INFILE> )
#{
    while( $record = <INFILE> )
    {
        $found_high_score = 0;
        @record_lines = split(/\n/,$record);

        foreach $line (@record_lines)
        {
            next if( $line =~/BLASTN|Query|Database|Fields/ );
            @line_fields = split(/\t/,$line);

            $query_id = $line_fields[0];
            $subject_id = $line_fields[1];
            $identity = $line_fields[2];
            $align_length = $line_fields[3];
            $mismatches = $line_fields[4];
            $gap_openings = $line_fields[5];
            $query_start = $line_fields[6];
            $query_end = $line_fields[7];
            $seq_start = $line_fields[8];
            $seq_end = $line_fields[9];
            $e_value = $line_fields[10];
            $bit_score = $line_fields[11];

            @query_fields = split(/\|/,$query_id);
#            print "Query Exon Name (without lcl): $query_fields[1]\n";

#
# if the query and subject exon id's do NOT match, then we 
# check the bit score...
# if the bit score is above the threshold, set the flag, 
# this is probably a copy exon... and just continue...
#
#sci_number($e_value);
#if( $bit_score>$min_bit_score  ) # || (0.5>sci_number($e_value)) )

            if( $query_fields[1] ne $subject_id )
            {
                if( $bit_score>$min_bit_score  )
                {
#                    if( $bit_score < $min_bit_score )
#                    {
#            print "Query ID: $query_id\n";
#            print "Subject ID: $subject_id\n";
#            print "Bit Score: $bit_score\n";
#            print "E-Value: $e_value\n";
#                    }
            
                    $found_high_score = 1;
                    last;
                }
            }

#my @query_fields = split(/:/,$query_id);
#my $query_strand_length = $query_fields[3]-$query_fields[2];
#my $coverage = ($align_length/$query_strand_length)*100;
#print "Query Sequence Length: $query_strand_length\n";
#print "Percent Coverage: $coverage\n";
#            print "Query ID: $query_id\n";
#            print "Subject ID: $subject_id\n";
#            print "Identity: $identity\n";
#            print "Alignment Length: $align_length\n";
#            print "Mismatches: $mismatches\n";
#            print "Gap Openings: $gap_openings\n";
#            print "Query Start: $query_start\n";
#            print "Query End: $query_end\n";
#            print "Sequence Start: $seq_start\n";
#            print "Sequence End: $seq_end\n";
#            print "E-Value: $e_value\n";
#            print "Bit Score: $bit_score\n";
#print "\n\n\n";
        }


#
# now, if we have gone through all of the lines in a particular record,
# we have the flag to look at and determine if we found other exons that
# had a bit score value over the minimum bit score.
# if there was nothing over the minimum score,
# then we consider the query a single copy...
#
        if( ! $found_high_score )
        {

            store_unique_entry( ">$query_id" );
#print OUTPUT "$query_id\n";
        }

        $loop_counter++;
#        last if($loop_counter>=10);
    }
#}



#
# Reset the input record separator....
#
local $/ = $old_record_separator;


#
# this is the exon sequences file.
#
open( EXON_FILE, "<$exon_file") or die "Cannot open file: $!\n\n";
my @exon_sequences = <EXON_FILE>;
close(EXON_FILE);

#
# now, we want to go through the files of exons, 
# and find what we have in our temporary output.
# then we will have the exons and the corresponding
# nucleotide strings that we are looking for.
#
foreach my $line (@exon_sequences)
{
    chomp $line; # get rid of newline at end of string.

    #$line =~ s/\r?\n$//;

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

#
# now we have a hash populated with the exon sequences, 
# with the exon header as the keys...
#

#foreach my $temp_key (keys(%exon_hash))
#{
#print "Key: $temp_key, Value: $exon_hash{$temp_key}\n\n";
#print OUTPUT "$temp_key\n$exon_hash{$temp_key}\n";
#}


#
# now, we go through the original exon input file so we can extract the 
# sequences that go with the headers...
#
foreach( @temp_output )
{
    #print ("ENTRY: " . $_ . "\n");
    #print OUTPUT $_;
    print OUTPUT ("$_\n$exon_hash{$_}\n");
    #print ("Exon: $_\n$exon_hash{$_}\n\n");

    $number_single_copy++;
}


print "\nWe parsed $loop_counter records!!\n";
print "Wrote $number_single_copy putative Single Copy exons to file.\n";

close (INFILE);
close (OUTPUT);



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





sub usage
{
    print "For now, look at the source code please!!\n";
    exit(-1);
}


sub sci_number
{
    my $sci_alpha = shift;
    my $rval = 0;

    if( $sci_alpha =~ /(\d+)(e-)(\d+)/ )
    {
#print "Found Match: $sci_alpha\n";
        
        $rval = "0.";

        for( 2..$3 )
        {
            $rval .= "0";
        }
        
        $rval .= "$1";
#        print "Rval: $rval\n";
    }
    else
    {
#        print "NOT A Match: $sci_alpha\n";
        $rval = $sci_alpha;
    }

    return $rval;
}


#
# End Of File.
#
