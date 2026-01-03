#!/usr/bin/perl

use strict;
use warnings;


use File::Basename;
use File::Path;
use Bio::Seq;
use Bio::SeqIO;
use Getopt::Long;

use Error qw(:try);
$Error::Debug = 1; # enables verbose stack trace 

#
# to help me see what is going on, 
# this produces a lot of output!!
#
my $debug = 0;

#
# the default minimum exon size to write to output,
# this is changed with --exon_size=800 or similiar.
#
my $minimum_exon_size = 500;

#
# the file to open, must be genbank format.
#
my $infile = "";
my $infiles_dir = "";

#
# this is the output file name.
#
my $outfile = "";



GetOptions( "outfile=s"=>\$outfile,
        "infile=s"=>\$infile,
        "infiles_dir=s"=>\$infiles_dir,
        "debug"=>\$debug,
        "exon_size=i"=> \$minimum_exon_size );

print "Minimum Exon Size is $minimum_exon_size\n";

if( $debug )
{
    print "Running in debug mode!\n";
}

if( ((!$infile)&&(!$infiles_dir) ) || (!$outfile) )
{
    usage();
}


#
# this gets the path from the input file...
# we add this path to the output filenames for now.
#
my $path;

#
# this is an array to map any '.gz' filenames to the gunzip command,
# all done automatically and on the fly for input...  cool, huh?
# 
my @input_files;

if( $infile )
{
    (undef,$path,undef) = fileparse($infile);

    @input_files = map { /\.(gz|Z)$/ ? "gzip -dc $_ |" : $_  } $infile;
}
else
{
    (undef,$path,undef) = fileparse($infiles_dir);

    opendir( DIR, $infiles_dir );
#
# if we have a directory for input, I am only accepting
# gzipped (*.gz) or Zipped (*.Z) files for now!!
#
    @input_files = grep( /\.(gz|Z)$/, readdir(DIR) );
    closedir(DIR);

    @input_files = map { /\.(gz|Z)$/ ? "gzip -dc $path$_ |" : $_  } @input_files;
}


foreach (@input_files)
{
    print "File: $_\n";
}

my $temp_outfile_name = $outfile;

$outfile = ">". $path . substr($temp_outfile_name,0,40) . ".fas";
open(OUTPUT,$outfile) || die "Can't open the output file $outfile: $! \n";


if( $debug )
{
    my $logfile = ">". $path . substr($temp_outfile_name,0,40) . ".log";
    open(LOGFILE,$logfile) || die "Can't open the output file $logfile: $! \n";
#
# set the autoflush flag for LOGFILE....
#
    my $saved_FH = select(LOGFILE);
    $| = 1;
    select($saved_FH);
}


#
# Global File Statistics that I will keep track of...
#
my $number_of_genes = 0;
my $number_of_exons = 0;
my $minimum_size_exons = 0;
my $skipped_exons = 0;


#
# variables used each time through the loop.
#
my $seq_in;
my $sequence = "";
my $gene_id = "";
my $exon_id = "";
my $exon_start = 0;
my $exon_end = 0;
my $reverse_strand = 1;
my $found_gene = 0;


foreach my $file (@input_files)
{
    if( $debug )
    {
        print "Opening genbank data file: $file\n";
        print LOGFILE "Opening genbank data file: $file\n";
    }

    $seq_in  = Bio::SeqIO->new( -format=>'genbank', -file=>$file );

    while( my $seq_object = $seq_in->next_seq() )
    {
        try
        {
            $sequence = "";
            $gene_id = "";
            $exon_id = "";
            $exon_start = 0;
            $exon_end = 0;
            $reverse_strand = 1;
            $found_gene = 0;

            foreach my $feature ($seq_object->get_SeqFeatures )
            {
                if( $feature->primary_tag eq 'gene' )
                {
                    my @gene_id = $feature->get_tag_values('gene');
                    $gene_id = $gene_id[0];
                    $found_gene = 1;
                    $number_of_genes++;
                }
                elsif( $feature->primary_tag eq 'exon' )
                {
                    foreach my $tag ($feature->get_all_tags) 
                    {
                        if( $tag eq 'note' )
                        {
                            foreach my $value ( $feature->get_tag_values($tag) )
                            {
                                if( $value =~ /^(exon_id=)(.*)$/ )
                                {
                                    $exon_id = $2;
                                }
                            }
                        }
                    }

                    $reverse_strand = $feature->strand;
                    $exon_start = $feature->start;
                    $exon_end = $feature->end;
                    $number_of_exons++;

                    if( $found_gene && ( ($exon_end-$exon_start) >= $minimum_exon_size) )
                    {
                        $sequence = $seq_object->subseq($feature->start,$feature->end);

                        my $out_sequence = Bio::PrimarySeq->new ( -seq => $sequence,
                                -id=>">lcl|$gene_id:$exon_id:$exon_start:$exon_end:$reverse_strand",
                                -moltype => 'dna' );

                        if( $reverse_strand < 0 )
                        { 
                            $out_sequence = $out_sequence->revcom;
                        }

                        print OUTPUT $out_sequence->id() . "\n";
                        print OUTPUT $out_sequence->seq() . "\n";

                        $minimum_size_exons++;
                    }
                    else
                    {
                        if( $debug )
                        {
                            print LOGFILE "Skipped Exon ",
                                  "$gene_id:$exon_id:$exon_start:$exon_end:$reverse_strand\n";    
                        }

                        $skipped_exons++;
                    }
                }
            }
        }
        catch Bio::Root::Exception with 
        {
            my $error = shift;
            print $error->stringify . "\n";
            print "Skipping output of Gene '$gene_id' Exon '$exon_id', length: ",
                  ($exon_end-$exon_start), "\n\n";

            if( $debug )
            {
                print LOGFILE ($error->stringify . "\n");
                print LOGFILE ("Skipping output of Gene '$gene_id' Exon '$exon_id', length: ",
                        ($exon_end-$exon_start), "\n\n");
            }
        }
        otherwise
        {
            my $error = shift;
            print "ERROR: 'otherwise' exception!!\n\n";
            print $error->stringify . "\n";

            if( $debug )
            {
                print LOGFILE ("ERROR: 'otherwise' exception!!\n\n");
                print LOGFILE ($error->stringify . "\n");
            }
        };
    }

    print "That is all, Buh-Bye now, have a spelndid day!!\n";
}


print "Overall Run Statistics:\n",
      "Number of Genes: $number_of_genes\n",
      "Total Number of Exons: $number_of_exons\n",
      "Exons with length > $minimum_exon_size: $minimum_size_exons\n",
      "Exons with length < $minimum_exon_size: $skipped_exons\n";


close(OUTPUT);

if( $debug )
{
    print LOGFILE "Overall Run Statistics:\n",
          "Number of Genes: $number_of_genes\n",
          "Total Number of Exons: $number_of_exons\n",
          "Exons with length > $minimum_exon_size: $minimum_size_exons\n",
          "Exons with length < $minimum_exon_size: $skipped_exons\n";


    close(LOGFILE);
}



sub usage
{
    print "\n\nUSAGE:\n",
          "./extract_genbank_exons.pl --infile=<infile> --outfile=<outfile>",
          " [--debug] [--exon_size=<minimum_exon_size]\n\n";
    print "A genbank format input file must be provided with --infile.\n";
    print "-OR- a complete directory of input files can be provided with ",
          "--infiles_dir=<directory>\n  Currently, the files MUST be *.gz or *.Z files!!\n";
    print "--outfile=<file> file will be written to in the same directory the infile is in.\n";
    print "--debug flag is optional, opens a logfile of the output file name,\n",
          " with some tracing/debugging information.\n";
    print "--exon_size sets the minimum exon size to write to output file.\n";

    exit( -1 );

}




#
# End of File.
#
