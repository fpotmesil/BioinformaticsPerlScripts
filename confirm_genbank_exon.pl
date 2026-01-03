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
# the file to open, must be genbank format.
#
my $infile = "";
my $infiles_dir = "";

my $locate_exon_id = "";

GetOptions( "locate_exon_id=s"=>\$locate_exon_id,
        "infile=s"=>\$infile,
        "infiles_dir=s"=>\$infiles_dir );


if( ( (!$infile)&&(!$infiles_dir) ) || (!$locate_exon_id) )
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
#
# Global File Statistics that I will keep track of...
#
my $number_of_genes = 0;
my $number_of_exons = 0;


#
# variables used each time through the loop.
#
my $seq_in;
my $gene_id = "";
my $exon_id = "";
my $found_gene = 0;

my $found_my_exon = 0;

foreach my $file (@input_files)
{
    print "Opening genbank data file: $file\n";

    $seq_in  = Bio::SeqIO->new( -format=>'genbank', -file=>$file );

    while( my $seq_object = $seq_in->next_seq() )
    {
        try
        {
            $gene_id = "";
            $exon_id = "";
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
                                    if( $locate_exon_id eq $2 )
                                    {
                                        print "Found a Match!!\n";
                                        print "File: '$file', Gene ID: $gene_id\n";

                                        print "Start: ", $feature->start,
                                              " End: ", $feature->end,
                                              " Strand: ", $feature->strand, "\n\n";

                                        $found_my_exon = 1;
                                    }

                                    $number_of_exons++;
                                    $exon_id = $2;

                                }
                            }
                        }
                    }
                }
            }
        }
        catch Bio::Root::Exception with 
        {
            my $error = shift;
            print $error->stringify . "\n";
            print "Skipping output of Gene '$gene_id' Exon '$exon_id'\n\n"
        }
        otherwise
        {
            my $error = shift;
            print "ERROR: 'otherwise' exception!!\n\n";
            print $error->stringify . "\n";
        };
    }

    print "Done Processing file '$file'\n";

    if( $found_my_exon )
    {
        last;
    }
}


print "Overall Run Statistics:\n",
      "Number of Genes: $number_of_genes\n",
      "Total Number of Exons: $number_of_exons\n";

print "That is all, Buh-Bye now, have a splendid day!\n";

sub usage
{
    print "\n\nUSAGE:\n",
          "./confirm_genbank_exons.pl --infile=<infile> --locate_exon_id=<exon_id_to_locate>\n\n";
    print "A genbank format input file must be provided with --infile.\n";
    print "-OR- a complete directory of input files can be provided with ",
          "--infiles_dir=<directory>\n  Currently, the files MUST be *.gz or *.Z files!!\n";
    print "--locate_exon_id is the exon id that we want to search for in the file(s)\n\n";

    exit( -1 );

}




#
# End of File.
#
