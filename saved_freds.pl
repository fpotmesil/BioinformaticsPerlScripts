#!/usr/bin/perl

use strict;
use warnings;


use File::Basename;
use File::Path;
use Bio::Seq;
use Bio::SeqIO;

#
# to help me see what is going on, 
# this produces a lot of output!!
#
my $debug = 1;

# Ask the user what they want to name the files
print STDOUT "Root name of the output file:\t";
my $input = <STDIN>;
chomp $input;

# Open up the various files we're going to use
(undef,my $path,undef) = fileparse($ARGV[0]);
my $outfile = ">". $path . substr($input,0,24) . ".fas";
open(OUTPUT,$outfile) || die "Can't open the output file $outfile: $! \n";

if ($debug )
{
    my $logfile = ">". $path . substr($input,0,24) . ".log";
    open(LOGFILE,$logfile) || die "Can't open the output file $logfile: $! \n";
}


#local $/ = "//\n";	# Enable paragraph mode. Slurp in each record

@ARGV = map { /\.(gz|Z)$/ ? "gzip -dc $_ |" : $_  } @ARGV;

my $seq_in;

foreach my $file (@ARGV) 
{
    eval
    {
        $seq_in  = Bio::SeqIO->new( -format=>'genbank', -file=>$file );
    };

    if( $@ )
    {
        print "Script FATAL ERROR!!\n$@\n\n";
        exit(-1);
    }

    while( my $seq_object = $seq_in->next_seq() )
    {
#        my $anno_collection = $seq_object->annotation;
#        for my $key ( $anno_collection->get_all_annotation_keys ) 
#        {
#            my @annotations = $anno_collection->get_Annotations($key);
#            for my $value ( @annotations ) 
#            {
#                print "tagname : ", $value->tagname, "\n";
## $value is an Bio::Annotation, and has an "as_text" method
#                print "  annotation value: ", $value->as_text, "\n";
#            }
#        }



        foreach my $feature ($seq_object->get_SeqFeatures )
        {



            for my $tag ($feature->get_all_tags) 
            {
                print "  tag: ", $tag, "\n";
                for my $value ($feature->get_tag_values($tag)) 
                {
                    print "    value: ", $value, "\n";
                }
            }




            if( $feature->primary_tag eq 'exon' )
            {
                if( $feature->strand < 0 )
                {
                    print "complement sequence!!\n";
                }
                else
                {
                    print "NO complement\n";
                }

#                print "Exon start: ",$feature->start," Exon end: ",
#                      $feature->end,".  \nFeature string:   GFF[",
#                      $feature->gff_string,"]\n";




                print "Exon start: ", $feature->start,
                      " Exon end: ", $feature->end,
                      ".  \nExon string: ",
                      $seq_object->subseq($feature->start, $feature->end),"\n";
            }
        }




#        foreach my $feature ( $seq_object->top_SeqFeatures() ) 
#        {
#            if( $feature->primary_tag() eq 'CDS' ) 
#            {
#                foreach my $exon ( $feature->sub_SeqFeature() ) 
#                {
#                    print "Exon start",$exon->start," Exon end ",$exon->end,"\n";
#                }
#            }
#        }


#        foreach my $feature ( $seq_object->top_SeqFeatures ) 
#        {
#            if ( $feature->primary_tag eq 'CDS' ) 
#            {
#                my $cds_object = $feature->spliced_seq;
#                print "CDS sequence is ",$cds_object->seq,"\n";
#            }
#        }


#print $seq_object->seq();

        print STDOUT "\n\nFile Record Break, hit 'Ctrl-D' to continue....\t";
        <STDIN>;
    }


    print "That is all now, have a spelndid day!!\n";



#open(INPUT,$file) || die "Can't open the input file $file: $! \n";
#
#    while( <INPUT> )
#    {
#        print;

## Ask the user what they want to name the files
#print STDOUT "File Record Break, hit 'Ctrl-D' to continue....\t";
#<STDIN>;
#
#    }

}







#
#
#
#
#
#

#!/usr/bin/perl

use strict;
use warnings;


use File::Basename;
use File::Path;
use Bio::Seq;
use Bio::SeqIO;

#
# to help me see what is going on, 
# this produces a lot of output!!
#
my $debug = 1;

# Ask the user what they want to name the files
print STDOUT "Root name of the output file:\t";
my $input = <STDIN>;
chomp $input;

# Open up the various files we're going to use
(undef,my $path,undef) = fileparse($ARGV[0]);
my $outfile = ">". $path . substr($input,0,24) . ".fas";
open(OUTPUT,$outfile) || die "Can't open the output file $outfile: $! \n";

if ($debug )
{
    my $logfile = ">". $path . substr($input,0,24) . ".log";
    open(LOGFILE,$logfile) || die "Can't open the output file $logfile: $! \n";
}


#local $/ = "//\n";	# Enable paragraph mode. Slurp in each record

@ARGV = map { /\.(gz|Z)$/ ? "gzip -dc $_ |" : $_  } @ARGV;

my $seq_in;

foreach my $file (@ARGV) 
{
    eval
    {
        $seq_in  = Bio::SeqIO->new( -format=>'genbank', -file=>$file );
    };

    if( $@ )
    {
        print "Script FATAL ERROR!!\n$@\n\n";
        exit(-1);
    }

    while( my $seq_object = $seq_in->next_seq() )
    {
        my $sequence = "";
        my $gene_id = "";
        my $exon_id = "";
        my $exon_start = 0;
        my $exon_end = 0;
        my $reverse_strand = 1;
        my $found_gene = 0;

        foreach my $feature ($seq_object->get_SeqFeatures )
        {
            if( $feature->primary_tag eq 'gene' )
            {
                my @gene_id = $feature->get_tag_values('gene');
                print "Gene: $gene_id[0]\n";
                $gene_id = $gene_id[0];
                $found_gene = 1;
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
                                print "Exon ID: $value\n";
                                $exon_id = $2;
                            }
                        }
                    }
                }

                $reverse_strand = $feature->strand;

                if( $feature->strand < 0 )
                {
                    print "complement sequence!!\n";
                }
                else
                {
                    print "NO complement\n";
                }

                print "Exon start: ", $feature->start,
                      " Exon end: ", $feature->end,
                      ".  \nExon string: ",
                      $seq_object->subseq($feature->start, $feature->end),"\n";

                $exon_start = $feature->start;
                $exon_end = $feature->end;

                $sequence = $seq_object->subseq($feature->start,$feature->end);

                if( $found_gene )
                {
                    print OUTPUT ">$gene_id|$exon_id|$exon_start|$exon_end|$reverse_strand\n$sequence\n";
                    print "ID: >$gene_id|$exon_id|$exon_start|$exon_end|$reverse_strand\nSequence: $sequence\n";
                }
            }
        }


        print STDOUT "\n\nFile Record Break, hit 'Ctrl-D' to continue....\t";
        <STDIN>;
    }

    print "That is all now, have a spelndid day!!\n";
}


#
# End of File.
#













#
# End of File.
#
