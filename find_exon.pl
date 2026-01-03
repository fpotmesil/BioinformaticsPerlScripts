#!/usr/local/bin/perl

# find_exon.pl
# This script will filter exons with length < min_length and count the number of exons selected and filtered out
# G. Lu, University of Nebraska at Omaha

use Bio::Perl;
use Bio::SeqIO;
use Getopt::Std;
use Getopt::Long;


if(!$ARGV[0]){
   print STDERR "
  ./find_exo.pl exon_size seqfile;

   Extracts exon sequences with length > s
  exon_size 	: exon size
  seqfile 	: entry file (or stdin)
";
}


 $exon_size=$ARGV[0];

 $infile = $ARGV[1];
 $output = $infile.".out";
 $in  = Bio::SeqIO->new ('-format' => 'Fasta' , '-file' => $infile);
 int i=0;

 open (outfile, ">$output");
 $out  = Bio::SeqIO->new ('-format' => 'Fasta', '-file' => ">$output");

while ($seqobj = $in->next_seq()){

 #$out  = Bio::SeqIO->new ('-format' => 'Fasta', '-file' => ">$output");
	@id=split(/\|/,$seqobj->display_id());
  		$exon_start = @id[1];
  		$exon_end = @id[2];

	if (($exon_end - $exon_start) > $exon_size)
   		{
    		#print outfile $seqobj;
		$out -> write_seq($seqobj);
print i++;
   		}
}

open (STATS, ">stats");
print STATS "the number of exons for ".$infile."\n".i."\n";

exit
