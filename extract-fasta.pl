#!/usr/local/bin/perl
use Bio::Perl;
use Bio::SeqIO;


#read exon ids from a file and extract exon sequences from another file

#input files

my @exonIds;
my @zebra;
my @fugu;
my @output;

my $exon = "SCNExon_25.txt";
my $fugu = "./input_files/fugu.out";
my $zebra = "./input_files/zebra.out";
my $output = "SCNExon_25.fasta";
$i=0;

open (INFILE, "<$exon") or die "Cannot open file: $!\n";
@exonIds = <INFILE>;
close (INFILE);

$in_zebra  = Bio::SeqIO->new ('-format' => 'Fasta' , '-file' => $zebra);
$in_fugu  = Bio::SeqIO->new ('-format' => 'Fasta' , '-file' => $fugu);

#open( OUTFILE, ">$output") or die "Cannot open file: $!\n\n";

$out  = Bio::SeqIO->new ('-format' => 'Fasta', '-file' => ">$output");

#zebra

while ($seqobj = $in_zebra->next_seq())
{

	foreach $exonId (@exonIds)
	{
		@ids = split(/\t/, $exonId);
	#while ($seqobj = $in_zebra->next_seq())
	#{	
print "ids".$ids[0]."\n";
		if (($ids[0] eq $seqobj->id)||($ids[1] eq $seqobj->id))
		{
			$out -> write_seq($seqobj);
		}
	}
}

while ($seqobj = $in_fugu->next_seq())
{
        foreach $exonId (@exonIds)
        {
                @ids = split(/\t/, $exonId);
  
      	 	if (($ids[0] eq $seqobj->id)||($ids[1] eq $seqobj->id))
         	{
        		$out -> write_seq($seqobj);
         	}
#	$i++;
#	print OUTFILE "/n the".$i."pair of exons: /n";
	}
}

#close(OUTFILE);

