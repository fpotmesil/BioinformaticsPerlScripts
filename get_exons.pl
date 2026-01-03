#!/usr/bin/perl

use strict;
use warnings;


use File::Basename;
use File::Path;
use Boulder::Genbank;
use Bio::Seq;
use Bio::SeqIO;


# My Postscript package
use PostScript::Figure;

#
# to help me see what is going on, 
# this produces a lot of output!!
#
my $debug = 1;

# Ask the user what they want to name the files
print STDOUT "Root name of the output file:\t";
my $input = <STDIN>;
chomp $input;

my $box_color = 'Transparent';


# We work from the bottom up in Postscript, so reverse sort the list of filenames
@ARGV = reverse sort @ARGV;

# Open up the various files we're going to use
(undef,my $path,undef) = fileparse($ARGV[0]);
my $outfile = ">". $path . substr($input,0,24) . ".fas";
open(OUTPUT,$outfile) || die "Can't open the output file $outfile: $! \n";

if ($debug )
{
	my $logfile = ">". $path . substr($input,0,24) . ".log";
	open(LOGFILE,$logfile) || die "Can't open the output file $logfile: $! \n";
}


# declare a new postscript element for the figure
my $genes = new PostScript::Figure;

# scale is 100 bp = 10 coordinate points (100/10 = 10 = scale_factor)
my $scale_factor = 10;
my $offset = 100;

# start at 0,0 and use a height of 20 for the boxes, with a spacing of 30 between genes.
my $x = 0;
my $height = 20;
my $spacing = 30;
my $label_font_size = 14;

# Count keeps track of the number of genes - we use this for vertical spacing
my $count = 1;

# Max_position is used to help draw the scale bar and scale the final
# figure to fit on one page
my $max_position = 0; 

local $/ = "//\n";	# Enable paragraph mode. Slurp in each record

@ARGV = map { /\.(gz|Z)$/ ? "gzip -dc $_ |" : $_  } @ARGV;

foreach my $file (@ARGV) 
{
	open(INPUT,$file) || die "Can't open the input file $file: $! \n";
	my $sequence = Boulder::Genbank->parse(<INPUT>);
	# Pull out the list of CDS and misc_features
	my @array = $sequence->Features->Cds;
	my @misc_features = $sequence->Features->Misc_feature;

	# this is for parsing the coding sequence 
	foreach my $j (0 .. $#array) 
    {
		# Just declare some variables to make 'use strict' happy;
		my @misc_feature_position = ();
		my @feature_length = ();
		my @features = undef;
		my $number_of_features = 0;

		my $position = $array[$j]->Position;
		my $name = $array[$j]->Gene;
		my $protein_id = $array[$j]->Protein_id;

		# Use the gene name, protein name, or both (if present) to name the CDS
        my $seq_id = undef;

		if (defined($protein_id)) 
        {
			$seq_id = $name . "|" . $protein_id;
		}
        else
        {
			$seq_id = $name;
		}
        
		# strip out the uncertainty
		$position =~ s/[<>]//g;


        if( $debug )
        {
		    print LOGFILE "Gene $name $protein_id $count ($position)\n";
        }

		my $rev = 0;
        
		if ($position =~m/complement/)
        {
            $rev = 1;
        }
        
		# this ignores everything but the basepair positions
		$position =~ m/(\d+.*\.\..*\d+)/;
		my $numbers = $1;
		# takes care of the case when the coding sequence has only one exon
		my @exons = ();

		if ($numbers =~ m/,/)
        {
            @exons = split(/,/,$numbers);
        }
		else
        {
            $exons[0] = $numbers;
        }

		# we keep track of the starting position of the exon, and how long it is
		(my $mRNA, my @position, my @exon_length);

		foreach my $i (0 .. $#exons) 
        {
			(my $start,my $end) = split(/\.\./,$exons[$i]);
			my $length = $end - $start +1;
			my $exon_seq = substr($sequence->Sequence,$start-1,$length);
			$mRNA .= $exon_seq;
			
			$position[$i] = $start;
			$exon_length[$i] = $length;
						
		}

		$mRNA =~ s/\/+$//;

		if( $debug )
		{
			print LOGFILE "New Bio::PrimarySeq with seq: $mRNA, seq_id: $seq_id\n";
		}
		
		# make a new sequence object and reverse complement it if necessary
		my $seqobj = Bio::PrimarySeq->new ( -seq => $mRNA,
		                                -id  => $seq_id,
		                            -moltype => 'dna'
		                                              );
		if ($rev == 1) 
        { 
			$rev = undef;
			$seqobj = $seqobj->revcom;
		}

		select OUTPUT;
		print ">" . $seqobj->id() . " mRNA sequence\n";
		print $seqobj->seq();
		print "\n";








		# now parse out misc features for this gene
        if( $debug )
        {
		    print LOGFILE "Misc Feature\tPosition\n";
        }

		foreach my $i (0 .. $#misc_features) 
        {
			my $gene_name = $misc_features[$i]->Gene;
			$gene_name =~ s/^ | $//; #strip off beginning or ending spaces
			my $feature_position = $misc_features[$i]->Position;

			if ($seq_id =~ m/$gene_name/) 
            {
                if( $debug )
                {
				    print LOGFILE "$gene_name\t$feature_position\n";
                }
			}

			# if you remove the =~ m/transmembrane/ clause below,
            #  the script will draw all misc_features
			# this can get messy though...
			if( ($seq_id =~ m/$gene_name/) and 
                    ($misc_features[$i]->Note =~ m/transmembrane/))
            {
				# using the same position parsing strategy as above
				$feature_position =~ s/[<>]//g;
				$feature_position =~ m/(\d+.*\.\..*\d+)/;
				my $feature_numbers = $1;
				my @feature_positions = ();

				if( $feature_numbers =~ m/,/ )
                {
                    @feature_positions = split(/,/,$feature_numbers);
                }
				else
                {
                    $feature_positions[0] = $feature_numbers;
                }

				foreach my $k (0 .. $#feature_positions) 
                {
					(my $start,my $end) = split(/\.\./,$feature_positions[$k]);
					my $length = $end - $start +1;
					$misc_feature_position[$number_of_features] = $start;
					$feature_length[$number_of_features] = $length;

                    if( $debug )
                    {
					    print LOGFILE "Feature $number_of_features starts at $start, length = $length\n";
                    }

					$number_of_features += 1;
				}
			}
		}
	
		# adjust feature number since we added an extra +1 in the loop above
		$number_of_features -= 1;

		# now take care of figure
		# first adjust positions by substracting first position from every position
		# the gene then starts at 0 and extends to the total length of the gene
		# and scale
		
		my $modifier = $position[0];

        if( $debug )
        {
		    print LOGFILE "Modifier: $modifier\n";
	    	print LOGFILE "Positions: $#position\n";
        }

		for my $i (0 .. $#position) 
        {
            if( $debug )
            {
			    print LOGFILE "Position $i: $position[$i]\n";
            }

			$position[$i] = int(($position[$i] - $modifier)/$scale_factor + $offset);

            if( $debug )
            {
			    print LOGFILE "ModPosition $i: $position[$i]\n";
            }

			$exon_length[$i] = int($exon_length[$i]/$scale_factor);

            if( $debug )
            {
			    print LOGFILE"Length: $exon_length[$i]\n";
            }

			my $temp_position = $position[$i] + $exon_length[$i];

			if ($temp_position > $max_position) 
            {
				$max_position = $temp_position;

                if( $debug )
                {
				    print LOGFILE "New maximum width position = $max_position\n";
                }
			}
		}
		
		# now scale the features too
		for my $i (0 .. $number_of_features) 
        {
			$misc_feature_position[$i] = int(($misc_feature_position[$i] - $modifier)/$scale_factor + $offset);
			$feature_length[$i] = int($feature_length[$i]/$scale_factor);
		}
		
		# calculate the vertical position on the page
		my $y1 = $count * $spacing;
		# y2 is for the lines - centered midpoint on boxes
		my $y2 = ($count * $spacing) + ($height/2);
		
		# let's add a name for the structure
		# we could use some improvement here - 
        # I'd like to be able to make it so long names
		# don't overlap the figure itself
		my $label_position;

		if ($label_font_size > $height) 
        {
			# font size is bigger than feature - not a great idea, 
            # so make it level with the box
			$label_position = $y1;
        }
		else
        {
			# center the label next to the box
			$label_position = $y1 + int(($height - $label_font_size)/2);
        }

		$genes->addText(x_coord      => 0,
		                y_coord      => $label_position,
		                font_name    => 'Myriad-Roman',
		                font_size    => $label_font_size,
		                text         => $seq_id,
		                fillcolor    => "Black");
		
		if ($#position > 0)
        {
			for my $i (0 .. $#position) 
            {
				# draw the exon
                if( $debug )
                {
				    print LOGFILE "Box:($position[$i],$y1,$exon_length[$i],$height)\n";
                }

				$genes->addBox(x_coord    => $position[$i],
				               y_coord    => $y1,
				               width      => $exon_length[$i],
				               height     => $height,
				               fillcolor  => $box_color);
				
				# if this isn't the last exon draw the line indicating the intron
				if ($i != $#position) 
                {
					my $x1 = $position[$i] + $exon_length[$i];
					my $x2 = $position[$i+1];

                    if( $debug )
                    {
					    print LOGFILE "Line: ($x1,$y2,$x2,$y2)\n";
                    }

					$genes->addLine(x_coord  => $x1,
					                y_coord  => $y2,
					                x_end    => $x2,
					                y_end    => $y2);
				}
				# otherwise, draw the box and the line indicating the exon and intron
			}
		}
		# takes care of the case where we only have one exon...
		else
        {
			$genes->addBox(x_coord    => $position[0],
                           y_coord    => $y1,
                           width      => $exon_length[0],
                           height     => $height);
		}

		# now draw the features as a series of blue boxes 
        # slightly bigger (+10% on either side)
		# than the exon boxes
		my $height_adjust = int($height * 0.10);
		$y1 = $y1 - $height_adjust;
		my $feature_height = $height + (2 * $height_adjust);

        if( $debug )
        {
		    print LOGFILE "Drawing Features\n";
        }

		for my $i (0 .. $number_of_features) 
        {
			$genes->addBox(x_coord      => $misc_feature_position[$i],
			               y_coord      => $y1,
			               width        => $feature_length[$i],
			               height       => $feature_height,
			               fillcolor    =>"Blue",
			               linecolor    =>"Black",
			               linewidth    => 0.5);

            if( $debug )
            {
			    print LOGFILE "Box: ($misc_feature_position[$i],$y1,$feature_length[$i],$feature_height)\n";
            }
		}

		# increment the counter, clear the variables and move on to the next gene
		$count +=1;
		#@exons = undef;
		#$mRNA = undef;
		#$modifier = undef;
		#@position = undef;
	}
}

# Draw a scale bar at the bottom of the figure, based on the maximum length gene
my $scalebar_center = 20; # this is the y coordinate of the scale bar
my $number_kb = int($max_position/100) - 1;

for my $i (0 .. $number_kb) 
{
	#draw crosshatch
	my $x1 = int(($i * 10 * $scale_factor) + $offset); #every kb
	my $x2 = int((($i + 1) * 10 * $scale_factor) + $offset);
	my $y1 = int($scalebar_center + 4);  # crosshatch is 8 points high
	my $y2 = int($scalebar_center - 4);

	$genes->addLine(x_coord  => $x1,
	                y_coord  => $y1,
	                x_end    => $x1,
	                y_end    => $y2);

	if ($i == $number_kb) 
    {
		# draw the final crosshatch
		$genes->addLine(x_coord  => $x2,
		                y_coord  => $y1,
		                x_end    => $x2,
		                y_end    => $y2);
	}

	# draw piece of scale bar
	$genes->addLine(x_coord  => $x1,
	                y_coord  => $scalebar_center,
	                x_end    => $x2,
	                y_end    => $scalebar_center);

	# finally add the text to indicate kb
	my $text = ($i + 1) . " kb";

	$genes->addText(x_coord      => $x2-10,
	                y_coord      => 1,
	                font_name    => 'Myriad-Roman',
	                font_size    => 12,
	                text         => $text,
	                fillcolor    => "Black");

}

# The size of a standard 8.5 x 11 letter page is 612 x 792 points
# With 1 inch margins this equals 540 by 720
# We'll autoscale the figure if it's larger than this

my ($fig_width,$fig_height) = $genes->getFigureSize;

if ($fig_width > 540 or $fig_height > 720) 
{
	print STDOUT "Autoscaling figure: width = $fig_width, height = $fig_height\n"; 
	$genes->autoScale;
}

my $ps_file = ">". $path . substr($input,0,24) . ".ps";
open (PSFILE, $ps_file) or die "Can not open Postscript output file: $!\n";
print PSFILE $genes->write;
close( PSFILE );




print STDOUT "All finished!\n";

close( LOGFILE );

exit(0);


#
#  End of File.
#
