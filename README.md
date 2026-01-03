# BioinformaticsPerlScripts
Perl Scripts for Bioinformatics 

These are some Bioinformatics scripts written in Perl that I started working on for Dr. Guoqing Lu in early June 2005.

I am getting rid of a super old 13 Gig drive that I used for a CSV and then SVN repo for development stuff torward the end of college

Dr Lu was attempting to find similiar exons between fugu and zebrafish, and I learned a lot about Bioinformatics from Dr Lu, this was a fun project!

The original README text file has the real instructions!!  

Initial data downloads were from 
ftp://ftp.ensembl.org/pub/fugu-32.2g/data/flatfiles/genbank/
ftp://ftp.ensembl.org/pub/zebrafish-32.5/data/flatfiles/genbank/

Notes from an old readme on Aug 15th 2005 after a meeting with Dr Lu:
when you calculate coverage, if one query has two subject sequences, you
need to either merge them or recalculate the alignment length. For example,
the first top two hits in your output file:

lcl|ENSDARG00000039412:ENSDARE00000424213:36184:37114:1
SINFRUG00000129985:SINFRUE00000771338:7739:8676:1    83.25 203  34    0
238  440  289  491  3e-31  133
lcl|ENSDARG00000039412:ENSDARE00000424213:36184:37114:1
SINFRUG00000129985:SINFRUE00000771338:7739:8676:1    79.43 316  65    0
526  841  607  922  1e-24  111

coverage = (alignment lenght_1 + alignement length_2) / min (query length,
subject length), i.e.,
(203 + 316)/Min [(37114-36184), (8676-7739)]
 
 Guoqing Lu
Assistant Professor
Department of Biology
University of Nebraska at Omaha
Omaha, NE 68182-0040
Tel: 402-5543195
Fax: 402-5543532

More notes:
I checked the output file and several things need to be improved.
1) just went though the first 5 pairs of hits, 2 pairs are duplicated.
Remove the redundancy and leave only the one with higher identity and
coverage.
2) need to calculate coverage for each blast hit
3) sort the blast hits in the descending order with coverage first and
identity second
4) provide an id for each pair of sequences, zf001, ...
e.g.,  >lcl|ENSDARG00000039412:ENSDARE00000424213:36184:37114:1 CHANGED TO
      >zf001|lcl|ENSDARG00000039412:ENSDARE00000424213:36184:37114:1
AND >lcl|SINFRUG00000129985:SINFRUE00000771338:7739:8676:1
      >zf001|lcl|SINFRUG00000129985:SINFRUE00000771338:7739:8676:1

You may want to replace the symbol ":" with "|" as NCBI does with its fasta
file.

Left a message on your cell. I will have a meeting tomorrow till 2pm. Will
be back at around 3pm. We can meet at around 4pm.

Talk to you tomorrow.


Guoqing Lu
Assistant Professor
Department of Biology
University of Nebraska at Omaha
Omaha, NE 68182-0040
Tel: 402-5543195
Fax: 402-5543532


Just want to remind you of checking your script and see if your script
removes duplicatded genes in Fugu

      lcl|ENSDARG00000023927:ENSDARE00000511417:308353:309232:-1
SINFRUG00000157902:SINFRUE00000784743:186382:187553:-1      82.35 204  36
0    569  772  861  1064  4e-27  119
      lcl|ENSDARG00000023927:ENSDARE00000511417:308353:309232:-1
SINFRUG00000157902:SINFRUE00000589925:186382:187568:-1      82.35 204  36
0    569  772  876  1079  4e-27  119

After you do blast, if one query has hits to different subject sequences,
it means there are duplication in the subject genome. In our case, fugu has
duplication.


