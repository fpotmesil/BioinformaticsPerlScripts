#!/usr/bin/perl

use warnings;
use strict;

use Getopt::Long;

print "SHIT!!\n";


my $database_file;
my $query_file;

GetOptions( "database_file=s"=>\$database_file,
        "query_file=s"=>\$query_file);


if( (!$database_file) || (!$query_file) )
{
    print "You'd better look at the code for input options!!\n";
    exit(-1);
}

print "Database File: $database_file, Query File: $query_file\n";

# Create databases directory
if (-d "databases")
{
    system ("rm -rf databases");
    system ("mkdir databases");
}
else
{
    system ("mkdir databases");
}

# Create output_files directory
if (-d "output_files")
{
    system ("rm -rf output_files");
    system ("mkdir output_files");
}
else
{
    system ("mkdir output_files");
}

print "Starting...\n";

chdir("databases");

print "Formatting the blast database for $database_file\n";
system ("formatdb -p F -i $database_file -o");


print "Running Blast Search against sequence $query_file, Please wait...\n";

system ("blastall -p blastn -d $database_file -i $query_file -e 0.00005 -F F -m 8 -o ../output_files/$database_file-$query_file.out");

print "$database_file-$query_file pairwise blast is done ...\n";



