#!/usr/bin/env perl

=head1 NAME

    convert_exonerate_gff_to_std_gff.pl

=head1 SYNOPSIS
 
    convert_exonerate_gff_to_std_gff.pl input_gff output_gff 
        where input_gff is the input exonerate gff file,
              output_gff is the name of the output gff file.

=head1 DESCRIPTION

    This script takes an input gff file from exonerate (<input_gff>) and converts it to
    a more standard gff file (<output_gff>).

=head1 VERSION
  
    Perl script last edited 12-Aug-2013.

=head1 CONTACT

    alc@sanger.ac.uk (Avril Coghlan)

=cut

# 
# Perl script convert_exonerate_gff_to_std_gff.pl
# Written by Avril Coghlan (alc@sanger.ac.uk)
# 12-Aug-13.
# Last edited 12-Aug-2013.
# SCRIPT SYNOPSIS: convert_exonerate_gff_to_std_gff.pl: converts a gff file from exonerate to more standard gff format.
#
#------------------------------------------------------------------#

use lib '/home/zoo/zool2500/tree_of_life/raw/gene_cluster_analysis/CS_genes/chemosensory_genes/';

# CHECK IF THERE ARE THE CORRECT NUMBER OF COMMAND-LINE ARGUMENTS:

use strict;
use warnings;

# xxx
# BEGIN {
#     unshift (@INC, '/nfs/users/nfs_a/alc/Documents/git/helminth_scripts/modules');
# }

#use HelminthGenomeAnalysis::AvrilGffUtils;
use AvrilGffUtils;

my $num_args               = $#ARGV + 1;
if ($num_args != 2)
{
    print "Usage of convert_exonerate_gff_to_std_gff.pl\n\n";
    print "perl convert_exonerate_gff_to_std_gff.pl <input_gff> <output_gff>\n";
    print "where <input_gff> is the input exonerate gff file,\n";
    print "      <output_gff> is the output gff file\n";
    print "For example, >perl convert_exonerate_gff_to_std_gff.pl PTRK_exonerate1 PTRK_exonerate1.gff\n";
    exit;
}

# FIND THE PATH TO THE INPUT EXONERATE GFF FILE:                     

my $input_gff              = $ARGV[0];

# FIND THE NAME OF THE OUTPUT GFF FILE:

my $output_gff             = $ARGV[1];

#------------------------------------------------------------------#

# RUN THE MAIN PART OF THE CODE:

&run_main_program($input_gff,$output_gff);

print STDERR "FINISHED.\n";

#------------------------------------------------------------------#

# RUN THE MAIN PART OF THE CODE:

sub run_main_program
{
    my $input_gff           = $_[0]; # INPUT GFF FILE
    my $output_gff          = $_[1]; # OUTPUT GFF FILE
    my $returnvalue;                 # RETURN VALUE FROM A FUNCTION 
  
   # CONVERT THE $input_gff FILE TO MORE STANDARD GFF FORMAT:
    $returnvalue = HelminthGenomeAnalysis::AvrilGffUtils::convert_exonerate_gff_to_standard_gff($input_gff,$output_gff);
}

#------------------------------------------------------------------#
