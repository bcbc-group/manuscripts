#!/usr/bin/perl

=head1 NAME

 snp_sliding_window.pl
 Tool to average genotype in low coverage sequencing data

=cut

=head1 SYPNOSIS

 snp_sliding_win.pl [-h] -i <input_vcf_file> -s <window_size> 

=head2 I<Flags:>

=over

=item -i

B<input_vcf_file>           vcf input file (mandatory)

=item -s

B<window_size>              window size to use to average

=item -h

B<help>                   print the help

=back

=cut

=head1 DESCRIPTION

 This script averages the genotype in a region
 
=cut

=head1 AUTHORS

srs57@cornell.edu

=cut

=head1 METHODS

 snp_sliding_window.pl

=cut

use strict;
use warnings;
use Getopt::Std;
use Data::Dumper;
use Math::Round;

our ($opt_i, $opt_s, $opt_h);
getopts("i:s:b:h");

if (!$opt_i && !$opt_s && !$opt_h) {
    print "There are n\'t any tags. Print help\n\n";
    help();
}

if ($opt_h) {
    help();
}

## Get the arguments and check them
my $in = $opt_i || 
    die("ARGUMENTS ERROR: -i <inputvcf> option was not supplied.\n");

my $window_len = $opt_s || 10000;

open my $ifh, '<', $in or die("ERROR openning file $in: $!\n");
my $chromo = 'null';
my $window_end = $window_len;
my @rows;

while(<$ifh>) {

    chomp($_);
    my $line = $_;

    if ($line =~ /^##/){
	next;
    }
    
    #print the column headers
    if ($line =~ /^#/){
	print $line . "\n";;
	next;
    }

    my @data = split(/\t/, $_);

    #for genotype columns, make an array of arrays
    #check the coordinate is in the current window  
    if ($data[0] eq $chromo && $data[1] < $window_end){

	my @lines;  

	for (my $i = 9; $i < @data; $i++){
	    my @format = split(/:/, $data[$i]);
	    push @lines, $format[0];
	}

	push @rows, [@lines];
    }

    #if coordinate is not in current window, make calculations for last window
    elsif ($data[0] eq $chromo && $data[1] >= $window_end){

	#pass parameters to the calculate genotype subroutine
	&calc_geno(\@rows,\@data,$window_end,$window_len);

	@rows = ();
        my @lines;

        for (my $i = 9; $i < @data; $i++){
            my @format = split(/:/, $data[$i]);
            push @lines, $format[0];
        }

        push @rows, [@lines];

	#set end point for window                                                                                                                                                                                                       
        $window_end = $data[1] + $window_len;
        $chromo = $data[0];
    }
    
    #Start calculation for a new window
    elsif ($data[0] ne $chromo && $chromo !~ /null/){

        &calc_geno(\@rows,\@data,$window_end,$window_len);

	@rows = ();
	my @lines;

	for (my $i = 9; $i < @data; $i++){
	    my @format = split(/:/, $data[$i]);
	    push @lines, $format[0];
	}

	push @rows, [@lines];

	#set end point for window                                                                                                                          
	$window_end = $data[1] + $window_len;
	$chromo = $data[0];

    }

    #Start first window calculation
    elsif ($chromo =~ /null/){
	my @lines;
	
        for (my $i = 9; $i < @data; $i++){
            my @format = split(/:/, $data[$i]);
            push @lines, $format[0];
        }
	
        push @rows, [@lines];
	
        #set end point for window                                                                                                                                 
        $window_end = $data[1] + $window_len;
	$chromo = $data[0];
    }

    else {
	print STDERR "Error in vcf file\n";
	die;
    }
}

#########################
#####subroutines#########
#########################

sub calc_geno{
    #pass @row and @data
    my @rows2 = @{$_[0]};
    my @data = @{$_[1]};
    my $window_end2 = $_[2];
    my $window_len2 = $_[3]; 
    my @transposed;
    
    #transpose genotype array                                                                                                                                  
    for my $row (@rows2) {
	
	for my $column (0 .. $#{$row}) {
	    push(@{$transposed[$column]}, $row->[$column]);
	}
    }
    
    #find average genotype in transposed array, give genotype coordinate as midpoint of window                                                                 
    my $midpoint = $window_end2 - ($window_len2 / 2);
    $midpoint = round($midpoint);
    print $chromo . "\t" . $midpoint . "\t";

    for (my $i = 2; $i < 8; $i++){
	print $data[$i] . "\t";
    }

    print "GT";
    
    #count the genotypes for each sample in the window                                                                                                         
    for my $new_row (@transposed) {
	my $var = 0;
	my $ref = 0;
	my $het = 0;
	my $missing = 0;
	
	for my $new_col (@{$new_row}) {

	    if ($new_col =~ /1\/1/){
		$var++;
	    }

	    elsif ($new_col =~ /0\/0/){
		$ref++;
	    }
	    
	    elsif ($new_col =~ /1\/0/ || $new_col =~ /0\/1/){
		$het++;
	    }
	    
	    elsif ($new_col =~ /\.\/\./){
		$missing++;
	    }
	}
	
	#Calculate best genotype for window                                                                                                                   
	my $num_snps = $het + $ref + $var;
	if ($num_snps > 10){

	    my $het_freq = $num_snps * 0.4;
	    my $homoz_freq = $num_snps * 0.6;

	    if ($het > $het_freq ){
		print "\t0/1";
	    }
	    
	    elsif ($var > $homoz_freq){
		print "\t1/1";
	    }
	    
	    elsif ($ref > $homoz_freq){
		print "\t0/0";
	    }

	    else {print "\t./.";}

	    ($var, $ref, $het, $missing) = 0;
	}
	
	else {print "\t./.";}
    }
    
    print "\n";

    return;
}

=head2 help

  Usage: help()
  Desc: print help of this script
  Ret: none
  Args: none
  Side_Effects: exit of the script
  Example: if (!@ARGV) {
               help();
           }

=cut

sub help {
  print STDERR <<EOF;
  $0:

    Description:

      snp_sliding_window.pl                                                                                                                                             
      Tool to average genotype in low coverage sequencing data                                                                                                          
    Usage:
       
      snp_window.pl [-h] -i <input_vcf_file> -s <window_size>

    Flags:

      -i <input_vcf_file>          vcf input file (mandatory)
      -s <window_size>             window size (1000 by default)
      -h <help>                    print the help
     
EOF
exit (1);
}

