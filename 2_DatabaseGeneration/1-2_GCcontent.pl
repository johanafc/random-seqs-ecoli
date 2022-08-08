#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
#use List::MoreUtils qw(uniq);
#use Array::Utils qw(:all);

my @OUT;
my @allSeq;
my %hashSeq;
my %PrFR;
my %PrFName;
my $name = $ARGV[1];
print  "$name\n";
open my $INPUT, "< $ARGV[0]";
open my $OUT, "> GCcontent-$name.fa";

while (my $lineaSeq = <$INPUT>){#2
chomp $lineaSeq;
push @allSeq, $lineaSeq;
}#2

#-------------------------------------------------------------------
#print Dumper \@allSeq;

foreach my $singleSeq (@allSeq){
	my $count = () = $singleSeq =~ /C|G/g;
#	print "$count\n";
	my $seqLen = length($singleSeq);
	my $GC = sprintf("%.3f", $count/$seqLen*100);
	print $OUT "$GC\n";
}

close $INPUT;

print "\nGC content calculated\n";
