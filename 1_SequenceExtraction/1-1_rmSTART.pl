#!/usr/bin/perl
#USAGE: perl 1-12_rmSTART.pl input.fa primers.txt

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
my $name = $ARGV[0];
my $diff = 162;
print  "$name\n";
open my $INPUT, "< $ARGV[0]";
open my $PRIMERS, "< $ARGV[1]";
open my $OUT, "> $ARGV[0].clean.fa"; #con un solo > es como >> en bash (es al reves)
#open my $OUTN, "> $ARGV[0].wrongLen.fa";
#open my $OUTN2, "> $ARGV[0].NOrev.fa";
open my $OUTN3, "> $ARGV[0].unmatched.fa";

while (my $lineaSeq = <$INPUT>){#2
chomp $lineaSeq;
push @allSeq, $lineaSeq;
}#2

my $count = 2;
for (my $i= 0 ; $i <= $#allSeq ; $i += $count){
$hashSeq{"$allSeq[$i-2]"} = "$allSeq[$i-1]";
}

#---------------------------------------------print Dumper \%hashSeq;

while (my $lineaPri = <$PRIMERS>){
chomp $lineaPri;
my @Col_allPri = (split /\s+/, $lineaPri);
#-------------------------------------------------------------------
print Dumper \@Col_allPri;
#[0] forward
#[1] reverse
#[2] name
$PrFR{"$Col_allPri[0]"} = "$Col_allPri[1]";
$PrFName{"$Col_allPri[0]"} = "$Col_allPri[2]";
}
#--------------------------------------------------print Dumper \%PrFR;
#%PrFR: Fr => Rv
#--------------------------------------------------print Dumper \%PrFName;
#%PrFName: Fr => name
#

my $countok = 0;
my $countwlen = 0;
my $countnorv = 0;
my $countnofw = 0;
my $countdump = 0;
foreach my $keySeq (keys %hashSeq){
	my $ValueSeq_Sequence = $hashSeq{$keySeq};
		foreach my $keyPrFR_PrimerForward (keys %PrFR){
			my $ValuePrFR_PrimerReverse = $PrFR{$keyPrFR_PrimerForward};
			my $lenF =  length($keyPrFR_PrimerForward);
			my $lenR =  length($ValuePrFR_PrimerReverse);
			my $lenSeq = length($ValueSeq_Sequence);
			if  ($ValueSeq_Sequence =~ m/$keyPrFR_PrimerForward[A-Z]{$diff}$ValuePrFR_PrimerReverse/){
				my $StartPosition_PrimerForward = index($ValueSeq_Sequence, $keyPrFR_PrimerForward)+18;
				my $StartPosition_PrimerReverse = index($ValueSeq_Sequence, $ValuePrFR_PrimerReverse);
				my $EndPosition_PrimerReverse = $StartPosition_PrimerReverse + $lenR;
				my $sequence = substr($ValueSeq_Sequence,$StartPosition_PrimerForward, $EndPosition_PrimerReverse - $StartPosition_PrimerForward);
				$countok++;
        print $OUT "$keySeq\n$sequence\n";
#				print ".";
			}
			elsif  ($ValueSeq_Sequence =~ m/$keyPrFR_PrimerForward.*$ValuePrFR_PrimerReverse/){
				my $StartPosition_PrimerForward = index($ValueSeq_Sequence, $keyPrFR_PrimerForward)+18;
				my $StartPosition_PrimerReverse = index($ValueSeq_Sequence, $ValuePrFR_PrimerReverse);
				my $EndPosition_PrimerReverse = $StartPosition_PrimerReverse + $lenR;
				my $sequence = substr($ValueSeq_Sequence,$StartPosition_PrimerForward, $EndPosition_PrimerReverse - $StartPosition_PrimerForward);
				$countwlen++;
        print $OUT "$keySeq\n$sequence\n";
#				print ".";
			}
			elsif  ($ValueSeq_Sequence =~ m/$keyPrFR_PrimerForward/){
				my $StartPosition_PrimerForward = index($ValueSeq_Sequence, $keyPrFR_PrimerForward)+18;
				my $EndPosition_PrimerForward = $StartPosition_PrimerForward + $lenF;
				my $sequence = substr($ValueSeq_Sequence,$StartPosition_PrimerForward);
				$countnorv++;
        print $OUTN3 "$keySeq\n$sequence\n";
#				print ".";
			}
			elsif ($ValueSeq_Sequence =~ m/$ValuePrFR_PrimerReverse/){
        my $Index_PrimerReverse = index($ValueSeq_Sequence, $ValuePrFR_PrimerReverse)+$lenR;
				my $sequence = substr($ValueSeq_Sequence,0,$Index_PrimerReverse);
				$countnofw++;
        print $OUTN3 "$keySeq\n$sequence\n";
#        print ".";
			}else{
				$countdump++;
				print $OUTN3 "$keySeq\n$ValueSeq_Sequence\n";
        print ".";
			}
		}
}

my $totalacc = $countok + $countwlen;
my $totalng = $countnorv + $countnofw + $countdump;
my $totaltotal = $totalng + $totalacc;
my $perok = $totalacc / $totaltotal *100;
print "\n\n$totaltotal TotalInput\n$totalacc Accepted\n$countok Correct(150)\n$countwlen WrongLength\n$countnorv NoReversePrimerFound\n$countnofw NoForwardPrimerFound\n$countdump Dumped\n$totalng Rejected\n$perok PercentageAccepted\n\n";

close $INPUT;
