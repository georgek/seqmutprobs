#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw(max);

print "NAME\tPOS\tCOV\tREF\tA\tC\tG\tT\tD\n";

my $prominence_thresh = 4;
sub filter {
    my (%nuc) = @_;
    my %filterednuc;

    my $fwd_max = max $nuc{"A"},$nuc{"C"},$nuc{"G"},$nuc{"T"};
    my $rev_max = max $nuc{"a"},$nuc{"c"},$nuc{"g"},$nuc{"t"};
    my $fwd_sum = $nuc{"A"} + $nuc{"C"} + $nuc{"G"} + $nuc{"T"};
    my $rev_sum = $nuc{"a"} + $nuc{"c"} + $nuc{"g"} + $nuc{"t"};
    my $fwd_prominence, my $rev_prominence;
    if ($fwd_sum > $prominence_thresh) {
        $fwd_prominence = $fwd_max / $fwd_sum;
    } else {
        $fwd_prominence = 0;
    }
    if ($rev_sum > $prominence_thresh) {
        $rev_prominence = $rev_max / $rev_sum;
    } else {
        $rev_prominence = 0;
    }
    my $ratio;
    if ($fwd_prominence > $rev_prominence) {
        $ratio = $fwd_sum > 0 ? $rev_sum/$fwd_sum : 0;
        $filterednuc{"A"} = int($nuc{"A"} * (1 + $ratio));
        $filterednuc{"C"} = int($nuc{"C"} * (1 + $ratio));
        $filterednuc{"G"} = int($nuc{"G"} * (1 + $ratio));
        $filterednuc{"T"} = int($nuc{"T"} * (1 + $ratio));
    } else {
        $ratio = $rev_sum > 0 ? $fwd_sum/$rev_sum : 0;
        $filterednuc{"A"} = int($nuc{"a"} * (1 + $ratio));
        $filterednuc{"C"} = int($nuc{"c"} * (1 + $ratio));
        $filterednuc{"G"} = int($nuc{"g"} * (1 + $ratio));
        $filterednuc{"T"} = int($nuc{"t"} * (1 + $ratio));
    }

    return %filterednuc;
}

while(<>) {
    my @list = split /\t/, $_;
    my $pileupnucleotides = $list[4];
    my %nuc;
    $nuc{"A"} = 0;
    $nuc{"C"} = 0;
    $nuc{"G"} = 0;
    $nuc{"T"} = 0;
    $nuc{"a"} = 0;
    $nuc{"c"} = 0;
    $nuc{"g"} = 0;
    $nuc{"t"} = 0;
    $nuc{"*"} = 0;
    while (length $pileupnucleotides) {
        my @matches = ($pileupnucleotides  =~ m/\.|,|\*|\+[0-9]+[ACGTNacgtn]+|-[0-9]+[ACGTNacgtn]+|[ACGTNacgtn]/g);
        my %ins = ();
        $pileupnucleotides = "";
        foreach my $m (@matches) {
            if ($m eq "." || $m eq ",") {
                $nuc{$list[2]}++;
            } elsif ($m =~ m/\+/) {
                #print "insertion here: " . $m . "\n";
                my @number = ($m =~ m/[0-9]+/g);
                my @sequence = ($m =~ m/[ACGTNacgtn]+/g);
                $pileupnucleotides .= substr $sequence[0], int($number[0]);
            } elsif ($m =~ m/\-/) {
                #print "deletion here: " . $m . "\n";
                my @number = ($m =~ m/[0-9]+/g);
                my @sequence = ($m =~ m/[ACGTNacgtn]+/g);
                $pileupnucleotides .= substr $sequence[0], int($number[0]);
            } else {
                $nuc{$m}++;
            }
        }
    }
    my %filterednuc = filter(%nuc);
    print $list[0] . "\t" . $list[1] . "\t" . $list[3] . "\t" .  $list[2] . "\t" .
        $filterednuc{"A"} . "\t" . $filterednuc{"C"} . "\t" .
        $filterednuc{"G"} . "\t" . $filterednuc{"T"} . "\t" . $nuc{"*"} . "\n";
}
