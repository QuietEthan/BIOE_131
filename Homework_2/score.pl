#!/usr/bin/env perl -w

die "Usage: $0 <program> <protein.fasta> <codonusage.txt> <alpha> <beta> <gamma> <delta>\n" unless @ARGV == 7;
my ($program, $protfile, $usagefile, $alpha, $beta, $gamma, $delta) = @ARGV;

my $codfile = "codons.fa";

my $t0 = time();
system "$program $protfile $usagefile >$codfile";
my $t1 = time();
my $runtime = $t1 - $t0;

sub readFasta {
    my ($file) = @_;
    open FILE, "<$file";
    my ($name, @name, %seq);
    while (<FILE>) {
	if (/^>(\S+)/) {
	    $name = $1;
	    $seq{$name} = '';
	    push @name, $name;
	} elsif (/(\S+)/) {
	    $seq{$name} .= $1;
	}
    }
    close FILE;
    die "No sequences in $file" unless @name;
    return $seq{$name[0]};  # select first sequence
}

my $protseq = readFasta ($protfile);
my $codseq = readFasta ($codfile);
$codseq =~ s/[tT]/u/g;

open COD, "<$usagefile";
my %codprob;
my %aa;
while (<COD>) {
    while (/([ACGU]{3}) (\S) ([\d\.]+)/g) {
	$codprob{lc $1} = $3;
	$aa{lc $1} = lc $2;
    }
}

if (length($codseq) != 3*length($protseq) + 3) {
    print "0 (length incompatible with translated sequence)\n";
    exit;
}

if ($aa{lc substr($codseq,length($codseq)-3,3)} ne '*') {
    print "0 (does not end with stop codon)\n";
    exit;
}

my $sum_log_q = 0;
my $sum_inv_dist = 0;
for (my $pos = 0; $pos < length($protseq); ++$pos) {
    my $codon = lc substr ($codseq, 3*$pos, 3);
    my $aa = lc substr ($protseq, $pos, 1);
    unless ($aa{$codon} eq $aa) {
	print "0 (doesn't correctly translate: at position $pos, $codon is not a codon for $aa)\n";
	exit;
    }
    $sum_log_q += log($codprob{$codon});
    for ($i = $pos - 1; $i >= 0; --$i) {
	my $prevcod = lc substr ($codseq, 3*$i, 3);
	if ($prevcod eq $codon) {
	    $sum_inv_dist += 1 / ($pos - $i);
	    last;
	}
    }
}

my $rnafold = "echo $codseq | RNAfold -T 23";
warn $rnafold;
open FOLD, "$rnafold|";
my $dummy = <FOLD>;
my $fold_energy = <FOLD>;
my $energy;
if ($fold_energy =~ /^\S+ +\( *([\d\+\-\.]+?)\)/) {
    $energy = $1;
}
close FOLD;

print $alpha*$sum_log_q - $beta*$sum_inv_dist - $gamma*$energy - $delta*$runtime, " = alpha*", $sum_log_q, " + beta*", $sum_inv_dist, " + gamma*", $energy, " - delta*", $runtime, "\n";

