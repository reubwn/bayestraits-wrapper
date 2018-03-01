#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Sort::Naturally;
use Scalar::Util qw(looks_like_number);

use Bio::Seq;
use Bio::SeqIO;

my $usage = "
Extracts representative sequences based on user-input LR threshold.
Picks the longest sequence from the OG as representative.

OPTIONS:
-i|--input       : table file from bayesTraitsWrapper_ML analysis [required]
-l|--LRStat      : LRStat threshold to define significance [required]
-b|--binary      : 1/0 presence/absence fasta file [required]
-g|--orthogroups : Orthogroups.txt file [required]
-f|--fasta       : sequences in fasta format, headers must match with Orthogroups.txt [required]
-o|--outdir      : outdir name [default 'output_seqdata']
-h|--help        : shows this message

USAGE: get_significant_seqs_from_LR.pl -i BayesTraitsML.table -l 5 -b Orthogroups.fasta -g Orthogroups.txt -f <(cat /path/to/seqs/*faa)
\n";

my ($inputfile,$lr,$binaryfile,$orthogroupsfile,$fastafile,$help,$debug);
my $outdir = "output_seqdata";

GetOptions (
	'i|input=s'       => \$inputfile,
	'l|LRStat=s'      => \$lr,
	'b|binary=s'      => \$binaryfile,
	'g|orthogroups=s' => \$orthogroupsfile,
	'f|fasta=s'       => \$fastafile,
	'o|outdir:s'      => \$outdir,
	'h|help'          => \$help,
	'd|debug'         => \$debug
);

die $usage if $help;
die $usage unless ($inputfile && $lr && $binaryfile && $orthogroupsfile);

if (-e $outdir && -d $outdir) {
	die "[ERROR] Directory '$outdir' already exists!\n";
} else {
	mkdir "$outdir";
	mkdir "$outdir/fasta";
}

## things we'll need
my ($n);
my (%table_hash,%binary_hash,%orthogroups_hash,%proteins_hash);

## confirm LR threshold
if (looks_like_number($lr)) {
	print STDERR "[INFO] LRStat threshold set to: $lr\n";
} else {
	die "[ERROR] Is LRStat threshold a number?\n";
}

## parse table $input
open (my $TABLE, $inputfile) or die $!;
while (<$TABLE>) {
	chomp;
	next if $. == 1;
	my @F = split (m/\s+/, $_);
	if ($F[3] > $lr) { ## ONLY IF LR > LRSTAT THRESHOLD
		$table_hash{$F[0]} = $F[3]; ## key= site number; val= LRStat
	}
	$n++;
}
close $TABLE;
print STDERR "[INFO] Number of sites parsed from $inputfile: ".commify($n)."\n";
print STDERR "[INFO] Number of sites with LRStat >= $lr: ".commify(scalar(keys %table_hash))." (".percentage(scalar(keys %table_hash),$n)."\%)\n";

## parse binary P/A file
my $in_b = Bio::SeqIO->new ( -file => $binaryfile, -format => 'fasta' );
while ( my $seq_obj = $in_b->next_seq() ) {
	my @a = split (m//, $seq_obj->seq()); ## split 010101 string to array
	$binary_hash{$seq_obj->display_name()} = \@a; ##key= seqid; val=$seq_obj
}
print STDERR "[INFO] Number of sequences in $binaryfile: ".commify(scalar(keys %binary_hash))."\n";
if ($debug) {
	foreach (keys %binary_hash) {
		my @a = @{$binary_hash{$_}};
		print STDOUT "$_\t@a\n"
	}
}

## parse Orthogroups.txt
open (my $GROUPS, $orthogroupsfile) or die $!;
while (<$GROUPS>) {
	chomp;
	my @a = split (m/:\s+/, $_);
	my @b = split (m/\s+/, $a[1]);
	$orthogroups_hash{$.} = \@b; ## key= OG name; val = @{list of proteins in OG}
}
close $GROUPS;
print STDERR "[INFO] Number of OGs in $orthogroupsfile: ".commify(scalar(keys %orthogroups_hash))."\n";

## get sequences from $fasta
my $in_f = Bio::SeqIO->new ( -file => $fastafile, -format => 'fasta' );
while ( my $seq_obj = $in_f->next_seq() ) {
	$proteins_hash{$seq_obj->display_name()} = $seq_obj; ##key= seqid; val=$seq_obj
}
print STDERR "[INFO] Number of sequences in $fastafile: ".commify(scalar(keys %proteins_hash))."\n";

## cycle through %table_hash
open (my $LOG, ">logfile.txt") or die $!;
foreach my $site (sort {$a<=>$b} keys %table_hash) {
	my @participant_names = @{$orthogroups_hash{$site}};
	my $representative = Bio::Seq->new( -length => 0 );

	print "[DEBUG] Site: $site\n[DEBUG] Participants: @participant_names\n[DEBUG] Representative: ".($representative->display_name())."\n" if $debug;

	foreach (@participant_names) {
		print "[DEBUG] ".$proteins_hash{$_}->display_name()."\t".$proteins_hash{$_}->length()."\n" if $debug;
		$representative = $proteins_hash{$_} if $proteins_hash{$_}->length() > $representative->length();
	}

	print $LOG "Site number: $site\n";
	print $LOG "Site LRStat: $table_hash{$site}\n";
	print $LOG "Participating sequences in OG: @participant_names\n";
	print $LOG "Representative sequence ID: ".$representative->display_name()."\n";
	print $LOG "Representative sequence length: ".$representative->length()."\n";
	print $LOG "Representative sequence sequence: ".$representative->seq()."\n";
	print $LOG "~~~\n";

	open (my $FAA, ">site_$site.faa") or die $!;
	print $FAA ">".$representative->display_name()."_".$site."\n".$representative->seq()."\n";
	close $FAA;

}
close $LOG;

## file cleanup
`mv logfile.txt $outdir`;
`mv *faa $outdir/fasta`;
`cat $outdir/fasta/*faa > $outdir/all_seqs.faa`;

print STDERR "[~~~~]\n";
print STDERR "[INFO] Representative sequnces written to $outdir/fasta\n";
print STDERR "[INFO] Multifasta written to $outdir/all_seqs.faa\n";
print STDERR "[INFO] Finished on ".`date`."\n";

################ SUBS

sub percentage {
    my $numerator = $_[0];
    my $denominator = $_[1];
    my $places = "\%.2f"; ## default is two decimal places
    if (exists $_[2]){$places = "\%.".$_[2]."f";};
    my $float = (($numerator / $denominator)*100);
    my $rounded = sprintf("$places",$float);
    return $rounded;
}

sub commify {
    my $text = reverse $_[0];
    $text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
    return scalar reverse $text;
}
