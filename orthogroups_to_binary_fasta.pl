#!/usr/bin/env perl

## reubwn Nov 2017

use strict;
use warnings;
use Getopt::Long;
use Sort::Naturally;

my $usage = "
Reformat Orthogroups.txt file to FASTA 1/0 matrix.
Requires orthogroup sequence IDs to be in the OrthoMCL-style format: \"genomeID\|proteinID\" (see README) 

OPTIONS:
-i|--in   : Orthogroups.txt file [required]
-o|--out  : output filename [default: Orthogroups.fasta]
-h|--help : shows this message

USAGE: orthogroups2fasta.pl -i Orthogroups.txt
\n";

my ($infile, $help);
my $outfile = "Orthogroups.fasta";

GetOptions (
	'i|in=s'  => \$infile,
	'o|out:s' => \$outfile,
	'h|help'  => \$help,
);

die $usage if $help;
die $usage unless $infile;

my (%tgenomes, %groups);

## parse Orthogroups.txt
open (my $INFILE, $infile) or die $!;
while (<$INFILE>) {
	chomp;
	my @a = split m/:\s+/;
	my @b = split (m/\s+/, $a[1]);
	my @c = map { s/\|.+//r } @b;
	my %sgenomes;
	@sgenomes{@c} = (); ## get non-redundant GIDS per OG
	@tgenomes{@c} = (); ## get all GIDs
	$groups{$a[0]} = \@{ [ nsort keys %sgenomes ] }; ## key= OGNAME; val= @[ GIDS ]
	#print qq{@{ [ nsort keys %sgenomes ] }}."\n"; ## uncomment to see non-redundant GID membership
}
close $INFILE;

print STDERR "[INFO] Number of OGs: ".(keys %groups)."\n";
print STDERR "[INFO] Number of genomes: ".(keys %tgenomes)."\n";
print STDERR "[INFO] Writing output to \'$outfile\'\n";

## reformat to fasta
open (my $OUTFILE, ">$outfile") or die $!;
foreach my $GID (nsort keys %tgenomes) {
	print $OUTFILE ">$GID\n";
	foreach (nsort keys %groups) {
		if (grep { $GID eq $_} @{ $groups{$_} }) {
			print $OUTFILE "1";
		} else {
			print $OUTFILE "0";
		}
	}
	print $OUTFILE "\n";
}
close $OUTFILE;

print STDERR "[INFO] Finished on ".`date`."\n";
