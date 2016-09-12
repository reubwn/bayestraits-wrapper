#!/usr/bin/env perl

use strict;
use warnings;

#use Pod::Usage;
use Getopt::Long;
use Sort::Naturally;
use Term::ANSIColor;
use List::Util qw (sum shuffle);

my $usage = "
".colored('TITLE','white bold')."

bayesTraitsWrapper_ML_NULL.pl

".colored('SYNOPSIS','white bold')."

A Perl wrapper around the BayesTraits program, to create a ".colored('null distribution','white bold')." of likelihood
ratio statistics across gene presence / absence data by randomly permuting the 1/0 gene data with
respect to the 1/0 trait data. This should break any existing associations and provide a suitable
null distribution against which to test the real data.

Takes as input:
  (a) Fasta formatted 1/0 gene presence/absence matrix file
  (b) Traits file of the form [GENOMENAME]\\t[1|0]\\n where 1 = trait presence, 0 = trait absence
  (c) Tree file in nexus format
  (d) Integer indicating how many times to permute the data for each site

Note that [GENOMENAME] should correspond between the matrix file, traits file and tree file
In addition, make sure BayesTraits is in \$PATH, and DEP.command and INDEP.command files are in the cwd

".colored('OPTIONS','white bold')."

--m|matrix  : matrix file [STR] (required)
--t|traits  : traits file [STR] (required)
--r|tree    : tree file in nexus format [STR] (required)
--p|prefix  : output prefix to write to *.BayesTraitsML.NULL.table [STR] (default: \"out\")
--n|nperms  : number of iterations to run for each orthologous group [INT] (default: 1)
--s|shuffle : type of permutation:
                'p|presence' = permute presence/absence data only
                't|traits'   = permute traits data only
                'b|both'     = permute both [STR] (default: both)
--l|lower   : skip sites (OGs) with <= this number of members [INT] (default: no lower limit)
--u|upper   : skip sites (Ogs) with >= this number of members [INT] (default: no upper limit)
--k|keep    : keep temp files (default: temp files are deleted)
--h|help    : displays this help and quits

".colored('USAGE','white bold')."

bayesTraitsWrapper_ML_NULL.pl \\
--matrix [MATRIXFILE] \\
--traits [TRAITSFILE] \\
--tree [TREEFILE] \\
--prefix [OUTFILEPREFIX] \\
--nperms [INT] \\
--shuffle traits > [LOGFILE]\n\n";

#############
## get input:
#############

## define variables:
my ($matrixfile,$traitsfile,$treefile,$keep,$help,$upper);
## set default values:
my $prefix = "out";
my $nperms = 1;
my $shuffle = "both";
my $lower = 0;
my $flag = 0;
GetOptions (
	'matrix|m=s'  => \$matrixfile,
	'traits|t=s'  => \$traitsfile,
	'tree|r=s'    => \$treefile,
	'prefix|p:s'  => \$prefix,
	'nperms|n:i'  => \$nperms,
	'shuffle|s:s' => \$shuffle,
	'lower|l:i'   => \$lower,
	'upper|u:i'   => \$upper,
	'keep|k'      => \$keep,
	'help|h'      => \$help,
);

die $usage if $help;
die $usage unless $matrixfile && $traitsfile && $treefile;

my (%binary,%traits,%results);
my (@genomes,@data);
my ($numberOfSites,$skipped);

###################################################
## get binary data from fasta-formatted matrixfile:
###################################################

open (my $fh, $matrixfile) or die "\n\tCannot open $matrixfile: $!\n\n";
chomp (my @FASTA=<$fh>);
close $fh;

foreach (@FASTA) {
	my @a;
	if ($_ =~ /^\>/) { ## fasta headers
		$_ =~ s/\>//;
		push (@genomes, $_)
	} else {
		@a = split (//, $_);
		$numberOfSites = @a;
		push (@data, \@a);
	}
}

## define $upper if not given by user:
$upper = scalar @genomes unless $upper;

##################
## define outfile:
##################

my $outfile = "$prefix.BayesTraitsML.NULL.table";

##############
## print info:
##############

print "
\t##~~~~~~~~~~~~~~~~~~~~~~~~~~##
\t# bayesTraitsWrapper_ML_NULL #
\t##~~~~~~~~~~~~~~~~~~~~~~~~~~##\n\n";

print "\tMatrix file: $matrixfile\n";
print "\tTraits file: $traitsfile\n";
print "\tTree file: $treefile\n";
print "\tOutput written to file: $outfile\n";
print "\tNumber of genomes: ".@genomes."\n";
print "\tNumber of sites: $numberOfSites\n";
print "\tNumber of permutations per site: $nperms\n";
print "\tSkip sites with fewer than $lower and more than $upper members\n";

###################
## get traits data:
###################

@binary{@genomes} = @data;

open (my $TRA, $traitsfile) or die "\n\tCannot open $traitsfile: $!\n\n";
chomp (my @TRAITS = <$TRA>);
close $TRA;

foreach (@TRAITS) {
	my @a = split (/\s+/, $_);
	$traits{$a[0]} = $a[1]; ## key = genome ID; val = trait value (1/0)
}

## sanity check - make sure the taxon names given in traitsfile are the SAME as those in the bins fasta:
foreach (keys %traits) {
	unless (defined $binary{$_}) {
		die "\n## ERROR: Genome ID's in $matrixfile are not the same as those in $traitsfile!\n## ERROR: Make sure they all match up before continuing.\n\n"
	}
}

foreach my $i (0 .. $numberOfSites - 1) {

	my $site = $i + 1;

	my @gene_bins;
	my @trait_bins;
	my @genome_names;

	foreach (nsort keys %binary) {
		my @row_bins = @{ $binary{$_} };

		## each should be ntax in length and in the correct relative order:
		push (@gene_bins, $row_bins[$i]);
		push (@trait_bins, $traits{$_});
		push (@genome_names, $_);
	}

	my $sum;
	$sum += $_ for @gene_bins;

	#######################################################################
	## don't test sites with fewer than $lower or more than $upper members:
	#######################################################################

	if ( ($sum <= $lower) or ($sum >= $upper) ) {
		## count number of skipped sites:
		$skipped++;

	} else {

		## do the calculation nperms times:
		for my $j ( 1 .. $nperms ) {

			## implement the shuffle to break any existing association between gene 1/0 and trait 1/0:
			my @shuff_geneBins = shuffle (@gene_bins);
			my @shuff_traitBins = shuffle (@trait_bins);

			######################
			## implement different shuffles depending on --shuffle flag:
			######################

			## shuffle presence / absence data only
			if ( ($shuffle eq 'p') or ($shuffle eq 'presence') ){
				print "\tShuffle mode: presence\n\n\t~~~\n\n" if $flag == 0;
				open (my $D, ">site_$site.$j.data") or die "\n\t$!\n\n";
					for my $i (0 .. $#genome_names) {
					print $D $genome_names[$i]."\t".$trait_bins[$i]."\t".$shuff_geneBins[$i]."\n";
				}
			close $D;
			## shuffle trait data only
			} elsif ( ($shuffle eq 't') or ($shuffle eq 'traits') ){
				print "\tShuffle mode: traits\n\n\t~~~\n\n" if $flag == 0;
				open (my $D, ">site_$site.$j.data") or die "\n\t$!\n\n";
					for my $i (0 .. $#genome_names) {
					print $D $genome_names[$i]."\t".$shuff_traitBins[$i]."\t".$gene_bins[$i]."\n";
				}
			## shuffle both p/a and trait data:
			} elsif ( ($shuffle eq 'b') or ($shuffle eq 'both') ){
				print "\tShuffle mode: both\n\n\t~~~\n\n" if $flag == 0;
				open (my $D, ">site_$site.$j.data") or die "\n\t$!\n\n";
					for my $i (0 .. $#genome_names) {
					print $D $genome_names[$i]."\t".$shuff_traitBins[$i]."\t".$shuff_geneBins[$i]."\n";
				}
			} else {
				die "## ERROR: Unrecognised option for shuffle: $shuffle\n\n";
			}

			print "\tSite: $site.$j\n";
			my ($LogL_DEP, $LogL_INDEP);

			###################################
			## run BayesTraits DEPENDENT model:
			###################################

			die "\n\t## ERROR: Problem running BayesTraits!\n\n" if (system ("BayesTraitsV2 $treefile site_$site.$j.data < DEP.command > data_$site.$j.DEP") != 0 );

			## get the logL from the output file:
			open (my $DEP, "data_$site.$j.DEP") or die "\n\t## ERROR: Cannot open 'data_$site.$j.DEP': $!\n\n";
			while (<$DEP>) {
				if ($_ =~ /Lh/) {
					my @a = split (/\s+/, <$DEP>); ## the <L> reads the next line of L
					$LogL_DEP = $a[1];
				}
			}
			close $DEP;

			#####################################
			## run BayesTraits INDEPENDENT model:
			#####################################

			die "\n\t## ERROR: Problem running BayesTraits!\n\n" if (system ("BayesTraitsV2 $treefile site_$site.$j.data < INDEP.command > data_$site.$j.INDEP") != 0 );

			## get the logL from the output file:
			open (my $INDEP, "data_$site.$j.INDEP") or die "\n\t## ERROR: Cannot open 'data_$site.$j.INDEP': $!\n\n";
			while (<$INDEP>) {
				if ($_ =~ /Lh/) {
					my @a = split (/\s+/, <$INDEP>); ## the <L> reads the next line of L
					$LogL_INDEP = $a[1];
				}
			}
			close $INDEP;

			##################################################
			## calculate LR statistic as: 2*(Log(D) - Log(I)):
			##################################################

			my $LRStat = 2*($LogL_DEP - $LogL_INDEP);

			print "\t\tLogL(DEP): ".$LogL_DEP."\n";
			print "\t\tLogL(INDEP): ".$LogL_INDEP."\n";
			print "\t\tLikelihood ratio statistic: ".$LRStat."\n\n";

			## push to %results:
			$results{"$site.$j"} = { 'LogL_DEP'   => $LogL_DEP,
									 'LogL_INDEP' => $LogL_INDEP,
									 'LRStat'     => $LRStat };

			## remove temp files to save space:
			unlink ("site_$site.$j.data", "site_$site.$j.data.log.txt", "data_$site.$j.DEP", "data_$site.$j.INDEP") or die "\n\tCannot unlink files: $!\n\n" unless $keep;
		}
	}
}

print "\t~~~\n\n";
print "\tNumber of below-threshold invariant sites skipped: ".$skipped."\n\n";

## print results to tab-delim file:
open (my $R, ">".$outfile) or die "\n\t$!\n\n";
print $R "site\tLogL_DEP\tLogL_INDEP\tLRStat\n";
foreach (sort {$a <=> $b} keys %results) {
	my %a = %{ $results{$_} };
	print $R $_."\t".$a{LogL_DEP}."\t".$a{LogL_INDEP}."\t".$a{LRStat}."\n";
}
close $R;

print "\tFinished.\n\n";

__END__
