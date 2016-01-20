#!/usr/bin/env perl 

use strict;
use warnings;

#use Pod::Usage;
use Getopt::Long;
use Sort::Naturally;
use Term::ANSIColor;
use List::Util qw (sum);

my $usage = "
".colored('TITLE','white bold')."

bayesTraitsWrapper_ML.pl

".colored('SYNOPSIS','white bold')."

A Perl wrapper around the BayesTraits program, running the ML version of the 
analysis across genome-wide gene presence/absence data.

Takes as input:
  (a) Fasta formatted 1/0 gene presence/absence matrix file
  (b) Traits file of the form [GENOMENAME]\\t[1|0]\\n where 1 = trait presence, 0 = trait absence
  (c) Tree file in nexus format
  (d) Note that [GENOMENAME] should correspond between the matrix file, traits file and tree file
  (e) In addition, make sure BayesTraits is in \$PATH, and DEP.command and INDEP.command files are in the cwd

Outputs a file that contains the ML estimate for both the dependent and independent 
models, as well as the calculated LR statistic: 2*(log(DEP) - log(INDEP)) which follows 
a Chi-square distribution with degrees of freedom = the difference in the number of 
parameters between models (in this case, df = 4).

".colored('OPTIONS','white bold')."

--m|matrix  : matrix file [STR] (required)
--t|traits  : traits file [STR] (required)
--r|tree    : tree file in nexus format [STR] (required)
--p|prefix  : output prefix to write to *.BayesTraitsML.table [STR] (default: \"out\")
--k|keep    : keep temp files (default: temp files are deleted)
--h|help    : displays this help and quits

".colored('USAGE','white bold')."

bayesTraitsWrapper_ML_NULL.pl \\
--matrix [MATRIXFILE] \\
--traits [TRAITSFILE] \\
--tree [TREEFILE] \\
--prefix [OUTFILEPREFIX] > [LOGFILE]

See Pagel (1994) Proc R Soc Lond B 255:37 http://www.evolution.rdg.ac.uk/BayesTraits.html for more
info on the BayesTraits program.\n\n";

#############
## get input:
#############

## define variables:
my ($matrixfile,$traitsfile,$treefile,$keep,$help);
## set default values:
my $prefix = "out";
GetOptions ( 
	'matrix|m=s'  => \$matrixfile,
	'traits|t=s'  => \$traitsfile,
	'tree|r=s'    => \$treefile,
	'prefix|p:s'  => \$prefix,
	'keep|k'      => \$keep,
	'help|h'      => \$help,
);

die $usage if $help;
die $usage unless $matrixfile && $traitsfile && $treefile;

####################
## define variables:
####################

my (%binary,%traits,%results);
my (@genomes,@data);
my ($numberOfSites,$skipped,$tested);

##################
## define outfile:
##################

my $outfile = "$prefix.BayesTraitsML.table";

###################################################
## get binary data from fasta-formatted matrixfile:
###################################################

open (my $FAS, $matrixfile) or die "\n\tCannot open $matrixfile: $!\n\n";
chomp (my @FASTA=<$FAS>);
close $FAS;

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

##############
## print info:
##############

print "
\t##~~~~~~~~~~~~~~~~~~~~~##
\t# bayesTraitsWrapper_ML #
\t##~~~~~~~~~~~~~~~~~~~~~##\n\n";

print "\tMatrix file: $matrixfile\n";
print "\tTraits file: $traitsfile\n";
print "\tTree file: $treefile\n";
print "\tOutput written to file: $outfile\n";
print "\tNumber of genomes: ".@genomes."\n";
print "\tNumber of sites: ".$numberOfSites."\n\n\t~~~\n\n";

@binary{@genomes} = @data;

## get traits data:
open (my $TRA, $traitsfile) or die "\n\tCannot open $traitsfile: $!\n\n";
chomp (my @TRAITS=<$TRA>);
close $TRA;

foreach (@TRAITS) {
	my @a = split (/\s+/, $_);
	$traits{$a[0]} = $a[1]; ## key = genome ID; val = trait value (1/0)
}

################
## sanity check:
################

## make sure the taxon names given in traitsfile are the SAME as those in the bins fasta:
foreach (keys %traits) {
	unless (defined $binary{$_}) {
		die "\n## ERROR: Genome ID's in $matrixfile are not the same as those in $traitsfile!\n## ERROR: Make sure they all match up before continuing.\n\n"
	} 
}

###################################
## iterate across all tested sites:
###################################

foreach my $i (0 .. $numberOfSites - 1) {
	
	my $site = $i + 1;
	# print "\tSite: ".$site."\n";
	
	my (@gene_bins,@trait_bins,@genome_names);
	
	foreach (nsort keys %binary) {
		my @row_bins = @{ $binary{$_} };
		
		## each should be ntax in length and in the correct relative order:
		push (@gene_bins, $row_bins[$i]);
		push (@trait_bins, $traits{$_});
		push (@genome_names, $_);
	}
	
	my $sum;
	$sum += $_ for @gene_bins;
	
	if ( ($sum <= 4) or ($sum >= 60) ) { ## 4 -- 60 variable sites represents the ~ 90% middle of the distribution of possible occurrence patterns
		$skipped++;
		
	} else {
		
		## create a (temporary) data file for input into BayesTraits
		open (my $D, ">site_$site.data") or die "\n\t$!\n\n";
		for my $i (0 .. $#genome_names) {
			print $D $genome_names[$i]."\t".$trait_bins[$i]."\t".$gene_bins[$i]."\n";
		}
		close $D;
		
		print "\tSite: ".$site."\n";
		my ($LogL_DEP, $LogL_INDEP);
		
		###################################
		## run BayesTraits DEPENDENT model:
		###################################

		die "\n\t## ERROR: Problem running BayesTraits!\n\n" if (system ("BayesTraits $treefile site_$site.data < DEP.command > data_$site.DEP") != 0 );
		
		## get the logL from the output file:
		open (my $DEP, "data_$site.DEP") or die "\n\tCannot open 'data_$site.DEP': $!\n\n";
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

		die "\n\t## ERROR: Problem running BayesTraits!\n\n" if (system ("BayesTraits $treefile site_$site.data < INDEP.command > data_$site.INDEP") != 0 );
		
		## get the logL from the output file:
		open (my $INDEP, "data_$site.INDEP") or die "\n\tCannot open 'data_$site.INDEP': $!\n\n";
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
		$results{$site} = { 'LogL_DEP'   => $LogL_DEP, 
							'LogL_INDEP' => $LogL_INDEP,
							'LRStat'     => $LRStat
						    };
		
		## remove temp files to save space:
		unlink ("site_$site.data", "site_$site.data.log.txt", "data_$site.DEP", "data_$site.INDEP") or die "\n\tCannot unlink files: $!\n\n" unless $keep;
		
		$tested++;
	}
}

print "\t~~~\n\n";
print "\tNumber of sites tested: ".$tested."\n\n";
print "\tNumber of below-threshold invariant sites skipped: ".$skipped."\n\n";

###################################
## print results to tab-delim file:
###################################

open (my $R, ">".$outfile) or die "\n\t$!\n\n";
print $R "site\tLogL_DEP\tLogL_INDEP\tLRStat\n";
foreach (sort {$a <=> $b} keys %results) {
	my %a = %{ $results{$_} };
	print $R $_."\t".$a{LogL_DEP}."\t".$a{LogL_INDEP}."\t".$a{LRStat}."\n";
}
close $R;

print "\tFinished doing it.\n\n";

__END__
