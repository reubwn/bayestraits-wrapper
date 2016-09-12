#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Sort::Naturally;

my $usage = "
Adds translate block to newick tree file.

OPTIONS:
-r|--tree   : tree file in nexus format [required]
-p|--prefix : prefix for output tree [\"translate_added\"]
-h|--help   : shows this message

USAGE: add_taxa_translate_to_nexus.pl -r [IN_TREE.nexus] -o [OUT]
\n";

my ($treefile,$help);
my $prefix = "translate_added";

GetOptions (
	'tree|r=s'    => \$treefile,
	'prefix|p:s'  => \$prefix,
	'help|h'      => \$help,
);

die $usage if $help;
die $usage unless $treefile;

my $translation = 1;
my %taxa;

## open treefile
open (my $TREE, $treefile) or die $!;
while (<$TREE>) {
  if ($_ =~ /begin trees/i){
    my $tree_line = <$TREE>;
    if ($tree_line =~ /translate/i) {
      die "File already has translate section\n\n";
    } else {
      my @matches = $tree_line =~ /(\w+)\:\d+/g; ## get taxa names
      @taxa{@matches}=(); ## make non-redundant
    }
  }
}
close $TREE;

## get a number for each taxa 1 .. n
foreach (nsort keys %taxa) {
  $taxa{$_} = $translation;
  $translation++;
}

## open treefile again...
open (my $TREE2, $treefile) or die $!;
open (my $OUT, ">$prefix.nexus") or die $!;
while (<$TREE2>) {
  if ($_ =~ /begin trees/i){
    print $OUT $_;
    print $OUT "\ttranslate\n";

    ## reformat translated taxa and print
    my @a;
    foreach (nsort keys %taxa) {
      push (@a, "\t\t$taxa{$_}\t$_");
    }
    print $OUT join ",\n", @a;
    print $OUT ";\n";

    ## get the tree string and do the substitution
    my $tree_line = <$TREE2>;
    my $check = join "|", keys %taxa;
    $tree_line =~ s/($check)/$taxa{$1}/g;
    print $OUT "\t$tree_line";

  } else { print $OUT $_ };
}
close $TREE2;
close $OUT;

__END__
