#!/usr/bin/env perl

use strict;
use warnings;
use Sort::Naturally;

my $translation = 1;
my %taxa;

open (my $TREE, $ARGV[0]) or die $!;
while (<$TREE>) {
  if ($_ =~ /begin trees/i){
#    print $_;
    my $tree_line = <$TREE>;
    #print "Tree: $tree_line\n";
    if ($tree_line =~ /translate/i) {
      die "File already has translate section\n\n";
    } else {
      my @matches = $tree_line =~ /(\w+)\:\d+/g; ## get taxa names
      @taxa{@matches}=(); ## make non-redundant
    }
  }
#  else { print $_; }
}
close $TREE;

foreach (nsort keys %taxa) {
  $taxa{$_} = $translation;
  $translation++;
  #print "$_\t$taxa{$_}\n";
}

open (my $TREE2, $ARGV[0]) or die $!;
open (my $OUT, ">$ARGV[0].translate_added.txt") or die $!;
while (<$TREE2>) {
  if ($_ =~ /begin trees/i){
    print $_;
    print "\ttranslate\n";

    my @a;
    foreach (nsort keys %taxa) {
      push (@a, "\t\t$taxa{$_}\t$_");
      #push (@a, "$_");
      #print "\n\t\t$_\t$taxa{$_}";
    }
    print join ",\n", @a;
    print ";\n";

    my $tree_line = <$TREE2>;
    my $check = join "|", keys %taxa;
    #print $tree_line;
    $tree_line =~ s/($check)/$taxa{$1}/g;
    print "\t$tree_line";

  } else { print $_ };
}
