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
