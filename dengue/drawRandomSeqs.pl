#!/usr/bin/perl

=head1 NAME

=head1 VERSION

=head1 DESCRIPTION

=head2 Overview

Use this file to parse blast reports in pairwise format (default, 
option -m 0 from command line), and produce an xml file for CGview,
with phylogenetically-based colors.

The script was originally written for the worm symbiont genome
project, and as such some of the variables and targets may be 
biased to values therein relevant.

=head2 Constructor and initialization

=head2 Outputs

=cut

use strict;
use Getopt::Long;
use Bio::SeqIO;

my ($sequences,$outfile,$count);
GetOptions('sequences=s' => \$sequences,
	   'out_file=s' => \$outfile,
	   'count=s' => \$count,);

my $seqbank = &loadSequences($sequences);

open my $out, '>', $outfile;
my $total = scalar keys %{$seqbank};
if ($total == $count) {
    &printAll($seqbank);
} else {
    &pickSeq($seqbank,$count);
}

exit;

#####

sub loadSequences {
    my $file = shift;

    my %hash = ();
    my $seqin = Bio::SeqIO->new('-file' => $file,
				'-format' => 'fasta');

    my $cnt = 0;
    while (my $seq_obj = $seqin->next_seq()) {
	my $seq = $seq_obj->seq();
	my $id = $seq_obj->display_id;
	$id =~ /gi\|(.+)\|gb\|(.+)\|/;
	my $gi = $1; my $gb = $2;
	
	my $desc = $seq_obj->description();
	$desc =~ /^\/(\d{1})\//;
	
	$hash{++$cnt} = {'acc' => $id,
			 'seq' => $seq,};
	
    }
    
    return \%hash;
}

sub printAll {
    my $hash = shift;
    
    my $tally = 0;
    foreach my $k (keys %{$hash}) {
	print $out qq(\>$$hash{$k}{'acc'}\n$$hash{$k}{'seq'}\n\n);
	++$tally;
    }

    &reportCount($tally);
}

sub pickSeq {
    my ($hash,$cnt,$tally) = @_;

    my $num = int(rand($total));
    print $num,"\n";
    if ($$hash{$num}) {
	print $out qq(\>$$hash{$num}{'acc'}\n$$hash{$num}{'seq'}\n\n);
	delete($$hash{$num});
	$tally++;
	&reportCount($tally) if ($cnt == $tally);
    }

    &pickSeq($hash,$cnt,$tally);    
}

sub reportCount {
    my $total = shift;
    print qq(created output fasta file with $total sequences\n\n);
    exit;
}
