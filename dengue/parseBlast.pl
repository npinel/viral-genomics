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
use Bio::SearchIO;

my $usage = "\nthe proper command is: perl blast2xmlPlasmid.pl --report=reportname --pl=plasmid \n\n".
    "\treport  => file name (and path)\n".
    "\tplasmid => name of plasmid molecule (for retrieving sequence length)\n\n";

my ($blastreport,$sequences,$outfile,$type);
GetOptions('blast_report=s' => \$blastreport,
           'sequences=s' => \$sequences,
	   'out_file=s' => \$outfile,
	   'type=s' => \$type,);

die $usage unless (defined($blastreport));

my $seqbank = &loadSequences($sequences);

my $blastin = new Bio::SearchIO(-format => 'blast',
				-file => $blastreport);

my $count = 0;
open my $out, '>', $outfile;
while (my $result = $blastin->next_result) {

    my $query   = $result->query_name;
    my $query_length = $result->query_length;
    my $qdesc = $result->query_description;
    
  HIT:while(my $hit = $result->next_hit) {
      my $hit_name = $hit->accession;

      if (!$type || $$seqbank{$hit_name}{'type'} == $type) { # if type is defined, output only sequences for the requested DENV type 
	HSP:while(my $hsp = $hit->next_hsp) {
	    if ($hsp->percent_identity >= 25) {
		my $description = $hit->description;
		my @coords = $hsp->range('hit');
		my $length = $coords[1] - $coords[0] + 1;
		my $subseq = substr($$seqbank{$hit_name}{'seq'}, $coords[0] - 1, $length);
		$hit_name .= "\_$$seqbank{$hit_name}{'type'}";
		print $out qq(\>$hit_name\n$subseq\n\n);
		$count++;
	    }   
	}
      }
  } # end of the while hit loop
} # end of the unless loop
undef($blastin);
close($out);

print qq(created output fasta file with $count sequences\n\n);

exit;

sub loadSequences {
    my $file = shift;

    my %hash = ();
    my $seqin = Bio::SeqIO->new('-file' => $file,
				'-format' => 'fasta');

	while (my $seq_obj = $seqin->next_seq()) {
	    my $seq = $seq_obj->seq();
	    my $id = $seq_obj->display_id;
	    $id =~ /gi\|(.+)\|gb\|(.+)\|/;
	    my $gi = $1; my $gb = $2;

	    my $desc = $seq_obj->description();
	    $desc =~ /^\/(\d{1})\//;

	    $hash{$gb} = {'seq' => $seq,
			  'type' => $1,};
	    
	}

    return \%hash;
}
