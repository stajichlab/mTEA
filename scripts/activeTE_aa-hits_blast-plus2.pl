#!/usr/bin/env perl

use strict;
use warnings;
use Bio::AlignIO;
use Bio::SearchIO;
use Bio::Seq;
use Bio::SeqIO;
use Bio::PrimarySeq;
use Bio::Tools::IUPAC;
use FindBin;
use File::Spec;
use Getopt::Long;
use Cwd;
use Cwd 'abs_path';
use Data::Dumper;

my $flank;
my $all;

GetOptions(
  'flank=i'  => \$flank,
  'all'      => \$all,
) or die "Please specify flank length\n";
if (@ARGV != 1) {
  print "Please specify one input file.\n";
  help();
}

my $infile = shift @ARGV;


if (!-f $infile) {
  print "\nSupplied filepath is not valid: $!";
  help();
}

if (!defined $flank) {
    print "\nPlease supply flank length!";
    help();
}

sub help {
  print "
usage:

activeTE_aa-hits.pl -f <int> -a <aa-hits_flank_file> <genome_file>

-f length of sequence flanks
-a keep all output files
";
  exit 1;
}

my $cwd = getcwd();
print "CWD: $cwd\n";
my $perl_path = abs_path($0);

#breakdown directory and filename, create output directory
my ($volume, $in_dir, $filename) = File::Spec->splitpath($infile);
my ($volume2, $perl_dir, $perl_filename) = File::Spec->splitpath($perl_path);

my @fname_fin =  split(/\.fa/, $filename);
#print "@fname_fin\n";
my $fname_fin    = $fname_fin[0];
my $out_dir_name = "aTE-peps_" . $fname_fin;
my $out_path     = File::Spec->catdir($in_dir, $out_dir_name);


if (!-d $out_path) {
  mkdir($out_path) or die "Can't create dir:$out_path $!\n"; 
}
my $tir_path = File::Spec->catpath($volume, $out_path, $fname_fin . "_tirs.fa");
open(my $tir_out, ">", $tir_path);

my $name_path = File::Spec->catpath($volume, $out_path, $fname_fin . "_ori-name_to_short-name.tab");
open(my $name_out, ">", $name_path);

print "opening file\n";
my $in = Bio::SeqIO->new(-file => $infile, -format => 'Fasta');
while ( my $seq_obj = $in->next_seq() ) {
    my $seq_id = $seq_obj->id();
    
    my $start;
    my $end; 
    my $contig;
    my $strand;
    if ( $seq_id =~ m/Sbjct:(.+)[_|\s]Length:.+Location:\(([0-9]*)[_|\s]*-[_|\s]*([0-9]*)\)[_|\s]Direction:(\w+)/ ) {
        $contig = $1;
        $start = $2;
        $end = $3;
        $strand = $4;
    }
    else {
       die "Could not find genomic locus positions in sequence header: $seq_id\n"; 
    }
    
    my $short_seq_id = (split /[_|\s]Length/, $seq_id)[0];
    $short_seq_id =~ s/Sbjct[-|:]//i;
    $short_seq_id =~ s/Query[-|:]//i;
    $short_seq_id =~ s/.fa|.fasta//i;
    $short_seq_id =~ s/_L1_clean//gi;
    $short_seq_id =~ s/[0-9]-//gi;
    $short_seq_id =~ s/-ClassII-TIR-/_/i;
    $short_seq_id =~ s/-/_/g;
    
    if (length($short_seq_id > 30) {
        $short_seq_id = substr($short_seq_id, 0, 30)
    }
    
    print "\nlong: $seq_id\nshort: $short_seq_id\n\n";
    print $name_out "$short_seq_id\t$seq_id\n";
    
    my $seq = $seq_obj->seq();
    if ($strand eq 'minus') {
        $seq =~ tr/ATGCatgc/TACGtacg/;
        $seq = reverse($seq);
    }
    
    my $seq_len = length($seq);
    my $first_path;
    my $last_path;
    
    if (defined $all) {
        $first_path = File::Spec->catpath($volume, $out_path, $short_seq_id . "first_half.txt");
        $last_path  = File::Spec->catpath($volume, $out_path, $short_seq_id . "last_half.txt");
    }
    else {
        $first_path = File::Spec->catpath($volume, $out_path, "first_half.txt");
        $last_path  = File::Spec->catpath($volume, $out_path, "last_half.txt");
    }
    my $first;
    my $last;
        
    if ($seq_len <= ($flank*2)) {
        die "Sequence length for $seq_id of $seq_len is less than required by the given flank length.\n";
    }
    $first = substr($seq, 0, $flank + 10);
    $last = substr($seq, -($flank+10));
    
    #save the two ends as files to use as inputs for a ggsearch search
    open(my $first_out, ">", $first_path) or die "Error creating $first_path. $!\n";
    open(my $last_out, ">", $last_path) or die "Error creating $last_path. $!\n";
    print $first_out ">first\n$first\n";
    print $last_out ">last\n$last\n";
    close($first_out);
    close($last_out);

    #create fasta output filename then call fasta search program
    my $out_opt;
    if (defined $all) {
        $out_opt = File::Spec->catpath($volume, $out_path, $short_seq_id . ".bl2seq.out");
    }
    else {
        $out_opt = File::Spec->catpath($volume, $out_path, $fname_fin . ".bl2seq.out");
    }
        
    system("blastn -task blastn -query $first_path -strand minus -subject $last_path -word_size 5 -gapopen 5 -outfmt 5 -gapextend 2 -penalty -3 -reward 2 -num_alignments 1000 -dust no -out $out_opt");

    my @tir_match_results = match_tirs($seq_obj, $out_opt, 3, $strand, $flank);
    #print "\n";
    #print Dumper(\@tir_match_results);
    #print "\n\n";
    my $c = 1;
    my $element_path = File::Spec->catpath($volume, $out_path, $short_seq_id . "_elements.fa");
    $element_path =~ s/:/_/g;
    open(my $element_out, ">", $element_path ) or die "Error creating $element_path. $!\n";
    
    LINE: foreach my $row_ref (@tir_match_results) {
        #print "\n\n";
        #print Dumper($row_ref);
        #print "\n\n";
        my @entry = @{$row_ref};
        if ($entry[0] == 1) {
            my %matches = %{$entry[1]};
            #print Dumper(\%matches);
            
            
            my $left_index = $matches{"query"}->[0];
            my $right_index = $matches{"hit"}->[0];
            print "left index: $left_index  right index: $right_index\n";
            my $left_tir = $matches{'query'}->[1];
            my $right_tir = $matches{'hit'}->[1];
            $right_tir =~ tr/ATGCatgc/TACGtacg/;
            $right_tir = reverse($right_tir);
            my $element = substr($seq, $left_index, -($right_index));
            
            my $title = $short_seq_id . "_tir" . $c;
            my $ltitle = $short_seq_id . "_ltir" . $c;
            my $rtitle = $short_seq_id . "_rtir" . $c;
            if ($strand eq 'minus') {
                $element =~ tr/ATGCatgc/TACGtacg/;
                $element = reverse($element);
                $left_tir =~ tr/ATGCatgc/TACGtacg/;
                my $new_right_tir = reverse($left_tir);
                $right_tir =~ tr/ATGCatgc/TACGtacg/;
                my $new_left_tir = reverse($right_tir);
                $left_tir = $new_left_tir;
                $right_tir = $new_right_tir;
            }
            my $test_left = substr($element, 0, 10);
            my $test_right = substr($element, -10);
            print "test left: $test_left  test right: $test_right\n";
            print "Left tir: $left_tir  Right tir: $right_tir\n";
            
            print $element_out ">$title\n$element\n";
            print $tir_out ">$ltitle\n$left_tir\n>$rtitle\n$right_tir\n";
            
            $c++;
        }
        else {
            
        }    
    }
}
close($tir_out);
close($name_out);
print "Finished searches!\n";



#--------------------------Subroutines---------------------------------#

sub match_tirs {
  my $self       = shift;    ## seq_obj
  my $input_path = shift;
  my $round      = shift;
  my $strand     = shift;
  my $flank      = shift;
  
  $self->throw("Need Bio::Seq argument")
    unless ref $self && $self->isa('Bio::Seq');
  if (!-f $input_path) {
    die "Supplied filepath is not valid";
  }
  my $seq = $self->seq();
  $seq =~ s/-//g;
  if ($strand eq 'minus') {
        $seq =~ tr/ATGCatgc/TACGtacg/;
        $seq = reverse($seq);
  }
  my $reverse_seq = reverse($seq);
  $reverse_seq =~ tr/ATGCatgc/TACGtacg/;
  
  my $seq_len = length($seq);
  my $seq_name = $self->id();
  my $fa_aln_obj = Bio::SearchIO->new(-format => 'blastxml', -file => $input_path);
  my @result;
  #print "ggsearch results path: $input_path\n";
  #my $last_match_seq;
  #my $last_query_seq;

  #go through the FASTA input object . . down to the HSP
  while (my $result = $fa_aln_obj->next_result) {
    #print Dumper($result);
    #print "numHits: ",$result->num_hits,"\n";
    if ($result->num_hits == 0) {
      push @result, [0, [0]];
      last;
    }
    while (my $hit = $result->next_hit) {
      my $hit_name = $hit->name();
      #print "hit name: $hit_name\n";
      while (my $hsp = $hit->next_hsp) {

        #grab the query, hit, and homology strings
        my $homo_string = $hsp->homology_string;
        my $query_str   = $hsp->query_string;
        my $hit_str     = $hsp->hit_string;
        my $len_homo    = length($homo_string);
        my $len_query   = length($query_str);
        my $len_hit     = length($hit_str);
        my $hit_start = $hsp->start('hit');
        my $hit_end = $hsp->end('hit');
        my $query_start = $hsp->start('query');
        my $query_end = $hsp->end('query');
        
        #print "homology string\n$homo_string\nQuery start: $query_start  hit_start: $hit_start\n";

        #initialize variables
        my $match_len     = 0;
        my $start_pos     = '';
        my $end_pos       = '';
        my $match_query   = '';
        my $match_hit     = '';
        my $match_mis_aln = 0;
        my $total_mis_aln = 0;
        my $last_good     = 0;
        my $hit_pos;
        my $query_pos;
        my $match_cutoff = 8;
        my $first_match = 0;

        #parse homology string, keeping track of match length and mismatches or gaps
        for (my $count = 0 ; $count < length($homo_string) ; $count++) {
          my $homo_char  = substr($homo_string, $count, 1);
          my $query_char = substr($query_str,   $count, 1);
          my $hit_char   = substr($hit_str,     $count, 1);
          if ($round == 1 and $count == 8 and $total_mis_aln >= 5) {
            if ($match_len < 3) {
              $match_len   = 0;
              $start_pos   = '';
              $match_query = '';
              $match_hit   = '';
              $end_pos     = '';
              #print "No TIRs found near start of sequences, resetting counts and ending\n";
              last;
            }
          }
          if ($round == 2) {
            if ($count == 6 and $total_mis_aln >= 4) {
              $match_len   = 0;
              $start_pos   = '';
              $match_query = '';
              $match_hit   = '';
              $end_pos     = '';
              last;
            }
          }
          ## skip any seqs that have 2 or more mismatches in the first 3 bases of the TIR
          #if ($round == 3 or $round == 1) {
          if ($round == 1) {
            if ($count == 3 and $total_mis_aln >= 2) {
              $match_len   = 0;
              $start_pos   = '';
              $match_query = '';
              $match_hit   = '';
              $end_pos     = '';
              last;
            }
          }

          if ($match_len == 0) {
            #if match length equals 0 and position is not a match, continue to next position
            if ($homo_char eq " ") {
                if ($round == 1 or $round == 2) {
                    $total_mis_aln++;
                    next;
                }
                else {
                    if ($first_match == 0) {
                        next;
                    }
                    else {
                        $total_mis_aln++;
                        next;
                    }
                }
            }

            #if position is a match, store info, continue to next position
            elsif ($homo_char eq "|") {
              $start_pos = $count;
              $last_good = $count;
              $match_len++;
              $match_query .= $query_char;
              $match_hit   .= $hit_char;
              if ($first_match == 0) {
                  $first_match = 1;
              }

              #print "Initial match at $start_pos\n";
              next;
            }
          }
          elsif ($match_len >= 1 and $match_len < ($match_cutoff - 1)) {

            #if match length is 1-3 and position is not a match, increment mismatch counter and check if more than one mismatch has occurred
            if ($homo_char eq " ") {
              $match_mis_aln++;
              $total_mis_aln++;

              #allow one mismatch, store info and continue
              if ($match_mis_aln <= 1) {
                $match_len++;
                $last_good = $count;
                $match_query .= $query_char;
                $match_hit   .= $hit_char;

                # print "First Mismatch at $count\n";
                next;
              }

              #more than one mismatch, reset counters and other info, continue
              elsif ($match_mis_aln > 1 and $match_len < 5) {
                $match_len     = 0;
                $start_pos     = '';
                $match_query   = '';
                $match_hit     = '';
                $match_mis_aln = 0;

                #print "Another Mismatch at $count, resetting counts\n";
                next;
              }
              elsif ($match_mis_aln < 3 and $match_len >= 5) {
                $match_len++;
                $last_good = $count;
                $match_query .= $query_char;
                $match_hit   .= $hit_char;

                # print "Another Mismatch at $count\n";
                next;
              }
              elsif ($total_mis_aln >= 3) {
                $match_len   = 0;
                $start_pos   = '';
                $match_query = '';
                $match_hit   = '';
                $end_pos     = '';

                #print "Another Mismatch at $count, resetting counts and ending\n";
                last;
              }
            }

            #position is a match, store info and continue
            elsif ($homo_char eq "|") {
              $last_good = $count;
              $match_len++;
              $match_query .= $query_char;
              $match_hit   .= $hit_char;

              #print "Another match at $count. Length is $match_len\n";
              next;
            }
          }
          elsif ($match_len >= $match_cutoff - 1) {

            #match length is $match_cutoff or higher. If position is not a match, increment mismatch counter and check if more than 2 mismatches have occurred. If a match, continue.
            if ($homo_char eq " ") {
              $match_mis_aln++;
              $total_mis_aln++;

              #mismatches under 3, store info and continue
              if ($match_mis_aln <= 3) {
                $match_len++;
                $last_good = $count;
                $match_query .= $query_char;
                $match_hit   .= $hit_char;

                # print "Another Mismatch at $count, proceeding\n";
                next;
              }

              #mismatches 3 or more, store final info for alignment match and end parsing
              elsif ($match_mis_aln >= 3) {
                $end_pos = $last_good;
                $match_query =~ s/-//g;
                $match_hit   =~ s/-//g;

                #reverse complement the match query sequence
                #$match_hit =~ tr/ATGCatgc/TACGtacg/;
                #$match_hit = reverse($match_hit);
                my $match_query_len = length($match_query);
                my $match_hit_len   = length($match_hit);
                
                                
                #find the position in the full sequence of the hit and query match sequences
                $hit_pos = index(uc($reverse_seq), uc($match_hit), ($flank-($hit_end+2)));
                $query_pos = index(uc($seq), uc($match_query), $query_start-2);

                #print "$seq_name hit_pos:$hit_pos query_pos:$query_pos $match_hit\n\n";
                #store sequence name and the hit and query info
                my %match = ("hit"   => [ $hit_pos, $match_hit, $match_hit_len ],
                    "query" => [ $query_pos, $match_query, $match_query_len ]
                  );
                ## add a catch for $hit_pos or $query_pos == 0 [index returns -1 when $match_hit was not found in $seq]
                push @result, [1, \%match];

                #print "Another Mismatch at $count. Match is long, pushing match info for output\n";
                last;
              }
            }

            #position is a match, store info and continue
            elsif ($homo_char eq "|") {
              $last_good = $count;
              $match_len++;
              $match_query .= $query_char;
              $match_hit   .= $hit_char;
              if ($count == (length($homo_string)-1)) {
                $end_pos = $last_good;
                $match_query =~ s/-//g;
                $match_hit   =~ s/-//g;

                #reverse complement the match query sequence
                #$match_hit =~ tr/ATGCatgc/TACGtacg/;
                #$match_hit = reverse($match_hit);
                my $match_query_len = length($match_query);
                my $match_hit_len   = length($match_hit);
                
                
                #find the position in the full sequence of the hit and query match sequences
                $hit_pos = index(uc($reverse_seq), uc($match_hit), ($flank-($hit_end+2)));
                $query_pos = index(uc($seq), uc($match_query),  $query_start-2);

                #print "$seq_name hit_pos:$hit_pos query_pos:$query_pos $match_hit\n\n";
                #store sequence name and the hit and query info
                my %match = ("hit"   => [ $hit_pos, $match_hit, $match_hit_len ],
                    "query" => [ $query_pos, $match_query, $match_query_len ]
                  );
                ## add a catch for $hit_pos or $query_pos == 0 [index returns -1 when $match_hit was not found in $seq]
                push @result, [1, \%match];

                #print "Another Mismatch at $count. Match is long, pushing match info for output\n";
                last;
              }
              #print "Another match at $count. Length is $match_len\n";
              next;
            }
          }
        }

        #add in check to see if TIRs were found.
        if ($end_pos eq '') {
          push @result, [0, [0]];
        }
      }
    }
  }
  return (@result);
}


