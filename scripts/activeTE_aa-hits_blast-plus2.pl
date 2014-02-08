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

activeTE_aa-hits.pl -f <int> -a <aa-hits_flank_file>

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
my $right_tir_path = File::Spec->catpath($volume, $out_path, $fname_fin . "_right-tirs.fa");
open(my $right_tir_out, ">", $right_tir_path);
my $left_tir_path = File::Spec->catpath($volume, $out_path, $fname_fin . "_left-tirs.fa");
open(my $left_tir_out, ">", $left_tir_path);
my $no_tir_path = File::Spec->catpath($volume, $out_path, $fname_fin . "_no-TIRs.txt");
open(my $no_tir_out, ">", $no_tir_path);

print "opening file\n";
my $in = Bio::SeqIO->new(-file => $infile, -format => 'Fasta');
while ( my $seq_obj = $in->next_seq() ) {
    my $seq_id = $seq_obj->id();
    
    my $start;
    my $end; 
    my $contig;
    my $strand;
    my $seq_id_ori = $seq_id;
    if ($seq_id =~ m/Chr_/) {
        $seq_id =~ s/Chr_/Chr/;
    }
    elsif ($seq_id =~ m/chr_/) {
        $seq_id =~ s/chr_/chr/;
    }
    print "Seq id ori: $seq_id_ori\nSeq id: $seq_id\n";
    if ( $seq_id =~ m/hit[0-9]*[_| ](.+)[_| ]([0-9]+)[_| ]([0-9]+)[_| ](\w+)/ ) {
        $contig = $1;
        $start = $2;
        $end = $3;
        $strand = $4;
    }
    else {
       die "Could not find genomic locus positions in sequence header: $seq_id\n"; 
    }
    
    
    my $seq = $seq_obj->seq();
    if ($strand eq 'minus') {
        $seq =~ tr/ATGCatgc/TACGtacg/;
        $seq = reverse($seq);
    }
    
    my $seq_len = length($seq);
    my $half = $seq_len/2.0;
    my $first_path;
    my $last_path;
    
    if (defined $all) {
        $first_path = File::Spec->catpath($volume, $out_path, $seq_id_ori . "first_half.txt");
        $last_path  = File::Spec->catpath($volume, $out_path, $seq_id_ori . "last_half.txt");
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
        $out_opt = File::Spec->catpath($volume, $out_path, $seq_id_ori . ".bl2seq.out");
    }
    else {
        $out_opt = File::Spec->catpath($volume, $out_path, $fname_fin . ".bl2seq.out");
    }
        
    system("blastn -task blastn -query $first_path -strand minus -subject $last_path -word_size 5 -gapopen 5 -outfmt 5 -gapextend 2 -penalty -3 -reward 2 -num_alignments 1000 -dust yes -out $out_opt");

    my @tir_match_results = match_tirs($seq_obj, $out_opt, 1, $strand, $flank);
    #print "\n";
    #print Dumper(\@tir_match_results);
    #print "\n\n";
    my $c = 1;
    my $element_path = File::Spec->catpath($volume, $out_path, $seq_id_ori . "_elements.fa");
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
            
            my $tsd1;
            my $tsd2;
            my $tsd3;
            my $tsd4;
            my $tsd4_catch = 0;
               
            my $left_2bp  = substr($matches{"hit"}->[1],  0, 2);
            my $right_2bp = substr($matches{"query"}->[1], -2);
            my $left_4bp  = substr($matches{"hit"}->[1],  0, 4);
            my $right_4bp = substr($matches{"query"}->[1], -4);
            
            my $left_tsd;
            my $alt_left_tsd;
            my $right_tsd;
            my $alt_right_tsd;
            
            if ($left_index >= 20 and ($seq_len - ($right_index + 20) >= 20)) {
                
                $left_tsd = substr($seq, ($left_index - 20),  20);
                $alt_left_tsd = substr($seq, ($left_index - 19), 20);
                
                $right_tsd = substr($seq, ($right_index), 20);
                $alt_right_tsd = substr($seq, ($right_index)-1, 20);
            }
                
            if ($left_2bp eq $right_2bp) {
                $tsd1 = $left_2bp;
                print "Round 1 TSD search - Yes! $fname_fin\n$tsd1\n";
            }
            if ($left_4bp eq $right_4bp) {
                $tsd2 = $left_4bp;
                print "Round 2 TSD search - Yes! $fname_fin\n$tsd2\n";
            }
            if (defined $left_tsd and defined $right_tsd) {
                for (my $i = 0 ; $i < 19 ; $i++) {
                    if (substr($left_tsd, $i) eq substr($right_tsd, 0, -($i))) {
                        $tsd3 = substr($left_tsd, $i);
                        print "Round 3 TSD search - Yes! $i $fname_fin\n$tsd3\n";
                        last;
                    }
                
                }
                for (my $i = 0 ; $i < 19 ; $i++) {
                    if (substr($alt_left_tsd, $i) eq substr($alt_right_tsd, 0, -($i))) {
                        $tsd4 = substr($alt_left_tsd, $i);
                        print "Round 4 TSD search - Yes! $i $fname_fin\n$tsd4\n";
                        $tsd4_catch++;
                        last;
                    }
                }
            }
            if ($tsd4) {
                if (!$tsd1 and !$tsd2 and !$tsd3) {
                    my $left_tir = substr($matches{'hit'}->[1], 1);
                    my $right_tir = $matches{'query'}->[1];
                    $right_tir =~ tr/ATGCatgc/TACGtacg/;
                    $right_tir = reverse($right_tir);
                    $right_tir = substr($right_tir, 0, -1);
                    $left_index++;
                    $right_index--;
                    my $difference = ($half - $left_index) + (($seq_len - $right_index) - $half);
                    my $element = substr($seq, $left_index, -($right_index));
                    
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
                    
                    print $tir_out ">$seq_id_ori\_ltir$c\_$left_index\-$right_index\-$difference\_tsd4\-$tsd4\n$left_tir\n";
                    print $tir_out ">$seq_id_ori\_rtir$c\_$left_index\-$right_index\-$difference\_tsd4\-$tsd4\n$right_tir\n";

                    print $left_tir_out ">$seq_id_ori\_ltir$c\_$left_index\-$right_index\-$difference\_tsd4\-$tsd4\n$left_tir\n";
                    print $right_tir_out ">$seq_id_ori\_rtir$c\_$left_index\-$right_index\-$difference\_tsd4\-$tsd4\n$right_tir\n";
                    
                    my $title = $seq_id_ori."_tir".$c."_".$left_index."-".$right_index . "-$difference" . "_tsd4-$tsd4";
                    print $element_out ">$title\n$element\n";
                    
                    
                    $c++;
                    next LINE;
                }
                elsif (!$tsd1 and !$tsd2 and $tsd3) {
                    if (length($tsd4) > length($tsd3)+1) {
                        my $left_tir = substr($matches{'hit'}->[1], 1);
                        my $right_tir = $matches{'query'}->[1];
                        $right_tir =~ tr/ATGCatgc/TACGtacg/;
                        $right_tir = reverse($right_tir);
                        $right_tir = substr($right_tir, 0, -1);
                        $left_index++;
                        $right_index--;
                        my $difference = ($half - $left_index) + (($seq_len - $right_index) - $half);
                        my $element = substr($seq, $left_index, -($right_index));
                        
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
                        print $tir_out ">$seq_id_ori\_ltir$c\_$left_index\-$right_index\-$difference\_tsd4\-$tsd4\n$left_tir\n";
                        print $tir_out ">$seq_id_ori\_rtir$c\_$left_index\-$right_index\-$difference\_tsd4\-$tsd4\n$right_tir\n";
                        
                        print $left_tir_out ">$seq_id_ori\_ltir$c\_$left_index\-$right_index\-$difference\_tsd4\-$tsd4\n$left_tir\n";
                        print $right_tir_out ">$seq_id_ori\_rtir$c\_$left_index\-$right_index\-$difference\_tsd4\-$tsd4\n$right_tir\n";
                        
                        my $title = $seq_id_ori."_tir".$c."_".$left_index."-".$right_index . "-$difference" . "_tsd4-$tsd4";
                        print $element_out ">$title\n$element\n";
                        
                        $c++;
                        next LINE;
                    }
                    else {
                        my $left_tir = $matches{'hit'}->[1];
                        my $right_tir = $matches{'query'}->[1];
                        $right_tir =~ tr/ATGCatgc/TACGtacg/;
                        $right_tir = reverse($right_tir);
                        my $difference = ($half - $left_index) + (($seq_len - $right_index) - $half);
                        my $element = substr($seq, $left_index, -($right_index));
                        
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
                        print $tir_out ">$seq_id_ori\_ltir$c\_$left_index\-$right_index\-$difference\_tsd3\-$tsd3\n$left_tir\n";
                        print $tir_out ">$seq_id_ori\_rtir$c\_$left_index\-$right_index\-$difference\_tsd3\-$tsd3\n$right_tir\n";
                        
                        print $left_tir_out ">$seq_id_ori\_ltir$c\_$left_index\-$right_index\-$difference\_tsd3\-$tsd3\n$left_tir\n";
                        print $right_tir_out ">$seq_id_ori\_rtir$c\_$left_index\-$right_index\-$difference\_tsd3\-$tsd3\n$right_tir\n";
                        
                        my $title = $seq_id_ori."_tir".$c."_".$left_index."-".$right_index . "-$difference" . "_tsd3-$tsd3";
                        print $element_out ">$title\n$element\n";
                        
                        $c++;
                        next LINE;
                    }
                }
                elsif ($tsd2) {
                    if (length($tsd4) > length($tsd2)+1) {
                        my $left_tir = substr($matches{'hit'}->[1], 1);
                        my $right_tir = $matches{'query'}->[1];
                        $right_tir =~ tr/ATGCatgc/TACGtacg/;
                        $right_tir = reverse($right_tir);
                        $right_tir = substr($right_tir, 0, -1);
                        $left_index++;
                        $right_index--;
                        my $difference = ($half - $left_index) + (($seq_len - $right_index) - $half);
                        my $element = substr($seq, $left_index, -($right_index));
                        
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
                        
                        print $tir_out ">$seq_id_ori\_ltir$c\_$left_index\-$right_index\-$difference\_tsd4\-$tsd4\n$left_tir\n";
                        print $tir_out ">$seq_id_ori\_rtir$c\_$left_index\-$right_index\-$difference\_tsd4\-$tsd4\n$right_tir\n";
                        
                        print $left_tir_out ">$seq_id_ori\_ltir$c\_$left_index\-$right_index\-$difference\_tsd4\-$tsd4\n$left_tir\n";
                        print $right_tir_out ">$seq_id_ori\_rtir$c\_$left_index\-$right_index\-$difference\_tsd4\-$tsd4\n$right_tir\n";
                        
                        my $title = $seq_id_ori."_tir".$c."_".$left_index."-".$right_index . "-$difference" . "_tsd4-$tsd4";
                        print $element_out ">$title\n$element\n";
                        
                        $c++;
                        next LINE;
                    }
                    else {
                        my $left_tir = substr($matches{'hit'}->[1], 4);
                        my $right_tir = $matches{'query'}->[1];
                        $right_tir =~ tr/ATGCatgc/TACGtacg/;
                        $right_tir = reverse($right_tir);
                        $right_tir = substr($right_tir, 0, -4);
                        $left_index = $left_index + 4;
                        $right_index = $right_index - 4;
                        my $difference = ($half - $left_index) + (($seq_len - $right_index) - $half);
                        my $element = substr($seq, $left_index, -($right_index));
                        
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
                        
                        print $tir_out ">$seq_id_ori\_ltir$c\_$left_index\-$right_index\-$difference\_tsd2\-$tsd2\n$left_tir\n";
                        print $tir_out ">$seq_id_ori\_rtir$c\_$left_index\-$right_index\-$difference\_tsd2\-$tsd2\n$right_tir\n";
                        
                        print $left_tir_out ">$seq_id_ori\_ltir$c\_$left_index\-$right_index\-$difference\_tsd2\-$tsd2\n$left_tir\n";
                        print $right_tir_out ">$seq_id_ori\_rtir$c\_$left_index\-$right_index\-$difference\_tsd2\-$tsd2\n$right_tir\n";
                        
                        my $title = $seq_id_ori."_tir".$c."_".$left_index."-".$right_index . "-$difference" . "_tsd2-$tsd2";
                        print $element_out ">$title\n$element\n";
                        
                        $c++;
                        next LINE;
                    }
                }
                elsif (!$tsd2 and !$tsd3 and $tsd1) {
                    if (length($tsd4) > length($tsd1)) {
                        my $left_tir = substr($matches{'hit'}->[1], 1);
                        my $right_tir = $matches{'query'}->[1];
                        $right_tir =~ tr/ATGCatgc/TACGtacg/;
                        $right_tir = reverse($right_tir);
                        $right_tir = substr($right_tir, 0, -1);
                        $left_index++;
                        $right_index--;
                        my $difference = ($half - $left_index) + (($seq_len - $right_index) - $half);
                        my $element = substr($seq, $left_index, -($right_index));
                        
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
                        
                        print $tir_out ">$seq_id_ori\_ltir$c\_$left_index\-$right_index\-$difference\_tsd4\-$tsd4\n$left_tir\n";
                        print $tir_out ">$seq_id_ori\_rtir$c\_$left_index\-$right_index\-$difference\_tsd4\-$tsd4\n$right_tir\n";
                        
                        print $left_tir_out ">$seq_id_ori\_ltir$c\_$left_index\-$right_index\-$difference\_tsd4\-$tsd4\n$left_tir\n";
                        print $right_tir_out ">$seq_id_ori\_rtir$c\_$left_index\-$right_index\-$difference\_tsd4\-$tsd4\n$right_tir\n";
                        
                        my $title = $seq_id_ori."_tir".$c."_".$left_index."-".$right_index . "-$difference" . "_tsd4-$tsd4";
                        print $element_out ">$title\n$element\n";
                        
                        $c++;
                        next LINE;
                    }
                    else {
                        my $left_tir = substr($matches{'hit'}->[1], 2);
                        my $right_tir = $matches{'query'}->[1];
                        $right_tir =~ tr/ATGCatgc/TACGtacg/;
                        $right_tir = reverse($right_tir);
                        $right_tir = substr($right_tir, 0, -2);
                        $left_index = $left_index + 2;
                        $right_index = $right_index - 2;
                        my $difference = ($half - $left_index) + (($seq_len - $right_index) - $half);
                        my $element = substr($seq, $left_index, -($right_index));
                        
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
                        
                        print $tir_out ">$seq_id_ori\_ltir$c\_$left_index\-$right_index\-$difference\_tsd1\-$tsd1\n$left_tir\n";
                        print $tir_out ">$seq_id_ori\_rtir$c\_$left_index\-$right_index\-$difference\_tsd1\-$tsd1\n$right_tir\n";
                        
                        print $left_tir_out ">$seq_id_ori\_ltir$c\_$left_index\-$right_index\-$difference\_tsd1\-$tsd1\n$left_tir\n";
                        print $right_tir_out ">$seq_id_ori\_rtir$c\_$left_index\-$right_index\-$difference\_tsd1\-$tsd1\n$right_tir\n";
                        
                        my $title = $seq_id_ori."_tir".$c."_".$left_index."-".$right_index . "-$difference" . "_tsd1-$tsd1";
                        print $element_out ">$title\n$element\n";
                        
                        $c++;
                        next LINE;
                    }
                }
            }
            
            if ($tsd1 and $tsd4_catch == 0) {
                if (!$tsd2 and !$tsd3) {
                    my $left_tir = substr($matches{'hit'}->[1], 2);
                    my $right_tir = $matches{'query'}->[1];
                    $right_tir =~ tr/ATGCatgc/TACGtacg/;
                    $right_tir = reverse($right_tir);
                    $right_tir = substr($right_tir, 0, -2);
                    $left_index = $left_index + 2;
                    $right_index = $right_index - 2;
                    my $difference = ($half - $left_index) + (($seq_len - $right_index) - $half);
                    my $element = substr($seq, $left_index, -($right_index));
                    
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
                    
                    print $tir_out ">$seq_id_ori\_ltir$c\_$left_index\-$right_index\-$difference\_tsd1\-$tsd1\n$left_tir\n";
                    print $tir_out ">$seq_id_ori\_rtir$c\_$left_index\-$right_index\-$difference\_tsd1\-$tsd1\n$right_tir\n";
                    
                    print $left_tir_out ">$seq_id_ori\_ltir$c\_$left_index\-$right_index\-$difference\_tsd1\-$tsd1\n$left_tir\n";
                    print $right_tir_out ">$seq_id_ori\_rtir$c\_$left_index\-$right_index\-$difference\_tsd1\-$tsd1\n$right_tir\n";
                    
                    my $title = $seq_id_ori."_tir".$c."_".$left_index."-".$right_index . "-$difference" . "_tsd1-$tsd1";
                    print $element_out ">$title\n$element\n";
                    
                    $c++;
                    next LINE;
                }
                elsif (!$tsd2 and $tsd3) {
                  if ( (length($tsd3) > length($tsd1)) and (substr($tsd3, 0, 2) eq $tsd1) and (substr($tsd3, -2) eq $tsd1)) {
                    my $left_tir = $matches{'hit'}->[1];
                    my $right_tir = $matches{'query'}->[1];
                    $right_tir =~ tr/ATGCatgc/TACGtacg/;
                    $right_tir = reverse($right_tir);
                    my $difference = ($half - $left_index) + (($seq_len - $right_index) - $half);
                    my $element = substr($seq, $left_index, -($right_index));
                    
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
                    
                    print $tir_out ">$seq_id_ori\_ltir$c\_$left_index\-$right_index\-$difference\_tsd3\-$tsd3\n$left_tir\n";
                    print $tir_out ">$seq_id_ori\_rtir$c\_$left_index\-$right_index\-$difference\_tsd3\-$tsd3\n$right_tir\n";
                    
                    print $left_tir_out ">$seq_id_ori\_ltir$c\_$left_index\-$right_index\-$difference\_tsd3\-$tsd3\n$left_tir\n";
                    print $right_tir_out ">$seq_id_ori\_rtir$c\_$left_index\-$right_index\-$difference\_tsd3\-$tsd3\n$right_tir\n";
                    
                    my $title = $seq_id_ori."_tir".$c."_".$left_index."-".$right_index . "-$difference" . "_tsd3-$tsd3";
                    print $element_out ">$title\n$element\n";
                    
                    $c++;
                    next LINE;
                  }
                  else {
                    my $left_tir = substr($matches{'hit'}->[1], 2);
                    my $right_tir = $matches{'query'}->[1];
                    $right_tir =~ tr/ATGCatgc/TACGtacg/;
                    $right_tir = reverse($right_tir);
                    $right_tir = substr($right_tir, 0, -2);
                    $left_index = $left_index + 2;
                    $right_index = $right_index - 2;
                    my $difference = ($half - $left_index) + (($seq_len - $right_index) - $half);
                    my $element = substr($seq, $left_index, -($right_index));
                    
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
                    
                    print $tir_out ">$seq_id_ori\_ltir$c\_$left_index\-$right_index\-$difference\_tsd1\-$tsd1\n$left_tir\n";
                    print $tir_out ">$seq_id_ori\_rtir$c\_$left_index\-$right_index\-$difference\_tsd1\-$tsd1\n$right_tir\n";
                    
                    print $left_tir_out ">$seq_id_ori\_ltir$c\_$left_index\-$right_index\-$difference\_tsd1\-$tsd1\n$left_tir\n";
                    print $right_tir_out ">$seq_id_ori\_rtir$c\_$left_index\-$right_index\-$difference\_tsd1\-$tsd1\n$right_tir\n";
                    
                    my $title = $seq_id_ori."_tir".$c."_".$left_index."-".$right_index . "-$difference" . "_tsd1-$tsd1";
                    print $element_out ">$title\n$element\n";
                    
                    $c++;
                    next LINE;
                  }
                }
                elsif ($tsd2 and !$tsd3 and $tsd4_catch == 0) {
                  if (substr($tsd2, 0, 2) eq $tsd1 and substr($tsd2, -2) eq $tsd1) {
                    my $left_tir = substr($matches{'hit'}->[1], 4);
                    my $right_tir = $matches{'query'}->[1];
                    $right_tir =~ tr/ATGCatgc/TACGtacg/;
                    $right_tir = reverse($right_tir);
                    $right_tir = substr($right_tir, 0, -4);
                    $left_index = $left_index + 4;
                    $right_index = $right_index - 4;
                    my $difference = ($half - $left_index) + (($seq_len - $right_index) - $half);
                    my $element = substr($seq, $left_index, -($right_index));
                    
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
                    
                    print $tir_out ">$seq_id_ori\_ltir$c\_$left_index\-$right_index\-$difference\_tsd2\-$tsd2\n$left_tir\n";
                    print $tir_out ">$seq_id_ori\_rtir$c\_$left_index\-$right_index\-$difference\_tsd2\-$tsd2\n$right_tir\n";
                    
                    print $left_tir_out ">$seq_id_ori\_ltir$c\_$left_index\-$right_index\-$difference\_tsd2\-$tsd2\n$left_tir\n";
                    print $right_tir_out ">$seq_id_ori\_rtir$c\_$left_index\-$right_index\-$difference\_tsd2\-$tsd2\n$right_tir\n";
                    
                    my $title = $seq_id_ori."_tir".$c."_".$left_index."-".$right_index . "-$difference" . "_tsd2-$tsd2";
                    print $element_out ">$title\n$element\n";
                    
                    $c++;
                    next LINE;
                  }
                  else {
                    my $left_tir = substr($matches{'hit'}->[1], 2);
                    my $right_tir = $matches{'query'}->[1];
                    $right_tir =~ tr/ATGCatgc/TACGtacg/;
                    $right_tir = reverse($right_tir);
                    $right_tir = substr($right_tir, 0, -2);
                    $left_index = $left_index + 2;
                    $right_index = $right_index - 2;
                    my $difference = ($half - $left_index) + (($seq_len - $right_index) - $half);
                    my $element = substr($seq, $left_index, -($right_index));
                    
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
                    
                    print $tir_out ">$seq_id_ori\_ltir$c\_$left_index\-$right_index\-$difference\_tsd1\-$tsd1\n$left_tir\n";
                    print $tir_out ">$seq_id_ori\_rtir$c\_$left_index\-$right_index\-$difference\_tsd1\-$tsd1\n$right_tir\n";
                    
                    print $left_tir_out ">$seq_id_ori\_ltir$c\_$left_index\-$right_index\-$difference\_tsd1\-$tsd1\n$left_tir\n";
                    print $right_tir_out ">$seq_id_ori\_rtir$c\_$left_index\-$right_index\-$difference\_tsd1\-$tsd1\n$right_tir\n";
                    
                    my $title = $seq_id_ori."_tir".$c."_".$left_index."-".$right_index . "-$difference" . "_tsd1-$tsd1";
                    print $element_out ">$title\n$element\n";
                    
                    $c++;
                    next LINE;
                  }
                }
                elsif ($tsd2 and $tsd3 and $tsd4_catch == 0) {
                  if ( (substr($tsd3, 0, 2) eq $tsd1) and (substr($tsd3, -2) eq $tsd1)) {
                    my $left_tir = $matches{'hit'}->[1];
                    my $right_tir = $matches{'query'}->[1];
                    $right_tir =~ tr/ATGCatgc/TACGtacg/;
                    $right_tir = reverse($right_tir);
                    my $difference = ($half - $left_index) + (($seq_len - $right_index) - $half);
                    my $element = substr($seq, $left_index, -($right_index));
                    
                    
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
                    
                    print $tir_out ">$seq_id_ori\_ltir$c\_$left_index\-$right_index\-$difference\_tsd3\-$tsd3\n$left_tir\n";
                    print $tir_out ">$seq_id_ori\_rtir$c\_$left_index\-$right_index\-$difference\_tsd3\-$tsd3\n$right_tir\n";
                    
                    print $left_tir_out ">$seq_id_ori\_ltir$c\_$left_index\-$right_index\-$difference\_tsd3\-$tsd3\n$left_tir\n";
                    print $right_tir_out ">$seq_id_ori\_rtir$c\_$left_index\-$right_index\-$difference\_tsd3\-$tsd3\n$right_tir\n";
                    
                    my $title = $seq_id_ori."_tir".$c."_".$left_index."-".$right_index . "-$difference" . "_tsd3-$tsd3";
                    print $element_out ">$title\n$element\n";
                    
                    $c++;
                    next LINE;
                  }
                  elsif ((substr($tsd2, 0, 2) eq $tsd1) and (substr($tsd2, -2) eq $tsd1)) {
                    my $left_tir = substr($matches{'hit'}->[1], 4);
                    my $right_tir = $matches{'query'}->[1];
                    $right_tir =~ tr/ATGCatgc/TACGtacg/;
                    $right_tir = reverse($right_tir);
                    $right_tir = substr($right_tir, 0, -4);
                    $left_index = $left_index + 4;
                    $right_index = $right_index - 4;
                    my $difference = ($half - $left_index) + (($seq_len - $right_index) - $half);
                    my $element = substr($seq, $left_index, -($right_index));
                    
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
                    
                    print $tir_out ">$seq_id_ori\_ltir$c\_$left_index\-$right_index\-$difference\_tsd2\-$tsd2\n$left_tir\n";
                    print $tir_out ">$seq_id_ori\_rtir$c\_$left_index\-$right_index\-$difference\_tsd2\-$tsd2\n$right_tir\n";
                    
                    print $left_tir_out ">$seq_id_ori\_ltir$c\_$left_index\-$right_index\-$difference\_tsd2\-$tsd2\n$left_tir\n";
                    print $right_tir_out ">$seq_id_ori\_rtir$c\_$left_index\-$right_index\-$difference\_tsd2\-$tsd2\n$right_tir\n";
                    
                    my $title = $seq_id_ori."_tir".$c."_".$left_index."-".$right_index . "-$difference" . "_tsd2-$tsd2";
                    print $element_out ">$title\n$element\n";
                    
                    $c++;
                    next LINE;
                  }
                  else {
                    my $left_tir = substr($matches{'hit'}->[1], 2);
                    my $right_tir = $matches{'query'}->[1];
                    $right_tir =~ tr/ATGCatgc/TACGtacg/;
                    $right_tir = reverse($right_tir);
                    $right_tir = substr($right_tir, 0, -2);
                    $left_index = $left_index + 2;
                    $right_index = $right_index - 2;
                    my $difference = ($half - $left_index) + (($seq_len - $right_index) - $half);
                    my $element = substr($seq, $left_index, -($right_index));
                    
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
                    
                    print $tir_out ">$seq_id_ori\_ltir$c\_$left_index\-$right_index\-$difference\_tsd1\-$tsd1\n$left_tir\n";
                    print $tir_out ">$seq_id_ori\_rtir$c\_$left_index\-$right_index\-$difference\_tsd1\-$tsd1\n$right_tir\n";
                    
                    print $left_tir_out ">$seq_id_ori\_ltir$c\_$left_index\-$right_index\-$difference\_tsd1\-$tsd1\n$left_tir\n";
                    print $right_tir_out ">$seq_id_ori\_rtir$c\_$left_index\-$right_index\-$difference\_tsd1\-$tsd1\n$right_tir\n";
                    
                    my $title = $seq_id_ori."_tir".$c."_".$left_index."-".$right_index . "-$difference" . "_tsd1-$tsd1";
                    print $element_out ">$title\n$element\n";
                    
                    $c++;
                    next LINE;
                  }
                }
            }
            elsif ($tsd3 and $tsd4_catch == 0) {
                if (!$tsd2) {
                    my $left_tir = $matches{'hit'}->[1];
                    my $right_tir = $matches{'query'}->[1];
                    $right_tir =~ tr/ATGCatgc/TACGtacg/;
                    $right_tir = reverse($right_tir);
                    my $difference = ($half - $left_index) + (($seq_len - $right_index) - $half);
                    my $element = substr($seq, $left_index, -($right_index));
                    
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
                    
                    print $tir_out ">$seq_id_ori\_ltir$c\_$left_index\-$right_index\-$difference\_tsd3\-$tsd3\n$left_tir\n";
                    print $tir_out ">$seq_id_ori\_rtir$c\_$left_index\-$right_index\-$difference\_tsd3\-$tsd3\n$right_tir\n";
                    
                    print $left_tir_out ">$seq_id_ori\_ltir$c\_$left_index\-$right_index\-$difference\_tsd3\-$tsd3\n$left_tir\n";
                    print $right_tir_out ">$seq_id_ori\_rtir$c\_$left_index\-$right_index\-$difference\_tsd3\-$tsd3\n$right_tir\n";
                    
                    my $title = $seq_id_ori."_tir".$c."_".$left_index."-".$right_index . "-$difference" . "_tsd3-$tsd3";
                    print $element_out ">$title\n$element\n";
                    
                    $c++;
                    next LINE;
                }
                else {
                  if ((length($tsd3) > length($tsd2)) and (substr($tsd3, 0, 4) eq $tsd2) and (substr($tsd3, -4) eq $tsd2)) {
                    my $left_tir = $matches{'hit'}->[1];
                    my $right_tir = $matches{'query'}->[1];
                    $right_tir =~ tr/ATGCatgc/TACGtacg/;
                    $right_tir = reverse($right_tir);
                    my $difference = ($half - $left_index) + (($seq_len - $right_index) - $half);
                    my $element = substr($seq, $left_index, -($right_index));
                    
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
                    
                    print $tir_out ">$seq_id_ori\_ltir$c\_$left_index\-$right_index\-$difference\_tsd3\-$tsd3\n$left_tir\n";
                    print $tir_out ">$seq_id_ori\_rtir$c\_$left_index\-$right_index\-$difference\_tsd3\-$tsd3\n$right_tir\n";
                    
                    print $left_tir_out ">$seq_id_ori\_ltir$c\_$left_index\-$right_index\-$difference\_tsd3\-$tsd3\n$left_tir\n";
                    print $right_tir_out ">$seq_id_ori\_rtir$c\_$left_index\-$right_index\-$difference\_tsd3\-$tsd3\n$right_tir\n";
                    
                    my $title = $seq_id_ori."_tir".$c."_".$left_index."-".$right_index . "-$difference" . "_tsd3-$tsd3";
                    print $element_out ">$title\n$element\n";
                    $c++;
                    next LINE;
                  }
                  else {
                    my $left_tir = substr($matches{'hit'}->[1], 4);
                    my $right_tir = $matches{'query'}->[1];
                    $right_tir =~ tr/ATGCatgc/TACGtacg/;
                    $right_tir = reverse($right_tir);
                    $right_tir = substr($right_tir, 0, -4);
                    $left_index = $left_index + 4;
                    $right_index = $right_index - 4;
                    my $difference = ($half - $left_index) + (($seq_len - $right_index) - $half);
                    my $element = substr($seq, $left_index, -($right_index));
                    
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
                    
                    print $tir_out ">$seq_id_ori\_ltir$c\_$left_index\-$right_index\-$difference\_tsd2\-$tsd2\n$left_tir\n";
                    print $tir_out ">$seq_id_ori\_rtir$c\_$left_index\-$right_index\-$difference\_tsd2\-$tsd2\n$right_tir\n";
                    
                    print $left_tir_out ">$seq_id_ori\_ltir$c\_$left_index\-$right_index\-$difference\_tsd2\-$tsd2\n$left_tir\n";
                    print $right_tir_out ">$seq_id_ori\_rtir$c\_$left_index\-$right_index\-$difference\_tsd2\-$tsd2\n$right_tir\n";
                    
                    my $title = $seq_id_ori."_tir".$c."_".$left_index."-".$right_index . "-$difference" . "_tsd2-$tsd2";
                    print $element_out ">$title\n$element\n";                  
                    
                    $c++;
                    next LINE;
                  }
                }
            }
            elsif ($tsd2 and $tsd4_catch == 0) {
                my $left_tir = substr($matches{'hit'}->[1], 4);
                my $right_tir = $matches{'query'}->[1];
                $right_tir =~ tr/ATGCatgc/TACGtacg/;
                $right_tir = reverse($right_tir);
                $right_tir = substr($right_tir, 0, -4);
                $left_index = $left_index + 4;
                $right_index = $right_index - 4;
                my $difference = ($half - $left_index) + (($seq_len - $right_index) - $half);
                my $element = substr($seq, $left_index, -($right_index));
                
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
                
                print $tir_out ">$seq_id_ori\_ltir$c\_$left_index\-$right_index\-$difference\_tsd2\-$tsd2\n$left_tir\n";
                print $tir_out ">$seq_id_ori\_rtir$c\_$left_index\-$right_index\-$difference\_tsd2\-$tsd2\n$right_tir\n";
                
                print $left_tir_out ">$seq_id_ori\_ltir$c\_$left_index\-$right_index\-$difference\_tsd2\-$tsd2\n$left_tir\n";
                print $right_tir_out ">$seq_id_ori\_rtir$c\_$left_index\-$right_index\-$difference\_tsd2\-$tsd2\n$right_tir\n";
                
                my $title = $seq_id_ori."_tir".$c."_".$left_index."-".$right_index . "-$difference" . "_tsd2-$tsd2";
                print $element_out ">$title\n$element\n";
                $c++;
                next LINE;
            }
            
            

            my $left_tir = $matches{'query'}->[1];
            my $right_tir = $matches{'hit'}->[1];
            $right_tir =~ tr/ATGCatgc/TACGtacg/;
            $right_tir = reverse($right_tir);
            my $element = substr($seq, $left_index, -($right_index));
            
            my $diff = ($half - $left_index) + (($seq_len - $right_index) - $half);
            
            my $title = $seq_id_ori . "_tir" . $c . "_" . $diff;
            my $ltitle = $seq_id_ori . "_ltir" . $c . "_" . $diff;
            my $rtitle = $seq_id_ori . "_rtir" . $c . "_" . $diff;
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
            
            print $left_tir_out ">$ltitle\n$left_tir\n";
            print $right_tir_out ">$rtitle\n$right_tir\n";
            
            $c++;
        }
        else {
            print $no_tir_out "$seq_id_ori\n";
        }    
    }
}
close($tir_out);
close($left_tir_out);
close($right_tir_out);
close($no_tir_out);
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
        my $match_cutoff = 10;
        my $first_match = 0;

        #parse homology string, keeping track of match length and mismatches or gaps
        for (my $count = 0 ; $count < length($homo_string) ; $count++) {
          my $homo_char  = substr($homo_string, $count, 1);
          my $query_char = substr($query_str,   $count, 1);
          my $hit_char   = substr($hit_str,     $count, 1);
          if ($round == 1 and $count == 8 and $total_mis_aln >= 4) {
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
          if ($round == 3 or $round == 1) {
            if (($count <= 3) and $total_mis_aln >= 1) {
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
                $total_mis_aln++;
                next;
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
          elsif ($match_len >= 1 and $match_len < $match_cutoff) {

            #if match length is 1-3 and position is not a match, increment mismatch counter and check if more than one mismatch has occurred
            if ($homo_char eq " ") {
              $match_mis_aln++;
              $total_mis_aln++;

              #allow one mismatch, store info and continue
              if ($match_mis_aln <= 2) {
                $match_len++;
                $last_good = $count;
                $match_query .= $query_char;
                $match_hit   .= $hit_char;

                # print "First Mismatch at $count\n";
                next;
              }

              else {
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
          elsif ($match_len >= $match_cutoff) {

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


