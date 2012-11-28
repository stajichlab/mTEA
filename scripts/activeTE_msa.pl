#!/usr/bin/perl

use strict;
use warnings;

use Bio::AlignIO;
use Bio::SearchIO;
use Data::Dumper;
use File::Basename;
use Getopt::Long;

#set default flank length and reset if option is given
my $flank = 100;
GetOptions ('flank:i' => \$flank);

#make sure only one input file is specified
my $len = @ARGV;
if ($len != 1) {
    die "Please specify one input file\n";
}

#this was to check that default or supplied flank length was set correctly
#print "$flank\n";

#open input file, create input object, and breakdown directory and filename
my $infile = shift @ARGV;
my $in_obj = Bio::AlignIO->new(-file => $infile, -format => 'fasta');
#print Dumper($in_obj);
my($filename, $out_dir, $suffix) = fileparse($infile);

#initialize a global array to store results later
my @good_aln = ();
my $initial_left_TSD;
my $initial_right_TSD;

#create an alignment object from input object
my $full_aln_obj = $in_obj->next_aln();

#grab all sequences from alignment
foreach my $seq_obj ( $full_aln_obj->each_seq() ) {
    
    #grab sequence name and shorten it for use in some output filenames
    my $seq_name = $seq_obj->id();
    my @seqn_part = split(":", $seq_name);
    my $fname_start = $seqn_part[0];
    
    #grab sequence, make a copy, strip leading and trailing hyphens
    my $seq_ori = $seq_obj->seq();
    my $seq = $seq_ori;    
    $seq =~ s/^-*//;
    $seq =~ s/-*$//;
    
    #initialize variables
    my $sub_seq;
    my $first;
    my $last;
    my $first_path = $out_dir . "first_half.txt";
    my $last_path = $out_dir . "last_half.txt";
    
    #strip dashes, cut off most of the flanks, and then get end  sequences
    $sub_seq = $seq;
    $sub_seq =~ s/-//g;
    $first = substr($sub_seq, $flank-5, 25);
    $last = substr($sub_seq, length($sub_seq)-($flank+19), 25);
    #print "first: $first\nlast: $last\n";
    
    #save the two ends as files to use as inputs for a FASTA search
    open(my $first_out, ">", $first_path) or die "Error creating $first_path. $!\n";
    open(my $last_out, ">", $last_path) or die "Error creating $last_path. $!\n";
    print $first_out ">first\n$first\n";
    print $last_out ">last\n$last\n";
    close($first_out);
    close($last_out);
    #print "$first\t$last\n";
    #print "$seq_name\n";
    
    #create fasta output filename then call FASTA
    my $out_opt = $out_dir . $fname_start . ".fasta.out";
    system("fasta36 -n -i -T 8 -d 1 $last_path $first_path > $out_opt");
    
    #open FASTA output file and creat input object
    my $fa_aln_obj = Bio::SearchIO->new(-format => 'fasta', -file   => $out_opt);
    
    #go through the FASTA input object . . down to the HSP
    while( my $result = $fa_aln_obj->next_result ) {
        while( my $hit = $result->next_hit ) {
            while( my $hsp = $hit->next_hsp ) {
                #grab the query, hit, and homology strings
                my $homo_string = $hsp->homology_string;
                my $query_str = $hsp->query_string;
                my $hit_str = $hsp->hit_string;
                $seq =~ s/-//g;
                #$seq =~ s/ /-/g;
                $hit_str =~ s/ /-/g;
                
                #initialize variables
                my $match_len = 0;
                my $start_pos = '';
                my $end_pos = '';
                my $match_query = '';
                my $match_hit = '';
                my $mis_aln = 0;
                my $last_good = 0;
                my $hit_pos;
                my $query_pos;
                
                
                #parse homology string, keeping track of match length and mismatches or gaps
                
                for (my $count = 0; $count < length($homo_string); $count++) {
                    my $homo_char =  substr($homo_string, $count, 1);
                    my $query_char =  substr($query_str, $count, 1);
                    my $hit_char =  substr($hit_str, $count, 1);
                    #print "$homo_char\n";
                    
                    if ($match_len == 0){
                        #if match length equals 0 and position is not a match, continue to next position
                        if ($homo_char eq " ") {
                            next;
                        }
                        #if position is a match, store info, continue to next position
                        elsif ($homo_char eq ":") {
                            $start_pos = $count;
                            $last_good = $count;
                            #print "1st If/1st elsif\t$mis_aln $match_len \n";
                            $match_len++;
                            $match_query .= $query_char;
                            $match_hit .= $hit_char;
                            next;
                        }
                    }
                           
                    elsif ($match_len >= 1 and $match_len < 4) {
                        #if match length is 1-3 and position is not a match, increment mismatch counter and check if more than one mismatch has occurred
                        if ($homo_char eq " ") {
                            $mis_aln++;
                            #allow one mismatch, store info and continue
                            if ($mis_aln <= 1) {
                                $match_len++;
                                $last_good = $count;
                                $match_query .= $query_char;
                                $match_hit .= $hit_char;
                                #print "1st elsif/1st If\t$mis_aln $match_len \n";
                                next;
                            }
                            #more than one mismatch, reset counters and other info, continue
                            elsif ($mis_aln > 1){ 
                                $match_len = 0;
                                $start_pos = '';
                                $match_query = '';
                                $match_hit = '';
                                $mis_aln = 0;
                                #print "1st elsif/1st elsif\t$mis_aln $match_len \n";
                                next;
                            }
                        }
                        #position is a match, store info and continue
                        elsif ($homo_char eq ":") {
                            #print "1st elsif/2nd elsif\t$mis_aln $match_len \n";
                            $last_good = $count;
                            $match_len++;
                            $match_query .= $query_char;
                            $match_hit .= $hit_char;
                            next;
                        }
                    }
                    
                    elsif ($match_len >= 4) {
                        #match length is 4 or higher. If position is not a match, increment mismatch counter and check if more than 2 mismatch has occurred. If a match, continue.
                        if ($homo_char eq " ") {
                            $mis_aln++;
                            #mismatches under 3, store info and continue
                            if ($mis_aln < 3) {
                                $match_len++;
                                $last_good = $count;
                                $match_query .= $query_char;
                                $match_hit .= $hit_char;
                                #print "Last-1st If\t$mis_aln $match_len \n";
                                next;
                            }
                            #mismatches 3 or more, store final info for alignment match and end parsing
                            elsif($mis_aln >= 3){
                                $end_pos = $last_good;
                                $match_query =~ s/-//g;
                                $match_hit =~ s/-//g;
                                
                                #reverse complement the match query sequence
                                $match_query =~ tr/ATGC/TACG/;
                                $match_query = reverse($match_query);
                                
                                my $match_query_len = length($match_query);
                                my $match_hit_len = length($match_hit);
                                
                                #find the position in the full sequence of the hit and query match sequences
                                $hit_pos = index($seq, $match_hit, 50)+1;
                                $query_pos = index($seq, $match_query, $hit_pos)+1;
                                #store sequence name and the hit and query info                                
                                my @match = ($seq_name, {"hit" =>[$hit_pos, $match_hit_len],"query" => [$query_pos, $match_query_len]});
                                push @good_aln, [@match];
                                last;
                            }
                        }
                        #position is a match, store info and continue
                        elsif ($homo_char eq ":") {
                            $last_good = $count;
                            $match_len++;
                            $match_query .= $query_char;
                            $match_hit .= $hit_char;
                            #print "Last-Very last elsif\t$mis_aln $match_len \n";
                            next;
                        }
                    }
                }
                #add in check to see if TIRs were found. if not, remove the sequence from the MSA, and print out to 'no TIRs list'
                if ($end_pos eq '') {
                    #print out to no TIRs list
                    $full_aln_obj->remove_seq($seq_obj);
                    next;
                }
            }
        }
    }
}

#find the TIRs in the staring alignment using the index positions of each TIR stored above

#initialize variables
my %hit_column_counts;
my %hit_match_len;
my %query_column_counts;
my %query_match_len;
my $query_aln_pos;
my $hit_aln_pos;
my @entry;

#go through each entry of @good_aln array
foreach my $row_ref (@good_aln) {
    @entry = @{$row_ref};
    #get column position in alignment of left (hit) TIR, using sequence name and index position of TIR in full sequence. Increment the count of that column position in %hit_column_counts. Also, increment the count of TIR length in %hit_match_len
    $hit_aln_pos = $full_aln_obj->column_from_residue_number( $entry[0], ${entry[1]{"hit"}}[0]);
    $hit_column_counts{$hit_aln_pos}++;
    $hit_match_len{${entry[1]{"hit"}}[1]}++;
    #do the same for the query (left) TIR
    $query_aln_pos = $full_aln_obj->column_from_residue_number( $entry[0], ${entry[1]{"query"}}[0]);
    $query_column_counts{$query_aln_pos}++;
    $query_match_len{${entry[1]{"query"}}[1]}++;
    #print "$entry[0]\nleft tsd: $entry[1]{'hit'}[2]\nright tsd: $entry[1]{'query'}[2]\n";
    #print "test: $query_match_len{${entry[1]{'query'}}[1]}\n";
}
#sort the hit and query hashes by largest count to smallest
my @sorted_hitcolumn_keys =  sort { $hit_column_counts{$b} <=> $hit_column_counts{$a} } keys (%hit_column_counts);
my @sorted_querycolumn_keys =  sort { $query_column_counts{$b} <=> $query_column_counts{$a} } keys (%query_column_counts);
my @sorted_hit_len_keys = sort { $hit_match_len{$b} <=> $hit_match_len{$a} } keys (%hit_match_len);
my @sorted_query_len_keys = sort { $query_match_len{$b} <=> $query_match_len{$a} } keys (%query_match_len);

#Extract the left and right TIRs as new alignment objects, starting at the most common starting column until that column plus the most common TIR length
my $left_TIR_aln_obj = $full_aln_obj->slice($sorted_hitcolumn_keys[0], $sorted_hitcolumn_keys[0]+$sorted_hit_len_keys[0]-1);
my $right_TIR_aln_obj = $full_aln_obj->slice($sorted_querycolumn_keys[0], $sorted_querycolumn_keys[0]+$sorted_query_len_keys[0]-1);

#calculate the overall percent identity of each TIR
my $right_tir_id1 = $right_TIR_aln_obj->percentage_identity;
my $left_tir_id1 = $left_TIR_aln_obj->percentage_identity;
#print "left start: $left_tir_id1\n";
#print "right start: $right_tir_id1\n";

#calculate the percent identity for indiviidual positions of both TIRs
my $left_pos = $left_TIR_aln_obj->length();
my $right_pos = $right_TIR_aln_obj->length();

for (my $i = 1;  $i <= $left_pos; $i++ ) {
    my $left_TIR_pos = $left_TIR_aln_obj->slice($i,$i);
    my $left_pos_id = $left_TIR_pos->percentage_identity;
    #print "left_pos: $left_pos_id\n";
}
for (my $ii = 1;  $ii <= $right_pos; $ii++ ) {
    my $right_TIR_pos = $right_TIR_aln_obj->slice($ii,$ii);
    my $right_pos_id = $right_TIR_pos->percentage_identity;
}



#find TSDs based off column # of TIRs. Must check for 2-bp TSD included in TIRS and then a longer TSD. Longer TSDs that start and finish with same 2-bp as the 2-bp TSD supercede the 2-bp TSD. Also, grab flanks and save them to files to be used by CD-Hit to cluster them.

my @putative_TSD;
my @TSD_info;

my $flanks_out_path = $out_dir . "flanks.fa";

open(my $flanks_out, ">", $flanks_out_path) or die "Error creating $flanks_out_path. $!\n";

foreach my $seq_obj ( $full_aln_obj->each_seq() ) {
    
    #grab sequence name and shorten it for use in some output filenames
    my $seq_name = $seq_obj->id();
    my @seqn_part = split(":", $seq_name);
    my $fname_start = $seqn_part[0];
    #print "$seq_name\n";
    
    #grab sequence, make a copy, strip hyphens
    my $seq_ori = $seq_obj->seq();
    my $seq = $seq_ori;
    $seq =~ s/-//g;
    
    #get seq_pos of the start of TIRs from column # in alignment. Adjust to include the 2-bp at outside ends of TIRs and then split into the 2-bp and 10-bp sequences to search for TSDs
    my $left_tsd_loc_obj = $seq_obj->location_from_column($sorted_hitcolumn_keys[0]);
    my $left_tsd_end_pos = $left_tsd_loc_obj->start();
    my $starting_left_flank = substr($seq, 0, $left_tsd_end_pos+1);
    my $left_tsd = substr($seq, $left_tsd_end_pos-21, 22);
    my $left_tsd_substr = substr($left_tsd, 10, 10);
    #print "left tsd end: $left_tsd_end_pos\t$left_tsd\n";
    
    my $right_tsd_loc_obj = $seq_obj->location_from_column(($sorted_querycolumn_keys[0]+$sorted_query_len_keys[0]-1));
    my $right_tsd_start_pos = $right_tsd_loc_obj->start();
    my $starting_right_flank = substr($seq, $right_tsd_start_pos-2);
    my $right_tsd = substr($seq, $right_tsd_start_pos-2, 22);
    my $right_tsd_substr = substr($right_tsd, 2, 10);
    #print "right tsd start: $right_tsd_start_pos\t$right_tsd\n";
    
    my $left_2bp = substr($left_tsd_substr,-2);
    my $right_2bp = substr($right_tsd_substr, 0, 2);
    #print "left 2bp: $left_2bp\tright 2bp: $right_2bp\n\n";
    
    my $tsd1;
    my $tsd2;
    for (my $i = 2; $i <= 12; $i++ ){
        if (substr($left_tsd,-($i)) eq substr($right_tsd, 0, $i)) {
            #print "Round 1 TSD search - Yes! $i\n";
            $tsd1 = substr($left_tsd,-($i));
            next;
        }
        else{
            #print "Round 1 TSD search - No, $i. Killing\n\n";
            last;
        }
    }
    for (my $i = 0; $i < 10; $i++ ){
        if (substr($left_tsd_substr, $i) eq substr($right_tsd_substr, 0, -($i))) {
            #print "Round 2 TSD search - Yes! $i\n";
            $tsd2 = substr($left_tsd_substr,-($i));
            next;
        }
        else{
            my $sub_print = substr($left_tsd_substr,$i);
            #print "$sub_print\n";
            #print "Round 2 TSD search - No, $i. Next . . .\n\n";
            next;
        }
    }   
    #Save found TSD to an array or report that a TSD wasn't found
    if ($tsd1) {
        #print "First tsd1 test\n";
        if (! $tsd2) {
            my $insertion_site = substr($left_tsd, -12) . substr($right_tsd, 2, 10);
            push @TSD_info, ($seq_name, $insertion_site);            
            push @putative_TSD, $tsd1;
            #print "For insertion site preference:  $insertion_site\n";
            my $left_flank = substr($starting_left_flank, 0, -2);
            my $right_flank = substr($starting_right_flank, 2);
            my $flanks = $left_flank . $right_flank;
            print $flanks_out ">$fname_start\n$flanks\n";
        }
        else {
            #if both a 2-bp TSD and a longer one is found, check if the longer one has the same 2-bp (as the found 2-bp TSD) at both ends. If so, use the longer TSD. Otherwise use the 2-bp TSD
            if ((length($tsd2) > length($tsd1)) and (substr($tsd2, 0, 2) eq $tsd1) and (substr($tsd2, -2) eq $tsd1))  {
                my $insertion_site = substr($left_tsd, (-2-length($tsd2)-10), (length($tsd2)+10)) . substr($right_tsd, 2+length($tsd2), 10);
                push @TSD_info, ($seq_name, $insertion_site);
                push @putative_TSD, $tsd2;
                my $left_flank = substr($starting_left_flank, 0, (-2-length($tsd2)));
                my $right_flank = substr($starting_right_flank, (2+length($tsd2)));
                my $flanks = $left_flank . $right_flank;
                print $flanks_out ">$fname_start\n$flanks\n";
            }
            else {
                my $insertion_site = substr($left_tsd, -12) . substr($right_tsd, 2, 10);
                push @TSD_info, ($seq_name, $insertion_site);
                push @putative_TSD, $tsd1;
                my $left_flank = substr($starting_left_flank, 0, -2);
                my $right_flank = substr($starting_right_flank, 2);
                my $flanks = $left_flank . $right_flank;
                print $flanks_out ">$fname_start\n$flanks\n";
            }
        }
    }
    else {
        if ($tsd2) {
            my $insertion_site = substr($left_tsd, (-2-length($tsd2)-10), (length($tsd2)+10)) . substr($right_tsd, 2+length($tsd2), 10);
            push @TSD_info, ($seq_name, $insertion_site);
            push @putative_TSD, [$tsd2];
            my $left_flank = substr($starting_left_flank, 0, (-2-length($tsd2)));
            my $right_flank = substr($starting_right_flank, (2+length($tsd2)));
            my $flanks = $left_flank . $right_flank;
            print $flanks_out ">$fname_start\n$flanks\n";
        }
        else {
            #store or print to report that TSD was not found for this copy
            # use $seq_name, etc
        }
    }
}

close($flanks_out);


#Determine the most common TSD by sequence first if possible or by length. If >80% of the TSDs are the same, then that sequence is stored as the TSD for output. Otherwise, look at the lengths of the TSDs and store the length is >80% of the TSDs have the same length.
my %TSD_counts;

my $TSD_array_length = @putative_TSD;
print "tsd array length: $TSD_array_length\n";

#count the occurances of a TSD sequence and sort with the highest value first in an array
foreach my $row (@putative_TSD) {
    #@entry = @{$row_ref};
    $TSD_counts{$row}++;
}
print Dumper(\%TSD_counts);
my @sorted_TSD_keys =  sort { $TSD_counts{$b} <=> $TSD_counts{$a} } keys (%TSD_counts);

#check whether the same TSD sequence was found in >80% of the copies. If not, count the occurances of TSD length, sort by highest occurances, and check if the length of the TSD is the same in >80% of the copies
my $final_TSD_length;
my $final_TSD_seq;
if (($TSD_counts{$sorted_TSD_keys[0]} / $TSD_array_length) > 0.8) {
    print "Check works!\n";
    $final_TSD_length = length($sorted_TSD_keys[0]);
    $final_TSD_seq = $sorted_TSD_keys[0];
}
else{
    my %TSD_length_counts;
    foreach my $row (@putative_TSD) {
        $TSD_length_counts{length($row)}++;
    }
    my @sorted_TSD_length_keys =  sort { $TSD_length_counts{$b} <=> $TSD_length_counts{$a} } keys (%TSD_length_counts);
    $final_TSD_length = $sorted_TSD_length_keys[0];
    #print Dumper(\%TSD_length_counts);
}
#rint "final TSD length: $final_TSD_length\n";

#get flanks to cluster by CD-Hit and also to get 10-bp on each side of the original TSD to look at target site preferences


my $flanks_cluster_path = $flanks_out_path . ".cluster";


system("cd-hit-est -i $flanks_out_path -o $flanks_cluster_path -c 0.8 -g 1 -T 0 -n 5 -d 0");
my $muscle_out_path = $flanks_out_path . ".msa";
system("muscle -in $flanks_out_path -out $muscle_out_path");
    



__END__

#grab the part of the sequence surrounded by dashes (from MSA) if possible and get ends
    if ($seq =~ m/-(\w{5,}-*\w{105,}[\w-]+?)-/) {
        $sub_seq = $1;
        $sub_seq =~ s/-//g;
        $first = substr($sub_seq, 0, 33);
        $last = substr($sub_seq, -33);
    }
    #or, grab the left part of the sequence with dashes (from MSA) and get right side by cutting off most of flank, get ends
    elsif ($seq =~ m/-(\w{5,}-*w{105,}[\w-]+?)/) {
        $sub_seq = $1;
        $sub_seq =~ s/-//g;
        $first = substr($sub_seq, 0, 33);
        $last = substr($seq, length($seq)-($flank+20), 33);
    }
    #or, grab the right part of the sequence with dashes (from MSA) and get left side by cutting off most of flank, get ends
    elsif($seq =~ m/(\w{5,}-*w{105,}[\w-]+?)-/) {
        $sub_seq = $1;
        $sub_seq =~ s/-//g;
        $first = substr($seq, $flank-5, 33);
        $last = substr($sub_seq, -33);
    }
