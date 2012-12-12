#!/usr/bin/perl

use strict;
use warnings;

use Bio::AlignIO;
use Bio::SearchIO;
use Data::Dumper;
use File::Spec;
use Getopt::Long;

#set default flank length and reset if option is given
#no longer needed but will keep here for now
#other options will likely be added
#my $flank = 100;
#GetOptions ('flank:i' => \$flank);

#make sure only one input file is specified
my $len = @ARGV;
if ($len != 1) {
    die "Please specify one input file\n";
}
my $infile = shift @ARGV;
if (!-f $infile) {
        die "Supplied filepath is not valid";
    }
#this was to check that default or supplied flank length was set correctly
#print "$flank\n";

#breakdown directory and filename, create output directory
my ($volume, $in_dir, $filename) = File::Spec->splitpath($infile);
my $out_dir = "activeTE_out_" . $filename;
my $out_path = File::Spec->catdir(($in_dir, $out_dir));
if (!-d $out_path) {
    system("mkdir $out_path");
}

#my $out_file_path = File::Spec->catpath($volume, $out_path, Put File name here!);
#generate trimal alignment gap summary file, using sed to remove the header info
my $gap_summary_out = File::Spec->catpath($volume, $out_path, $filename . ".gap_summary");
system("trimal -in $infile -sgc | sed -n '4,\$p' > $gap_summary_out");

#initialize a global array to store results later
my @good_aln = ();
my $initial_left_TSD;
my $initial_right_TSD;

#create input object and an alignment object from it
my $in_obj = Bio::AlignIO->new(-file => $infile, -format => 'fasta');
my $full_aln_obj = $in_obj->next_aln();
my $full_aln_len = $full_aln_obj->length();
my @full_id_array;
#calculate the % identity of each position in the full alignment and print to file
for (my $i = 1;  $i <= $full_aln_len; $i++ ) {
    my $pos = $full_aln_obj->slice($i,$i,1);
    my $pos_id = $pos->percentage_identity;
    my @info = ($i, $pos_id);
    push @full_id_array, [@info];
}
my $full_id_out_path = File::Spec->catpath($volume, $out_path, $filename . ".full_id");
open(my $id_out, ">", $full_id_out_path);

foreach my $row_ref (@full_id_array) {
    my @pos = @{$row_ref};
    print $id_out "$pos[0]\t$pos[1]\n";
}
close($full_id_out_path);

#generate a gap column matric fro the full alignment
my $gap_cols = $full_aln_obj->gap_col_matrix();

#open and parse the trimal output file with the % gap ID of each position in the full alignment
open(my $gap_id, "<", $gap_summary_out);
my @gap_id_array;
while (my $line = <$gap_id>){
    chomp $line;
    #print "$line\n";
    my @info = $line =~ /\S+\s+(\S+)\s+\S+/;
    #my @info = split(/\s+/, $line);
    push @gap_id_array, @info
}

#initialize variable to be used in determining the sequences to be removed from the full alignment based on high % gap ID and low sequence % ID
my %gap_seq_remove;
my %gap_seq_count;
my @gap_seq_pos_remove;

for (my $i = 0; $i < $full_aln_len; $i++){
    #print "$id_array[$i]\n";
    my $id_row_ref = $full_id_array[$i];
    my @id_info = @{$id_row_ref};
    my @gap_col_array = @{$gap_cols};
    my $gap_col_hashref = $gap_col_array[$i];
    my %gap_col_hash = %{$gap_col_hashref};
    if ($gap_id_array[$i] >= 97.0 and $id_info[1] <= 34.0){
        foreach my $key (keys %gap_col_hash) {
            #print "KEY: $key\n";
            if ($gap_col_hash{$key} != 1){
                #print "$key\n\n";
                my $seq_obj = $full_aln_obj->get_seq_by_id($key);
                #print "$seq_obj\n";
                my $seq = $seq_obj->seq();
                #$seq =~ s/'-'//g;
                my $seq_pos = $seq_obj->location_from_column($i+1);
                my @info = ($key, $i+1, $seq_obj);
                push @gap_seq_pos_remove, [@info];
                $gap_seq_count{$key}++;
                $gap_seq_remove{$key} = $seq_obj;
            }
        }
    }
}
foreach my $key (keys %gap_seq_remove) {
    my $seq_obj = $gap_seq_remove{$key};
    #print Dumper(\$seq_obj);
    $full_aln_obj->remove_seq($seq_obj);
}

$full_aln_obj = $full_aln_obj->remove_gaps('-',1);
my $trim_aln_out = File::Spec->catpath($volume, $out_path, $filename . ".trim");
my $out = Bio::AlignIO->new(-file => ">$trim_aln_out", -format => 'fasta');
$out->write_aln($full_aln_obj);

my @trim_id_array;
my $trim_aln_len = $full_aln_obj->length();
for (my $i = 1;  $i <= $trim_aln_len; $i++ ) {
    my $pos = $full_aln_obj->slice($i,$i,1);
    my $pos_id = $pos->percentage_identity;
    my @info = ($i, $pos_id);
    push @trim_id_array, [@info];
}
my $trim_id_out_path = File::Spec->catpath($volume, $out_path, $filename . ".trim_id");
open(my $trim_id_out, ">", $trim_id_out_path);

foreach my $row_ref (@trim_id_array) {
    my @pos = @{$row_ref};
    print $trim_id_out "$pos[0]\t$pos[1]\n";
}
close($trim_id_out_path);

my $left_tir_start;
my $left_tir_count = 0;
for (my $i = 1;  $i < $trim_aln_len; $i++ ) {
    my $pos = $full_aln_obj->slice($i,$i,1);
    my $pos_id = $pos->percentage_identity;
    #print "left %ID: $pos_id\t Step: $i\tCount:$left_tir_count\n";
    if ($left_tir_count == 0) {
        if ($pos_id >= 90.0) {
            $left_tir_start = $i;
            $left_tir_count++;
            next;
        }
        else {
            next;
        }
    }
    elsif ($left_tir_count >= 1 and $left_tir_count <= 4) {
        if ($pos_id >= 90.0) {
            $left_tir_count++;
            next;
        }
        else {
            $left_tir_count = 0;
            undef($left_tir_start);
            next;
        }    
    }
    elsif ($left_tir_count >=5) {
        last;
    }
}

my $right_tir_start;
my $right_tir_count = 0;
for (my $i = 0;  $i < $trim_aln_len; $i++ ) {
    my $pos = $full_aln_obj->slice($trim_aln_len-$i,$trim_aln_len-$i,1);
    my $pos_id = $pos->percentage_identity;
    if ($right_tir_count == 0) {
        if ($pos_id >= 90.0) {
            $right_tir_start = $trim_aln_len-$i;
            $right_tir_count++;
            next;
        }
        else {
            next;
        }
    }
    elsif ($right_tir_count >= 1 and $right_tir_count <= 4) {
        if ($pos_id >= 90.0) {
            $right_tir_count++;
            next;
        }
        else {
            $right_tir_count = 0;
            undef($right_tir_start);
            next;
        }    
    }
    elsif ($right_tir_count >= 5) {
        last;
    }
}

#need to store the column positions of the TIRs in the trimmed alignment and get the coulumns of them in the original alignminent to remove only sequences that cause gaps in the tirs
my %trim_pos_hash;
my %trim_left_pos_hash;
my %trim_right_pos_hash;

foreach my $seq_obj ($full_aln_obj->each_seq()) {
    my $left_res_pos_obj = $seq_obj->location_from_column($left_tir_start);
    my $right_res_pos_obj = $seq_obj->location_from_column($right_tir_start);
    my $seq = $seq_obj->seq();
    my $left_res_pos = $left_res_pos_obj->start();
    my $left_res = substr($seq, $left_res_pos, 1);
    my $right_res_pos = $right_res_pos_obj->start();
    my $right_res = substr($seq, $right_res_pos, 1);
    my @trim_pos_array = ($left_res_pos, $right_res_pos);
    my $seq_name = $seq_obj->id();
    $trim_pos_hash{$seq_name} = @trim_pos_array;
    $trim_left_pos_hash{$seq_name} = $left_res_pos;
    $trim_right_pos_hash{$seq_name} = $right_res_pos;
}

#grab all sequences from trimmed alignment
foreach my $seq_obj ( $full_aln_obj->each_seq() ) {
    #grab sequence name and shorten it for use in some output filenames
    my $seq_name = $seq_obj->id();
    my @seqn_part = split(":", $seq_name);
    my $fname_start = $seqn_part[0];
    
    #grab sequence, make a copy, strip leading and trailing hyphens
    my $seq = $seq_obj->seq();
    $seq =~ s/-//g;
    
    #initialize variables
    my $sub_seq;
    my $first;
    my $last;
    my $first_path = File::Spec->catpath($volume, $out_path, $filename . "first_half.txt");
    my $last_path = File::Spec->catpath($volume, $out_path, $filename . "last_half.txt");
    
    #get end  sequences
    $first = substr($seq, $trim_left_pos_hash{$seq_name}-1, 35);
    $last = substr($seq, $trim_right_pos_hash{$seq_name}-35,35);
    
    #save the two ends as files to use as inputs for a FASTA search
    open(my $first_out, ">", $first_path) or die "Error creating $first_path. $!\n";
    open(my $last_out, ">", $last_path) or die "Error creating $last_path. $!\n";
    print $first_out ">first\n$first\n";
    print $last_out ">last\n$last\n";
    close($first_out);
    close($last_out);
    
    #create fasta output filename then call FASTA
    my $out_opt = File::Spec->catpath($volume, $out_path, $fname_start . ".ggsearch.out");
    system("ggsearch36 -n -i -T 8 -d 1 $last_path $first_path > $out_opt");
    
    my @tir_match_result = match_tirs($seq_obj, $out_opt);
    
    if ($tir_match_result[0] == 1) {
        push @good_aln, $tir_match_result[1];
    }
    else {
        $full_aln_obj->remove_seq($seq_obj);
    }
}

#find the TIRs in the trimmed alignment using the index positions of each TIR stored above
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
}

#sort the hit and query column and match length hashes by largest count to smallest
my @sorted_hitcolumn_keys =  sort { $hit_column_counts{$b} <=> $hit_column_counts{$a} } keys (%hit_column_counts);
my @sorted_querycolumn_keys =  sort { $query_column_counts{$b} <=> $query_column_counts{$a} } keys (%query_column_counts);
my @sorted_hit_len_keys = sort { $hit_match_len{$b} <=> $hit_match_len{$a} } keys (%hit_match_len);
my @sorted_query_len_keys = sort { $query_match_len{$b} <=> $query_match_len{$a} } keys (%query_match_len);

#get the end residue positions for the tirs from one sequence in the trimmed alignment that will be used to find the tir interval in the full alignment. Then sequences that were excluded can be added back if they don't cause gaps in the tirs
my $trim_seq_obj = $full_aln_obj->get_seq_by_pos(1);
my $trim_seq_name = $trim_seq_obj->id();
my $trim_seq = $trim_seq_obj->seq();
my $trim_seq_left_tir_end_pos =  $trim_left_pos_hash{$trim_seq_name}+$sorted_hit_len_keys[0]-1;
my $trim_seq_right_tir_end_pos_obj = $trim_seq_obj->location_from_column($sorted_querycolumn_keys[0]);
my $trim_seq_right_tir_end_pos = $trim_seq_right_tir_end_pos_obj->start();

#reread the original alignment file
my $in_obj2 = Bio::AlignIO->new(-file => $infile, -format => 'fasta');
my $ori_aln_obj = $in_obj2->next_aln();
my $ori_aln_len = $ori_aln_obj->length();

#get the column positions of tirs in the original alignment
my $left_tir_start_full = $ori_aln_obj->column_from_residue_number($trim_seq_name, $trim_left_pos_hash{$trim_seq_name});
my $left_tir_end_full = $ori_aln_obj->column_from_residue_number($trim_seq_name, $trim_seq_left_tir_end_pos);
my $right_tir_start_full = $ori_aln_obj->column_from_residue_number($trim_seq_name, $trim_right_pos_hash{$trim_seq_name});
my $right_tir_end_full = $ori_aln_obj->column_from_residue_number($trim_seq_name, $trim_seq_right_tir_end_pos);

#go back through the full alignment and remove only the sequences that create gaps in the tirs
my %gap_seq_remove2;
my @search_tirs;

foreach my $row_ref (@gap_seq_pos_remove) {
    @entry = @{$row_ref};
    if (($entry[1] >= $left_tir_start_full and $entry[1] <= $left_tir_end_full) or ($entry[1] >= $right_tir_end_full and $entry[1] <= $right_tir_start_full)){
        $gap_seq_remove2{$entry[0]}++;
    }
    else {
        push @search_tirs, $entry[0];
    }
}

foreach my $key (keys %gap_seq_remove2) {
    my $remove = $ori_aln_obj->get_seq_by_id($key);
    $ori_aln_obj->remove_seq($remove);
}
my $final_aln_obj = $ori_aln_obj->remove_gaps('-',1);

#get the column positions of tirs in the intermediate-final alignment
my $left_tir_start_final = $final_aln_obj->column_from_residue_number($trim_seq_name, $trim_left_pos_hash{$trim_seq_name});
my $left_tir_end_final = $final_aln_obj->column_from_residue_number($trim_seq_name, $trim_seq_left_tir_end_pos);
my $right_tir_start_final = $final_aln_obj->column_from_residue_number($trim_seq_name, $trim_right_pos_hash{$trim_seq_name});
my $right_tir_end_final = $final_aln_obj->column_from_residue_number($trim_seq_name, $trim_seq_right_tir_end_pos);

foreach my $seq_name (@search_tirs) {
    if (exists $gap_seq_remove2{$seq_name}) {
        next;
    }
    my $seq_obj = $final_aln_obj->get_seq_by_id($seq_name);
    my @seqn_part = split(":", $seq_name);
    my $fname_start = $seqn_part[0];
    my $seq = $seq_obj->seq();
    my $seq_len = length($seq);
    $seq =~ s/-//g;
    my $first;
    my $last;
    my $first_path = File::Spec->catpath($volume, $out_path, "first_half.txt");
    my $last_path = File::Spec->catpath($volume, $out_path, "last_half.txt");
    my $left_tir_start = $seq_obj->location_from_column($left_tir_start_final);
    my $left_tir_start_pos = $left_tir_start->start();
    my $right_tir_start = $seq_obj->location_from_column($right_tir_start_final);
    my $right_tir_start_pos = $right_tir_start->start();
        
    $first = substr($seq, $left_tir_start_pos-1, 35);
    $last = substr($seq, $right_tir_start_pos-35,35);
    
    open(my $first_out, ">", $first_path) or die "Error creating $first_path. $!\n";
    open(my $last_out, ">", $last_path) or die "Error creating $last_path. $!\n";
    print $first_out ">first\n$first\n";
    print $last_out ">last\n$last\n";
    close($first_out);
    close($last_out);
    my $out_opt = File::Spec->catpath($volume, $out_path, $fname_start . ".fasta.out");
    system("ggsearch36 -n -i -T 8 -d 1 $last_path $first_path > $out_opt");
    
    my @tir_match_result = match_tirs($seq_obj, $out_opt);
    if ($tir_match_result[0] == 0) {
        $full_aln_obj->remove_seq($seq_obj);
    }
}

$final_aln_obj = $ori_aln_obj->remove_gaps('-',1);
my $final_aln_out = File::Spec->catpath($volume, $out_path, $filename . ".final");
$out = Bio::AlignIO->new(-file => ">$final_aln_out", -format => 'fasta');
$out->write_aln($final_aln_obj);

#get the column positions of tirs in the final alignment
$left_tir_start_final = $final_aln_obj->column_from_residue_number($trim_seq_name, $trim_left_pos_hash{$trim_seq_name});
$left_tir_end_final = $final_aln_obj->column_from_residue_number($trim_seq_name, $trim_seq_left_tir_end_pos);
$right_tir_start_final = $final_aln_obj->column_from_residue_number($trim_seq_name, $trim_right_pos_hash{$trim_seq_name});
$right_tir_end_final = $final_aln_obj->column_from_residue_number($trim_seq_name, $trim_seq_right_tir_end_pos);


#Extract the left and right TIRs as new alignment objects, starting at the most common starting column until that column plus the most common TIR length
my $left_TIR_aln_obj = $final_aln_obj->slice($left_tir_start_final, $left_tir_end_final,1);
my $right_TIR_aln_obj = $final_aln_obj->slice($right_tir_end_final, $right_tir_start_final,1);

#calculate the overall percent identity of each TIR
my $right_tir_id1 = $right_TIR_aln_obj->percentage_identity;
my $left_tir_id1 = $left_TIR_aln_obj->percentage_identity;
#print "left start: $left_tir_id1\n";
#print "right start: $right_tir_id1\n";

#calculate the percent identity for indiviidual positions of both TIRs
my $left_pos = $left_TIR_aln_obj->length();
my $right_pos = $right_TIR_aln_obj->length();

for (my $i = 1;  $i <= $left_pos; $i++ ) {
    my $left_TIR_pos = $left_TIR_aln_obj->slice($i,$i,1);
    my $left_pos_id = $left_TIR_pos->percentage_identity;
    #print "left_pos: $left_pos_id\n";
}
for (my $ii = 1;  $ii <= $right_pos; $ii++ ) {
    my $right_TIR_pos = $right_TIR_aln_obj->slice($ii,$ii,1);
    my $right_pos_id = $right_TIR_pos->percentage_identity;
}

#find TSDs based off column # of TIRs. Must check for 2-bp TSD included in TIRS and then a longer TSD. Longer TSDs that start and finish with same 2-bp as the 2-bp TSD supersede the 2-bp TSD. Also, grab flanks and save them to files to be used by CD-Hit to cluster them.
my @putative_TSD;
my @TSD_info;
my $flanks_out_path = File::Spec->catpath($volume, $out_path, $filename . "flanks.fa");

open(my $flanks_out, ">", $flanks_out_path) or die "Error creating $flanks_out_path. $!\n";

foreach my $seq_obj ( $final_aln_obj->each_seq() ) {
    
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
    my $left_tsd_loc_obj = $seq_obj->location_from_column($left_tir_start_final);
    my $left_tsd_end_pos = $left_tsd_loc_obj->start();
    my $starting_left_flank = substr($seq, 0, $left_tsd_end_pos+1);
    my $left_tsd = substr($seq, $left_tsd_end_pos-21, 22);
    my $left_tsd_substr = substr($left_tsd, 10, 10);
    
    my $right_tsd_loc_obj = $seq_obj->location_from_column(($right_tir_start_final));
    my $right_tsd_start_pos = $right_tsd_loc_obj->start();
    my $starting_right_flank = substr($seq, $right_tsd_start_pos-2);
    my $right_tsd = substr($seq, $right_tsd_start_pos-2, 22);
    my $right_tsd_substr = substr($right_tsd, 2, 10);
    
    my $left_2bp = substr($left_tsd_substr,-2);
    my $right_2bp = substr($right_tsd_substr, 0, 2);
    
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
            $tsd2 = substr($left_tsd_substr,-($i));
            next;
        }
        else{
            my $sub_print = substr($left_tsd_substr,$i);
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
            #store or print to report that TSD was not found for this copy use $seq_name, etc
            #use the end position of tirs to grab the flanking sequences
            my $left_flank = substr($seq, 0, $left_tsd_end_pos-1);
            my $right_flank = substr($seq, $right_tsd_start_pos+1);
            my $flanks = $left_flank . $right_flank;
            print $flanks_out ">$fname_start\n$flanks\n";
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
    $TSD_counts{$row}++;
}
my @sorted_TSD_keys =  sort { $TSD_counts{$b} <=> $TSD_counts{$a} } keys (%TSD_counts);

#check whether the same TSD sequence was found in >80% of the copies. If not, count the occurances of TSD length, sort by highest occurances, and check if the length of the TSD is the same in >80% of the copies
my $final_TSD_length;
my $final_TSD_seq;
if (($TSD_counts{$sorted_TSD_keys[0]} / $TSD_array_length) > 0.8) {
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
}
#print "final TSD length: $final_TSD_length\n";

#get flanks to cluster by CD-Hit and also to get 10-bp on each side of the original TSD to look at target site preferences

my $flanks_cluster_path = $flanks_out_path . ".cluster";

system("cd-hit-est -i $flanks_out_path -o $flanks_cluster_path -c 0.8 -g 1 -T 0 -n 5 -d 0");
my $muscle_out_path = $flanks_out_path . ".msa";
system("muscle -in $flanks_out_path -out $muscle_out_path");
    

#-----------------------------------#
sub match_tirs {
    my $self = shift;
    #my $seq_obj = shift;
    my $input_path = shift;
    
    $self->throw("Need Bio::LocatableSeq argument")
        unless ref $self && $self->isa( 'Bio::LocatableSeq');
    if (!-f $input_path) {
        die "Supplied filepath is not valid";
    }
        
    my $seq = $self->seq();
    $seq =~ s/-//g;
    my $seq_name = $self->id();
    
    my $fa_aln_obj = Bio::SearchIO->new(-format => 'fasta', -file   => $input_path);
    my @result;
    
    
    #go through the FASTA input object . . down to the HSP
    while( my $result = $fa_aln_obj->next_result ) {
        while( my $hit = $result->next_hit ) {
            while( my $hsp = $hit->next_hsp ) {
                #grab the query, hit, and homology strings
                my $homo_string = $hsp->homology_string;
                my $query_str = $hsp->query_string;
                my $hit_str = $hsp->hit_string;
                my $len_homo = length($homo_string);
                my $len_query = length($query_str);
                my $len_hit = length($hit_str);
                #print "$query_str\t$len_query\n$homo_string\t$len_homo\n$hit_str\t$len_hit\n";
                
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
                    #print "count = $count\n";
                    my $homo_char =  substr($homo_string, $count, 1);
                    my $query_char =  substr($query_str, $count, 1);
                    my $hit_char =  substr($hit_str, $count, 1);
                    #print "$homo_char\n";
                    
                    if ($match_len == 0){
                        #if match length equals 0 and position is not a match, continue to next position
                        if ($homo_char eq " ") {
                            $mis_aln++;
                            if ($mis_aln >= 4) {
                                last;
                            }                            
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
                                push @result, (1, [@match]);
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
                    push @result, (0);
                }
            }
        }
    }
    return (@result);
}
