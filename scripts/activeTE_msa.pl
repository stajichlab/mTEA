#!/usr/bin/env perl

use strict;
use warnings;
use Bio::AlignIO;
use Bio::SearchIO;
use Bio::Seq;
use Bio::PrimarySeq;
use Bio::Tools::IUPAC;
use FindBin;
use File::Spec;
use Getopt::Long;

#set default flank length and reset if option is given
my $flank = 100;
my $all;
GetOptions (
'flank:i' => \$flank,
'all' => \$all
);
if (defined $all) {
    print "Intermediate files will be kept\n";
}
else{
    print "Intermediate files will be removed\n";
}

#make sure only one valid input file is specified
my $len = @ARGV;
if ($len != 1) {
    die "Please specify one input file\n";
}
my $infile = shift @ARGV;
if (!-f $infile) {
        die "Supplied filepath is not valid";
}
if ($flank <= 25) {
    die "Please use alignments with flanks  >25 nucleotides long";
}
 
#breakdown directory and filename, create output directory
my ($volume, $in_dir, $filename) = File::Spec->splitpath($infile);
my $out_dir_name = "activeTE_out_" . $filename;
my $out_path = File::Spec->catdir(($in_dir, $out_dir_name));
if (!-d $out_path) {
    system("mkdir $out_path");
}
my @fname_split = split('\\.', $filename);
my $fname_fin = $fname_split[0];
my %element_info;

#generate trimal alignment gap summary file, using sed to remove the header info
my $gap_summary_out = File::Spec->catpath($volume, $out_path, $filename . ".full_gap_summary");
system("trimal -in $infile -sgc | sed -n '4,\$p' > $gap_summary_out");

#initialize a global array to store results later
my @good_aln = ();
my $initial_left_TSD;
my $initial_right_TSD;

my @in_array;
open(my $in, "<", $infile);
while (my $line = <$in>){
    chomp $line;
    if ($line =~ /^([0-9]*).+ Query:(.*) Sbjct:(.*) Length.+Location:\(([0-9]*) \- ([0-9]*)\) Direction:(.+)/) {
        $line =~ s/ /_/g;
    }
    push @in_array, $line;
}
close($in);
open(my $fix_out, ">", $infile);
foreach my $line (@in_array) {
    print $fix_out "$line\n";
}
close($fix_out);

#create input object and an alignment object from it
my $in_obj = Bio::AlignIO->new(-file => $infile, -format => 'fasta');
my $full_aln_obj = $in_obj->next_aln();
my $full_aln_len = $full_aln_obj->length();

#generate a gap column matrix from the full alignment
my $gap_cols = $full_aln_obj->gap_col_matrix();

#Get the length of each sequence without gaps in the alignment and store in hash for later use
my %copy_lengths;
foreach my $seq_obj ($full_aln_obj->each_seq()) {
    my $seq = $seq_obj->seq();
    my $seq_id = $seq_obj->id();
    $seq =~ s/-//g;
    $copy_lengths{$seq_id} = length($seq);
}

#calculate the % nucleotide identity and the fraction of copies with sequence at each position the full alignment and print to file
my @full_id_array;
print "Starting Full ID calculation\n";
for (my $i = 1;  $i <= $full_aln_len; $i++ ) {
    my $pos = $full_aln_obj->slice($i,$i,1);
    my $pos_id = $pos->percentage_identity;
    my @gap_col_array = @{$gap_cols};
    my $gap_col_hashref = $gap_col_array[$i-1];
    my %gap_col_hash = %{$gap_col_hashref};
    my $total_count;
    my $base_count;
    foreach my $key (keys %gap_col_hash) {
        if ($gap_col_hash{$key} != 1) {
            $base_count++;
        }
        $total_count++;
    }
    my $pos_present = $base_count/$total_count;
    my @info = ($i, $pos_id, $pos_present);
    push @full_id_array, [@info];
}

my $full_id_out_path = File::Spec->catpath($volume, $out_path, $filename . ".full_id");
open(my $id_out, ">", $full_id_out_path);
foreach my $row_ref (@full_id_array) {
    my @pos = @{$row_ref};
    print $id_out "$pos[0]\t$pos[1]\t$pos[2]\n";
}
close($full_id_out_path);

#open and parse the trimal output file with the gap %  of each position in the full alignment
open(my $gap_id, "<", $gap_summary_out);
my @gap_id_array;
while (my $line = <$gap_id>){
    chomp $line;
    my @info = $line =~ /\S+\s+(\S+)\s+\S+/;
    push @gap_id_array, @info
}

#remove gap causing sequences from the full alignment based on high percentage of gaps at a position and either low sequence %ID or low fraction of sequences with sequence at that position
my %gap_seq_remove;
my %gap_seq_count;
my @gap_seq_pos_remove;

for (my $i = 0; $i < $full_aln_len; $i++){
    my $id_row_ref = $full_id_array[$i];
    my @id_info = @{$id_row_ref};
    my @gap_col_array = @{$gap_cols};
    my $gap_col_hashref = $gap_col_array[$i];
    my %gap_col_hash = %{$gap_col_hashref};
    if ($gap_id_array[$i] >= 90.0 and ($id_info[1] <= 40 or $id_info[2] <= 30)) {
        foreach my $key (keys %gap_col_hash) {
            if ($gap_col_hash{$key} != 1){
                my $seq_obj = $full_aln_obj->get_seq_by_id($key);
                my $seq = $seq_obj->seq();
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

my $trimmed_aln_obj = $full_aln_obj;
my @trim_id_array;
my $trim_gap_cols = $trimmed_aln_obj->gap_col_matrix();
my $trim_aln_len = $trimmed_aln_obj->length();
for (my $i = 1;  $i <= $trim_aln_len; $i++ ) {
    my $pos = $trimmed_aln_obj->slice($i,$i,1);
    my $pos_id = $pos->percentage_identity;
    my @gap_col_array = @{$trim_gap_cols};
    my $gap_col_hashref = $gap_col_array[$i-1];
    my %gap_col_hash = %{$gap_col_hashref};
    my $total_count;
    my $base_count;
    foreach my $key (keys %gap_col_hash) {
        if ($gap_col_hash{$key} != 1) {
            $base_count++;
        }
        $total_count++;
    }
    my $pos_present = $base_count/$total_count;
    my @info = ($i, $pos_id, $pos_present);
    push @trim_id_array, [@info];
}

my $trim_id_out_path = File::Spec->catpath($volume, $out_path, $filename . ".trim_id");
open(my $trim_id_out, ">", $trim_id_out_path);
foreach my $row_ref (@trim_id_array) {
    my @pos = @{$row_ref};
    print $trim_id_out "$pos[0]\t$pos[1]\t$pos[2]\n";
}
close($trim_id_out_path);

my $left_tir_start1 = 0;
my $left_tir_count = 0;
my $left_mm = 0;
for (my $i = 1;  $i < $trim_aln_len; $i++ ) {
    my $pos = $trimmed_aln_obj->slice($i,$i,1);
    my $pos_id = $pos->percentage_identity;
    if ($left_tir_count == 0) {
        if ($pos_id >= 80.0) {
            $left_tir_start1 = $i;
            $left_tir_count++;
            next;
        }
        else {
            next;
        }
    }
    elsif ($left_tir_count >= 1 and $left_tir_count <= 4) {
        if ($pos_id >= 80.0) {
            $left_tir_count++;
            next;
        }
        else {
            if ($left_mm <= 1) {
                $left_mm++;
                $left_tir_count++;
                next;
            }
            elsif ($left_mm > 1) {
                $left_tir_count = 0;
                $left_tir_start1 = 0;
                $left_mm = 0;
                next;
            }
        }    
    }
    elsif ($left_tir_count >=5) {
        last;
    }
}

my $right_tir_start1 = 0;
my $right_tir_count = 0;
my $right_mm = 0;
for (my $i = 0;  $i < $trim_aln_len; $i++ ) {
    my $pos = $trimmed_aln_obj->slice($trim_aln_len-$i,$trim_aln_len-$i,1);
    my $pos_id = $pos->percentage_identity;
    if ($right_tir_count == 0) {
        if ($pos_id >= 80.0) {
            $right_tir_start1 = $trim_aln_len-$i;
            $right_tir_count++;
            next;
        }
        else {
            next;
        }
    }
    elsif ($right_tir_count >= 1 and $right_tir_count <= 4) {
        if ($pos_id >= 80.0) {
            $right_tir_count++;
            next;
        }
        else {
            if ($right_mm <= 1) {
                $right_mm++;
                $right_tir_count++;
                next;
            }
            elsif ($right_mm > 1) {
                $right_tir_count = 0;
                $right_tir_start1 = 0;
                $right_mm = 0;
                next;
            }
        }    
    }
    elsif ($right_tir_count >= 5) {
        last;
    }
}

#check whether any sequences have gaps in the TIRSs and remove them
my %trim_gap_seq_remove;
my %trim_gap_seq_count;
my @trim_gap_seq_pos_remove;

for (my $i = 0; $i < $trim_aln_len; $i++){
    if (($i >= $left_tir_start1-1 and $i <= $left_tir_start1 + 49) or  ($i >= $right_tir_start1-51 and $i <= $right_tir_start1-1)) {
        my @trim_gap_col_array = @{$trim_gap_cols};
        my $trim_gap_col_hashref = $trim_gap_col_array[$i];
        my %trim_gap_col_hash = %{$trim_gap_col_hashref};
        my $pos_gaps = 0;
        my $pos_base = 0;
        foreach my $key (keys %trim_gap_col_hash) {
            if ($trim_gap_col_hash{$key} == 1){
                $pos_gaps++;
            }
            else {
                $pos_base++;
            }
        }
        if ($pos_gaps/($pos_gaps+$pos_base)>= 0.6) {
            foreach my $key (keys %trim_gap_col_hash) {
                if ($trim_gap_col_hash{$key} != 1){
                    my $seq_obj = $trimmed_aln_obj->get_seq_by_id($key);
                    my $seq = $seq_obj->seq();
                    my $seq_pos = $seq_obj->location_from_column($i+1);
                    my @info = ($key, $i+1, $seq_obj);
                    push @trim_gap_seq_pos_remove, [@info];
                    $trim_gap_seq_count{$key}++;
                    $trim_gap_seq_remove{$key} = $seq_obj;
                }
            }
        }
        else {
            foreach my $key (keys %trim_gap_col_hash) {
                if ($trim_gap_col_hash{$key} == 1){
                    my $seq_obj = $trimmed_aln_obj->get_seq_by_id($key);
                    my $seq = $seq_obj->seq();
                    my $seq_pos = $seq_obj->location_from_column($i+1);
                    my @info = ($key, $i+1, $seq_obj);
                    push @trim_gap_seq_pos_remove, [@info];
                    $trim_gap_seq_count{$key}++;
                    $trim_gap_seq_remove{$key} = $seq_obj;
                }
            }
        }
    }
}

foreach my $key (keys %trim_gap_seq_remove) {
    my $seq_obj = $trim_gap_seq_remove{$key};
    #print Dumper(\$seq_obj);
    $trimmed_aln_obj->remove_seq($seq_obj);
}
my $test_len = $trimmed_aln_obj->length();
if ($test_len == 0) {
    #open a file to store info on why analysis of an element was aborted
    my $bad_out_path = File::Spec->catpath($volume, $out_path, $filename . ".bad");
    open(my $bad_out, ">", $bad_out_path);
    print $bad_out "$filename\tTIRs not found in MSA\n";
    close($bad_out);
    print "TIRs not found in MSA, cleaning up files then exiting\n";
    if (!defined $all) {
        clean_files($out_path);
    }
    exit 0;
}

my $trim_aln_out2 = File::Spec->catpath($volume, $out_path, $filename . ".trim2");
my $out2 = Bio::AlignIO->new(-file => ">$trim_aln_out2", -format => 'fasta');
$out2->write_aln($trimmed_aln_obj);

#Store the column positions of the potential TIRs in the trimmed alignment and get the coulumns of them in the original alignminent to remove only sequences that cause gaps in the tirs
my %trim_pos_hash;
my %trim_left_pos_hash;
my %trim_right_pos_hash;
my %left_tir_start_check_counts;
my %right_tir_start_check_counts;

foreach my $seq_obj ($trimmed_aln_obj->each_seq()) {
    my $left_res_pos_obj = $seq_obj->location_from_column($left_tir_start1);
    my $right_res_pos_obj = $seq_obj->location_from_column($right_tir_start1);
    my $seq = $seq_obj->seq();
    my $seq_len = length($seq);
    my $seq_name = $seq_obj->id();
    my $left_res_pos = $left_res_pos_obj->start();
    my $left_res = substr($seq, $left_res_pos, 1);
    my $right_res_pos = $right_res_pos_obj->start();
    my $right_res = substr($seq, $right_res_pos, 1);
    my @trim_pos_array = ($left_res_pos, $right_res_pos);
    my $right_flank_pos = $seq_len - $right_res_pos;
    $trim_pos_hash{$seq_name} = @trim_pos_array;
    $trim_left_pos_hash{$seq_name} = $left_res_pos;
    $trim_right_pos_hash{$seq_name} = $right_res_pos;
    $left_tir_start_check_counts{$left_res_pos}++;
    $right_tir_start_check_counts{$right_flank_pos}++;
}

#sort the left and right TIR starts by largest count to smallest
my @sorted_left_tir_start_keys =  sort { $left_tir_start_check_counts{$b} <=> $left_tir_start_check_counts{$a} } keys (%left_tir_start_check_counts);
my @sorted_right_tir_start_keys =  sort { $right_tir_start_check_counts{$b} <=> $right_tir_start_check_counts{$a} } keys (%right_tir_start_check_counts);

#check if TIRs start in flanks, indicating the flanks are highly similar
my $left_flank_catch = 0;
my $right_flank_catch = 0;

if (!defined $sorted_left_tir_start_keys[0]) {
    if (!defined $sorted_right_tir_start_keys[0]) {
        #open a file to store info on why analysis of an element was aborted
        my $bad_out_path = File::Spec->catpath($volume, $out_path, $filename . ".bad");
        open(my $bad_out, ">", $bad_out_path);
        print $bad_out "$filename\tTIRs not found or both flanks similar\n";
        close($bad_out);
        if (!defined $all) {
            print "Cleaning up files\n";
            clean_files($out_path);
        }
        exit 0;
    }
    $left_flank_catch++;
    #open a file to store info on why analysis of an element was aborted
    my $bad_out_path = File::Spec->catpath($volume, $out_path, $filename . ".bad");
    open(my $bad_out, ">", $bad_out_path);
    print $bad_out "$filename\tLeft flank similar\n";
    close($bad_out);
    if (!defined $all) {
        print "Cleaning up files\n";
        clean_files($out_path);
    }
    exit 0;
}
elsif (!defined $sorted_right_tir_start_keys[0]) {
    #open a file to store info on why analysis of an element was aborted
    my $bad_out_path = File::Spec->catpath($volume, $out_path, $filename . ".bad");
    open(my $bad_out, ">", $bad_out_path);
    print $bad_out "$filename\tRight flank similar\n";
    close($bad_out);
    if (!defined $all) {
        print "Cleaning up files\n";
        clean_files($out_path);
    }
    exit 0;
}

if ($sorted_left_tir_start_keys[0] <= $flank-25) {
    $left_flank_catch++;
}
if ($sorted_right_tir_start_keys[0] <= $flank-25) {
    $right_flank_catch++;
}

if ($left_flank_catch != 0) {
    if ($right_flank_catch != 0) {
        #open a file to store info on why analysis of an element was aborted
        my $bad_out_path = File::Spec->catpath($volume, $out_path, $filename . ".bad");
        open(my $bad_out, ">", $bad_out_path);
        print $bad_out "$filename\tBoth flanks similar\n";
        close($bad_out);
        if (!defined $all) {
            print "Cleaning up files\n";
            clean_files($out_path);
        }
        exit 0;
    }
    else {
        #open a file to store info on why analysis of an element was aborted
        my $bad_out_path = File::Spec->catpath($volume, $out_path, $filename . ".bad");
        open(my $bad_out, ">", $bad_out_path);
        print $bad_out "$filename\tLeft flank similar\n";
        close($bad_out);
        if (!defined $all) {
            print "Cleaning up files\n";
            clean_files($out_path);
        }
        exit 0;
    }
}
elsif ($right_flank_catch != 0) {
    #open a file to store info on why analysis of an element was aborted
    my $bad_out_path = File::Spec->catpath($volume, $out_path, $filename . ".bad");
    open(my $bad_out, ">", $bad_out_path);
    print $bad_out "$filename\tRight flank similar\n";
    close($bad_out);
    if (!defined $all) {
        print "Cleaning up files\n";
        clean_files($out_path);
    }
    exit 0;
}

my @bad_remove;
my @bad_aln;
#grab all sequences from trimmed alignment
foreach my $seq_obj ( $trimmed_aln_obj->each_seq() ) {
    #grab sequence name and shorten it for use in some output filenames
    my $seq_name = $seq_obj->id();
    my @seqn_part = split(":", $seq_name);
    my $fname_start = $seqn_part[0];
    
    #grab sequence & strip hyphens
    my $seq = $seq_obj->seq();
    $seq =~ s/-//g;
    
    #initialize variables
    my $sub_seq;
    my $first;
    my $last;
    my $first_path = File::Spec->catpath($volume, $out_path, $filename . "first_half.txt");
    my $last_path = File::Spec->catpath($volume, $out_path, $filename . "last_half.txt");
    
    #get end  sequences
    $first = substr($seq, $trim_left_pos_hash{$seq_name}-1, 50);
    $last = substr($seq, $trim_right_pos_hash{$seq_name}-50,50);
    
    #save the two ends as files to use as inputs for a ggsearch search
    open(my $first_out, ">", $first_path) or die "Error creating $first_path. $!\n";
    open(my $last_out, ">", $last_path) or die "Error creating $last_path. $!\n";
    print $first_out ">first\n$first\n";
    print $last_out ">last\n$last\n";
    close($first_out);
    close($last_out);
    
    #create fasta output filename then call ggsearch
    my $out_opt = File::Spec->catpath($volume, $out_path, $fname_start . ".ggsearch.out");
    system("ggsearch36 -n -i -T 8 -d 1 $last_path $first_path > $out_opt");
    
    my @tir_match_result = match_tirs($seq_obj, $out_opt, 1);
    
    if ($tir_match_result[0] == 1) {
        push @good_aln, $tir_match_result[1];
    }
    else {
        push @bad_aln, $seq_name;
        push @bad_remove, $seq_obj;
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
my $good_aln_len = @good_aln;

if ($good_aln_len == 0) {
    print "No TIRs found from MSA: First round\n";
}
else {
    print "Found TIRs in MSA: 1st round\n";
}

#go through each entry of @good_aln array
foreach my $row_ref (@good_aln) {
    @entry = @{$row_ref};
    #get column position in alignment of left (hit) TIR, using sequence name and index position of TIR in full sequence. Increment the count of that column position in %hit_column_counts. Also, increment the count of TIR length in %hit_match_len
    $hit_aln_pos = $trimmed_aln_obj->column_from_residue_number( $entry[0], ${entry[1]{"hit"}}[0]);
    $hit_column_counts{$hit_aln_pos}++;
    $hit_match_len{${entry[1]{"hit"}}[1]}++;
    #do the same for the query (left) TIR
    $query_aln_pos = $trimmed_aln_obj->column_from_residue_number( $entry[0], ${entry[1]{"query"}}[0]);
    $query_column_counts{$query_aln_pos}++;
    $query_match_len{${entry[1]{"query"}}[1]}++;
}

#sort the hit and query column and match length hashes by largest count to smallest
my @sorted_hitcolumn_keys =  sort { $hit_column_counts{$b} <=> $hit_column_counts{$a} } keys (%hit_column_counts);
my @sorted_querycolumn_keys =  sort { $query_column_counts{$b} <=> $query_column_counts{$a} } keys (%query_column_counts);
my @sorted_hit_len_keys = sort { $hit_match_len{$b} <=> $hit_match_len{$a} } keys (%hit_match_len);
my @sorted_query_len_keys = sort { $query_match_len{$b} <=> $query_match_len{$a} } keys (%query_match_len);

#repeat above but shifted out 3bp on each end to catch cases where found TIR is slightly off
my %trim_pos_hash2;
my %trim_left_pos_hash2;
my %trim_right_pos_hash2;
my %left_tir_start_check_counts2;
my %right_tir_start_check_counts2;
my @good_aln2;

foreach my $seq_obj ($trimmed_aln_obj->each_seq()) {
    my $left_res_pos_obj = $seq_obj->location_from_column($left_tir_start1);
    my $right_res_pos_obj = $seq_obj->location_from_column($right_tir_start1);
    my $seq = $seq_obj->seq();
    my $seq_len = length($seq);
    my $left_res_pos = ($left_res_pos_obj->start()) - 3;
    my $left_res = substr($seq, $left_res_pos, 1);
    my $right_res_pos = ($right_res_pos_obj->start()) + 3;
    my $right_res = substr($seq, $right_res_pos, 1);
    my @trim_pos_array = ($left_res_pos, $right_res_pos);
    my $seq_name = $seq_obj->id();
    my $right_flank_pos = $seq_len - $right_res_pos;
    
    $trim_pos_hash2{$seq_name} = @trim_pos_array;
    $trim_left_pos_hash2{$seq_name} = $left_res_pos;
    $trim_right_pos_hash2{$seq_name} = $right_res_pos;
    $left_tir_start_check_counts2{$left_res_pos}++;
    $right_tir_start_check_counts2{$right_flank_pos}++;
}

#sort the left and right TIR starts by largest count to smallest
my @sorted_left_tir_start_keys2 =  sort { $left_tir_start_check_counts2{$b} <=> $left_tir_start_check_counts2{$a} } keys (%left_tir_start_check_counts2);
my @sorted_right_tir_start_keys2 =  sort { $right_tir_start_check_counts2{$b} <=> $right_tir_start_check_counts2{$a} } keys (%right_tir_start_check_counts2);

my @bad_remove2;
my @bad_aln2;
#grab all sequences from trimmed alignment
foreach my $seq_obj ( $trimmed_aln_obj->each_seq() ) {
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
    $first = substr($seq, $trim_left_pos_hash2{$seq_name}-1, 50);
    $last = substr($seq, $trim_right_pos_hash2{$seq_name}-50,50);
    
    #save the two ends as files to use as inputs for a ggsearch search
    open(my $first_out, ">", $first_path) or die "Error creating $first_path. $!\n";
    open(my $last_out, ">", $last_path) or die "Error creating $last_path. $!\n";
    print $first_out ">first\n$first\n";
    print $last_out ">last\n$last\n";
    close($first_out);
    close($last_out);
    
    #create fasta output filename then call ggsearch
    my $out_opt = File::Spec->catpath($volume, $out_path, $fname_start . ".ggsearch.out2");
    system("ggsearch36 -n -i -T 8 -d 1 $last_path $first_path > $out_opt");
    
    my @tir_match_result = match_tirs($seq_obj, $out_opt, 2);
    
    if ($tir_match_result[0] == 1) {
        push @good_aln2, $tir_match_result[1];
    }
    else {
        push @bad_aln2, $seq_name;
        push @bad_remove2, $seq_obj;
    }
}

#find the TIRs in the trimmed alignment using the index positions of each TIR stored above
#initialize variables
my %hit_column_counts2;
my %hit_match_len2;
my %query_column_counts2;
my %query_match_len2;
my $query_aln_pos2;
my $hit_aln_pos2;
my @entry2;

my $good_aln2_len = @good_aln2;

if ($good_aln2_len == 0) {
    print "No TIRs found from MSA: 2nd round\n";
}
else {
    print "Found TIRs in MSA: 2nd round\n";
}

#go through each entry of @good_aln2 array
foreach my $row_ref (@good_aln2) {
    @entry2 = @{$row_ref};
    #get column position in alignment of left (hit) TIR, using sequence name and index position of TIR in full sequence. Increment the count of that column position in %hit_column_counts. Also, increment the count of TIR length in %hit_match_len
    $hit_aln_pos2 = $trimmed_aln_obj->column_from_residue_number( $entry2[0], ${entry2[1]{"hit"}}[0]);
    $hit_column_counts2{$hit_aln_pos}++;
    $hit_match_len2{${entry2[1]{"hit"}}[1]}++;
    #do the same for the query (left) TIR
    $query_aln_pos2 = $trimmed_aln_obj->column_from_residue_number( $entry2[0], ${entry2[1]{"query"}}[0]);
    $query_column_counts2{$query_aln_pos}++;
    $query_match_len2{${entry2[1]{"query"}}[1]}++;
}

#sort the hit and query column and match length hashes by largest count to smallest
my @sorted_hitcolumn_keys2 =  sort { $hit_column_counts2{$b} <=> $hit_column_counts2{$a} } keys (%hit_column_counts2);
my @sorted_querycolumn_keys2 =  sort { $query_column_counts2{$b} <=> $query_column_counts2{$a} } keys (%query_column_counts2);
my @sorted_hit_len_keys2 = sort { $hit_match_len2{$b} <=> $hit_match_len2{$a} } keys (%hit_match_len2);
my @sorted_query_len_keys2 = sort { $query_match_len2{$b} <=> $query_match_len2{$a} } keys (%query_match_len2);

my $good_aln_len2 = @good_aln2;
my $bad_aln_len = @bad_aln;
my $bad_aln_len2 = @bad_aln2;

print "Good aln len:  $good_aln_len  Good aln len2: $good_aln_len2\n";
print "Bad aln len:  $bad_aln_len  Bad aln len2: $bad_aln_len2\n";

#check which run generated the longer TIRs
if ( ($good_aln_len2 > $good_aln_len) and ($bad_aln_len2 < $bad_aln_len) ){
    if ( ($sorted_hit_len_keys2[0] > $sorted_hit_len_keys[0] and $sorted_query_len_keys2[0] > $sorted_query_len_keys[0]) and ($sorted_hitcolumn_keys2[0] > $sorted_hitcolumn_keys[0] and $sorted_querycolumn_keys2[0] > $sorted_querycolumn_keys[0]) ) {
        print "Alignment set 2 better\n";
        @sorted_hitcolumn_keys = @sorted_hitcolumn_keys2;
        @sorted_querycolumn_keys = @sorted_querycolumn_keys2;
        @sorted_hit_len_keys = @sorted_hit_len_keys2;
        @sorted_query_len_keys = @sorted_query_len_keys2;
        %trim_pos_hash = %trim_pos_hash2;
        %trim_left_pos_hash = %trim_left_pos_hash2;
        %trim_right_pos_hash = %trim_right_pos_hash2;
        %left_tir_start_check_counts = %left_tir_start_check_counts2;
        %right_tir_start_check_counts = %right_tir_start_check_counts2;
        @good_aln = @good_aln2;
        @bad_aln = @bad_aln2;
        @bad_remove = @bad_remove2;
    }
    else {
        print "Alignment set 1 better\n\n";
    }
}
else {
    print "Alignment set 1 better\n\n";
}

#open a file to store info on why analysis of an element was aborted
my $bad_out_path = File::Spec->catpath($volume, $out_path, $filename . ".bad");
open(my $bad_out, ">", $bad_out_path);
foreach my $seq_obj (@bad_remove) {
    my $seq_name = $seq_obj->id();
    print $bad_out "$seq_name\tNo TIRs found\n";
    $trimmed_aln_obj->remove_seq($seq_obj);
}

#check whether potential TIRs were found
if ($left_tir_start1 == 0) {
    if ($right_tir_start1 == 0) {
        print $bad_out, "$filename\tTIRs not found in MSA";
        close($bad_out);
        if (!defined $all) {
            print "Cleaning up files\n";
            clean_files($out_path);
        }
        exit 0;
    }
    print $bad_out, "$filename\tLeft TIR not found";
    close($bad_out);
    if (!defined $all) {
        print "Cleaning up files\n";
        clean_files($out_path);
    }
    exit 0;
}
elsif ($right_tir_start1 == 0) {
    print $bad_out, "$filename\tRight TIR not found";
    close($bad_out);
    if (!defined $all) {
        print "Cleaning up files\n";
        clean_files($out_path);
    }
    exit 0;
}

#get the end residue positions for the tirs from one sequence in the trimmed alignment that will be used to find the tir interval in the full alignment. Then sequences that were excluded can be added back if they don't cause gaps in the tirs
my %tir_positions;
foreach my $seq_obj ($trimmed_aln_obj->each_seq()) {
    my $seq_name = $seq_obj->id();
    my $trim_seq = $seq_obj->seq();
    my $left_tir_start_pos = $trim_left_pos_hash{$seq_name};
    my $left_tir_end_pos =  $trim_left_pos_hash{$seq_name}+$sorted_hit_len_keys[0]-1;
    my $right_tir_start_pos = $trim_right_pos_hash{$seq_name};
    my $right_tir_end_pos_obj = $seq_obj->location_from_column($sorted_querycolumn_keys[0]);
    my $right_tir_end_pos = $right_tir_end_pos_obj->start();
    $tir_positions{$seq_name}{'left_tir_start'} = $left_tir_start_pos;
    $tir_positions{$seq_name}{'left_tir_end'} = $left_tir_end_pos;
    $tir_positions{$seq_name}{'right_tir_start'} = $right_tir_start_pos;
    $tir_positions{$seq_name}{'right_tir_end'} = $right_tir_end_pos;
}

#reread the original alignment file
my $in_obj2 = Bio::AlignIO->new(-file => $infile, -format => 'fasta');
my $ori_aln_obj = $in_obj2->next_aln();
my $ori_aln_len = $ori_aln_obj->length();

#get the column positions of tirs in the original alignment
my $left_tir_start;
my $left_tir_end;
my $right_tir_start;
my $right_tir_end;

foreach my $seq_obj ($ori_aln_obj->each_seq()) {
    my $seq_name = $seq_obj->id();
    if (!defined $tir_positions{$seq_name}) {
        next;
    }
    else {
        $left_tir_start = $ori_aln_obj->column_from_residue_number($seq_name, $tir_positions{$seq_name}{'left_tir_start'});
        $left_tir_end = $ori_aln_obj->column_from_residue_number($seq_name, $tir_positions{$seq_name}{'left_tir_end'});
        $right_tir_start = $ori_aln_obj->column_from_residue_number($seq_name, $tir_positions{$seq_name}{'right_tir_start'});
        $right_tir_end = $ori_aln_obj->column_from_residue_number($seq_name, $tir_positions{$seq_name}{'right_tir_end'});
        last;
    }
}
my $tir_length = $left_tir_end - $left_tir_start;

#go back through the full alignment and remove only the sequences that create gaps in the tirs or that failed the first tir search
my %gap_seq_remove2;
my %search_tirs;

for (my $i = 0; $i < $ori_aln_len; $i++){
    if (($i >= $left_tir_start-1 and $i <= $left_tir_end-1) or  ($i >= $right_tir_start-1 and $i <= $right_tir_end-1)) {
        my $id_row_ref = $full_id_array[$i];
        my @id_info = @{$id_row_ref};
        my @gap_col_array = @{$gap_cols};
        my $gap_col_hashref = $gap_col_array[$i];
        my %gap_col_hash = %{$gap_col_hashref};
        if ($gap_id_array[$i] >= 60.0 and $id_info[1] <= 40){
            foreach my $key (keys %gap_col_hash) {
                if ($gap_col_hash{$key} != 1){
                    $gap_seq_remove2{$key}++;
                }
            }
        }
        elsif ($gap_id_array[$i] <= 40.0 ) {
            foreach my $key (keys %gap_col_hash) {
                if ($gap_col_hash{$key} == 1){
                    $gap_seq_remove2{$key}++;
                }
            }
        }
    }
}

foreach my $row_ref (@gap_seq_pos_remove) {
    @entry = @{$row_ref};
    if (exists $gap_seq_remove2{$entry[0]}) {
        next;
    }
    else {
        $search_tirs{$entry[0]}++;
    }
}
foreach my $key (keys %gap_seq_remove) {
    if (exists $gap_seq_remove2{$key}) {
        next;
    }
    else {
        $search_tirs{$key}++;
    }
}

#keep track of removed sequences to print to file
my @removed_sequences;

foreach my $key (keys %gap_seq_remove2) {
    my $remove = $ori_aln_obj->get_seq_by_id($key);
    my $seq_id = $remove->id();
    push @removed_sequences, $seq_id;
    $ori_aln_obj->remove_seq($remove);
}

foreach my $seq_name (@bad_aln) {
    my $seq_obj = $ori_aln_obj->get_seq_by_id($seq_name);
    if (defined $seq_obj) {
        push @removed_sequences, $seq_name;
        $ori_aln_obj->remove_seq($seq_obj);
    }
}

#print list of removed sequences to file
my $removed_seq_path = File::Spec->catpath($volume, $out_path, $filename . ".removed_sequences.txt");
open(my $removed_seq, ">", $removed_seq_path);
my %removed_seq_hash;
foreach my $seq_id (@removed_sequences) {
    $removed_seq_hash{$seq_id}++;
    print $removed_seq "$seq_id\n";
}

my $final_aln_obj = $ori_aln_obj->remove_gaps('-',1);
my $int_aln_out = File::Spec->catpath($volume, $out_path, $filename . ".intermediate");
$out = Bio::AlignIO->new(-file => ">$int_aln_out", -format => 'fasta');
$out->write_aln($final_aln_obj);

#get the column positions of tirs in the intermediate alignment
foreach my $seq_obj ($final_aln_obj->each_seq()) {
    my $seq_name = $seq_obj->id();
    if (!defined $tir_positions{$seq_name}) {
        next;
    }
    else {
        $left_tir_start = $final_aln_obj->column_from_residue_number($seq_name, $tir_positions{$seq_name}{'left_tir_start'});
        $left_tir_end = $final_aln_obj->column_from_residue_number($seq_name, $tir_positions{$seq_name}{'left_tir_end'});
        $right_tir_start = $final_aln_obj->column_from_residue_number($seq_name, $tir_positions{$seq_name}{'right_tir_start'});
        $right_tir_end = $final_aln_obj->column_from_residue_number($seq_name, $tir_positions{$seq_name}{'right_tir_end'});
        last;
    }
}

my $new_gap_cols = $final_aln_obj->gap_col_matrix();
my $final_len = $final_aln_obj->length();
my %new_gap_seq_remove;

for (my $i = 0; $i < $final_len; $i++){
    if (($i >= $left_tir_start-1 and $i <= $left_tir_end-1) or  ($i >= $right_tir_start-1 and $i <= $right_tir_end-1)) {
        my @trim_gap_col_array = @{$new_gap_cols};
        my $trim_gap_col_hashref = $trim_gap_col_array[$i];
        my %trim_gap_col_hash = %{$trim_gap_col_hashref};
        my $pos_gaps = 0;
        my $pos_base = 0;
        foreach my $key (keys %trim_gap_col_hash) {
            if ($trim_gap_col_hash{$key} == 1){
                $pos_gaps++;
            }
            else {
                $pos_base++;
            }
        }
        if ($pos_gaps/($pos_gaps+$pos_base)> 0.5) {
            foreach my $key (keys %trim_gap_col_hash) {
                if ($trim_gap_col_hash{$key} != 1){
                    my $seq_obj = $final_aln_obj->get_seq_by_id($key);
                    $new_gap_seq_remove{$key} = $seq_obj;
                }
            }
        }
        else {
            foreach my $key (keys %trim_gap_col_hash) {
                if ($trim_gap_col_hash{$key} == 1){
                    my $seq_obj = $final_aln_obj->get_seq_by_id($key);
                    $new_gap_seq_remove{$key} = $seq_obj;
                }
            }
        }
    }
}

foreach my $key (keys %new_gap_seq_remove) {
    my $seq_obj = $new_gap_seq_remove{$key};
    $final_aln_obj->remove_seq($seq_obj);
}
my $check_len = $final_aln_obj->length();
$final_aln_obj = $final_aln_obj->remove_gaps('-',1);
my %left_tsd_info;
my %right_tsd_info;

foreach my $seq_name (keys %search_tirs) {
    if (exists $gap_seq_remove2{$seq_name}) {
        next;
    }
    if (exists $removed_seq_hash{$seq_name}) {
        next;
    }
    if (exists $new_gap_seq_remove{$seq_name}) {
        next;
    }
    my $seq_obj = $final_aln_obj->get_seq_by_id($seq_name);
    my @seqn_part = split(":", $seq_name);
    my $fname_start = $seqn_part[0];
    my $seq = $seq_obj->seq();
    my $seq_len = length($seq);
    if ($seq =~ m/NNNN/g) {
        $final_aln_obj->remove_seq($seq_obj);
        next;
    }
    my $left_test_pos = substr($seq, $left_tir_start-1, 1);
    my $right_test_pos = substr($seq, $right_tir_start-1, 1);
    if ($left_test_pos eq "-"  or $right_test_pos eq "-") {
        $final_aln_obj->remove_seq($seq_obj);
        next;
    }
    
    my $left_tir_start_obj = $seq_obj->location_from_column($left_tir_start);
    my $left_tir_start_pos = $left_tir_start_obj->start();
    my $right_tir_start_obj = $seq_obj->location_from_column($right_tir_start);
    my $right_tir_start_pos = $right_tir_start_obj->start();
    
    $seq =~ s/-//g;
    my $first;
    my $last;
    my $first_path = File::Spec->catpath($volume, $out_path, "first_half.txt");
    my $last_path = File::Spec->catpath($volume, $out_path, "last_half.txt");    
    $first = substr($seq, $left_tir_start_pos-1, 50);
    $last = substr($seq, $right_tir_start_pos-50,50);
    
    #grab TSD and initial flanks now too
    my $starting_left_flank = substr($seq, 0, $left_tir_start_pos);
    my $left_tsd = substr($seq, $left_tir_start_pos-21, 22);
    my $left_tsd_substr = substr($left_tsd, 10, 10);
    my $starting_right_flank = substr($seq, $right_tir_start_pos-2);
    my $right_tsd = substr($seq, $right_tir_start_pos-2, 22);
    my $right_tsd_substr = substr($right_tsd, 2, 10);
    
    #continue TIR search
    open(my $first_out, ">", $first_path) or die "Error creating $first_path. $!\n";
    open(my $last_out, ">", $last_path) or die "Error creating $last_path. $!\n";
    print $first_out ">first\n$first\n";
    print $last_out ">last\n$last\n";
    close($first_out);
    close($last_out);
    my $out_opt = File::Spec->catpath($volume, $out_path, $fname_start . ".fasta.out");
    system("ggsearch36 -n -i -T 8 -d 1 $last_path $first_path > $out_opt");
    
    my @tir_match_result = match_tirs($seq_obj, $out_opt, 1);
    if ($tir_match_result[0] == 0) {
        $final_aln_obj->remove_seq($seq_obj);
    }
    else {
        $left_tsd_info{$seq_name} = [$left_tsd, $left_tsd_substr, $starting_left_flank];
        $right_tsd_info{$seq_name} = [$right_tsd, $right_tsd_substr, $starting_right_flank];
    }
}

$final_aln_obj = $final_aln_obj->remove_gaps('-',1);
my $int_aln_out2 = File::Spec->catpath($volume, $out_path, $filename . ".intermediate2");
$out = Bio::AlignIO->new(-file => ">$int_aln_out2", -format => 'fasta');
$out->write_aln($final_aln_obj);
my $last_len = $final_aln_obj->length();

#get the column positions of tirs in the final alignment
foreach my $seq_obj ($final_aln_obj->each_seq()) {
    my $seq_name = $seq_obj->id();
    if (!defined $tir_positions{$seq_name}) {
        next;
    }
    else {
        $left_tir_start = $final_aln_obj->column_from_residue_number($seq_name, $tir_positions{$seq_name}{'left_tir_start'});
        $left_tir_end = $final_aln_obj->column_from_residue_number($seq_name, $tir_positions{$seq_name}{'left_tir_end'});
        $right_tir_start = $final_aln_obj->column_from_residue_number($seq_name, $tir_positions{$seq_name}{'right_tir_start'});
        $right_tir_end = $final_aln_obj->column_from_residue_number($seq_name, $tir_positions{$seq_name}{'right_tir_end'});
        last;
    }
}

#Extract the left and right TIRs as new alignment objects, starting at the most common starting column until that column plus the most common TIR length
my $left_TIR_aln_obj = $final_aln_obj->slice($left_tir_start, $left_tir_end,1);
my $right_TIR_aln_obj = $final_aln_obj->slice($right_tir_end, $right_tir_start,1);
my $element_aln_obj = $final_aln_obj->slice($left_tir_start,$right_tir_start,1);

my $gap_summary_out2 = File::Spec->catpath($volume, $out_path, $filename . ".final_gap_summary");
system("trimal -in $int_aln_out2 -sgc | sed -n '4,\$p' > $gap_summary_out2");

#open and parse the trimal output file with the % gap ID of each position in the full alignment
open(my $gap_id2, "<", $gap_summary_out2);
my @gap_id_array2;
while (my $line = <$gap_id2>){
    chomp $line;
    my @info = $line =~ /\S+\s+(\S+)\s+\S+/;
    #my @info = split(/\s+/, $line);
    push @gap_id_array2, @info
}
#generate a gap column matrix from the full alignment
my $final_gap_cols = $final_aln_obj->gap_col_matrix();

my @final_id_array;
$final_len = $final_aln_obj->length();
#calculate the % identity of each position in the full alignment and print to file
print "Starting Final ID calculation\n";
for (my $i = 1;  $i <= $final_len; $i++ ) {
    my $pos = $final_aln_obj->slice($i,$i,1);
    my $pos_id = $pos->percentage_identity;
    my @info = ($i, $pos_id);
    push @final_id_array, [@info];
}
print "Finished Final ID calculation\n";

#remove gap causing sequences in the TIRs from the full alignment based on high % gap ID and low sequence % ID
my %final_gap_seq_remove;
my %final_gap_seq_count;
my @final_gap_seq_pos_remove;

for (my $i = 0; $i < $last_len; $i++){
    for (($i >= $left_tir_start-1 and $i <= $left_tir_end-1) or ($i >= $right_tir_end-1 and $i <= $right_tir_start-1)) {
        my $id_row_ref = $final_id_array[$i];
        my @id_info = @{$id_row_ref};
        my @gap_col_array = @{$final_gap_cols};
        my $gap_col_hashref = $gap_col_array[$i];
        my %gap_col_hash = %{$gap_col_hashref};
        if ($gap_id_array2[$i] >= 90.0 and $id_info[1] <= 40){
            foreach my $key (keys %gap_col_hash) {
                if ($gap_col_hash{$key} != 1){
                    my $seq_obj = $final_aln_obj->get_seq_by_id($key);
                    my $seq = $seq_obj->seq();
                    my $seq_pos = $seq_obj->location_from_column($i+1);
                    my @info = ($key, $i+1, $seq_obj);
                    push @final_gap_seq_pos_remove, [@info];
                    $final_gap_seq_count{$key}++;
                    $final_gap_seq_remove{$key} = $seq_obj;
                }
            }
        }
        elsif ($gap_id_array2[$i] <= 5 and $id_info[1] >= 90) {
            foreach my $key (keys %gap_col_hash) {
                if ($gap_col_hash{$key} == 1){
                    my $seq_obj = $final_aln_obj->get_seq_by_id($key);
                    my $seq = $seq_obj->seq();
                    my $seq_pos = $seq_obj->location_from_column($i+1);
                    my @info = ($key, $i+1, $seq_obj);
                    push @final_gap_seq_pos_remove, [@info];
                    $final_gap_seq_count{$key}++;
                    $final_gap_seq_remove{$key} = $seq_obj;
                }
            }
        }
    }
}
foreach my $key (keys %final_gap_seq_remove) {
    my $seq_obj = $final_gap_seq_remove{$key};
    $final_aln_obj->remove_seq($seq_obj);
}
my $test_len2 = $final_aln_obj->length();
$final_aln_obj = $final_aln_obj->remove_gaps('-',1);
my $int_aln_out3 = File::Spec->catpath($volume, $out_path, $filename . ".intermediate3");
$out = Bio::AlignIO->new(-file => ">$int_aln_out3", -format => 'fasta');
$out->write_aln($final_aln_obj);

#get the column positions of tirs in the final alignment
foreach my $seq_obj ($final_aln_obj->each_seq()) {
    my $seq_name = $seq_obj->id();
    if (!defined $tir_positions{$seq_name}) {
        next;
    }
    else {
        $left_tir_start = $final_aln_obj->column_from_residue_number($seq_name, $tir_positions{$seq_name}{'left_tir_start'});
        $left_tir_end = $final_aln_obj->column_from_residue_number($seq_name, $tir_positions{$seq_name}{'left_tir_end'});
        $right_tir_start = $final_aln_obj->column_from_residue_number($seq_name, $tir_positions{$seq_name}{'right_tir_start'});
        $right_tir_end = $final_aln_obj->column_from_residue_number($seq_name, $tir_positions{$seq_name}{'right_tir_end'});
        last;
    }
}

#find TSDs based off column # of TIRs. Must check for 2-bp TSD included in TIRS and then a longer TSD. Longer TSDs that start and finish with same 2-bp as the 2-bp TSD supersede the 2-bp TSD. Also, grab flanks and save them to files to be used by CD-Hit to cluster them.
my @putative_TSD;
my @TSD_info;
my @no_TSD_found;
my @put_TSD_names;
my $flanks_out_path = File::Spec->catpath($volume, $out_path, $filename . "_flanks.fa");
open(my $flanks_out, ">", $flanks_out_path) or die "Error creating $flanks_out_path. $!\n";
my $tsd1_count = 0;
my $tsd2_count = 0;

print "Starting code to find TSDs\n";
foreach my $seq_obj ($final_aln_obj->each_seq()) {
    my $seq_name = $seq_obj->id();
    my @seqn_part = split(":", $seq_name);
    my $fname_start = $seqn_part[0];
    my $seq_ori = $seq_obj->seq();
    my $seq = $seq_ori;
    $seq =~ s/-//g;
    my $left_tsd;
    my $left_tsd_substr;
    my $starting_left_flank;
    my $right_tsd;
    my $right_tsd_substr;
    my $starting_right_flank;
    my $last = 0;
    
    #get tsd and flank sequences from hashes if exists
    if (exists $left_tsd_info{$seq_name} and exists $right_tsd_info{$seq_name}) {
        my @left_tsd_info = @{$left_tsd_info{$seq_name}};
        $left_tsd = $left_tsd_info[0];
        $left_tsd_substr = $left_tsd_info[1];
        $starting_left_flank = $left_tsd_info[2];
    
        my @right_tsd_info = @{$right_tsd_info{$seq_name}};
        $right_tsd = $right_tsd_info[0];
        $right_tsd_substr = $right_tsd_info[1];
        $starting_right_flank = $right_tsd_info[2];
    }
    else {
        #get seq_pos of the start of TIRs from column # in alignment. Adjust to include the 2-bp at outside ends of TIRs and then split into the 2-bp and 10-bp sequences to search for TSDs
        
        my $left_tsd_loc_obj = $seq_obj->location_from_column($left_tir_start);
        my $left_tsd_end_pos = $left_tsd_loc_obj->start();
        $starting_left_flank = substr($seq, 0, $left_tsd_end_pos);
        $left_tsd = substr($seq, $left_tsd_end_pos-21, 22);
        $left_tsd_substr = substr($left_tsd, 10, 10);
        
        my $right_tsd_loc_obj = $seq_obj->location_from_column($right_tir_start);
        my $right_tsd_start_pos = $right_tsd_loc_obj->start();
        $starting_right_flank = substr($seq, $right_tsd_start_pos-2);
        $right_tsd = substr($seq, $right_tsd_start_pos-2, 22);
        $right_tsd_substr = substr($right_tsd, 2, 10);
    }
    
    my $left_2bp = substr($left_tsd,-2);
    my $right_2bp = substr($right_tsd, 0, 2);
    my $tsd1;
    my $tsd2;
    if ($left_2bp eq $right_2bp) {
        print "Round 1 TSD search - Yes!\n";
        $tsd1 = $left_2bp;
    }
    
    for (my $i = 0; $i < 9; $i++ ){
        if (substr($left_tsd_substr, $i) eq substr($right_tsd_substr, 0, -($i))) {
            print "Round 2 TSD search - Yes! $i\n";
            $tsd2 = substr($left_tsd_substr, $i);
            next;
        }
        else{
            my $sub_print = substr($left_tsd_substr,$i);
            next;
        }
    }   
    #Save found TSD to an array or report that a TSD wasn't found
    if ($tsd1) {
        if (! $tsd2) {
            my $insertion_site = substr($left_tsd, -12) . substr($right_tsd, 2, 10);
            push @TSD_info, [$seq_name, $insertion_site];            
            push @putative_TSD, $tsd1;
            push @put_TSD_names, [$seq_name, $tsd1];
            my $left_flank = substr($starting_left_flank, 0, -2);
            my $right_flank = substr($starting_right_flank, 2);
            my $flanks = $left_flank . $right_flank;
            print $flanks_out ">$fname_start\n$flanks\n";
            $tsd1_count++;
        }
        else {
            #if both a 2-bp TSD and a longer one is found, check if the longer one has the same 2-bp (as the found 2-bp TSD) at both ends. If so, use the longer TSD. Otherwise use the 2-bp TSD
            if ((length($tsd2) > length($tsd1)) and (substr($tsd2, 0, 2) eq $tsd1) and (substr($tsd2, -2) eq $tsd1))  {
                my $insertion_site = substr($left_tsd, (-2-length($tsd2)-10), (length($tsd2)+10)) . substr($right_tsd, 2+length($tsd2), 10);
                push @TSD_info, [$seq_name, $insertion_site];
                push @putative_TSD, $tsd2;
                push @put_TSD_names, [$seq_name, $tsd2];
                my $left_flank = substr($starting_left_flank, 0, (-2-length($tsd2)));
                my $right_flank = substr($starting_right_flank, (2+length($tsd2)));
                my $flanks = $left_flank . $right_flank;
                print $flanks_out ">$fname_start\n$flanks\n";
                $tsd2_count++;
            }
            else {
                my $insertion_site = substr($left_tsd, -12) . substr($right_tsd, 2, 10);
                push @TSD_info, [$seq_name, $insertion_site];
                push @putative_TSD, $tsd1;
                push @put_TSD_names, [$seq_name, $tsd1];
                my $left_flank = substr($starting_left_flank, 0, -2);
                my $right_flank = substr($starting_right_flank, 2);
                my $flanks = $left_flank . $right_flank;
                print $flanks_out ">$fname_start\n$flanks\n";
                $tsd1_count++;
            }
        }
    }
    else {
        if ($tsd2) {
            my $insertion_site = substr($left_tsd, (-2-length($tsd2)-10), (length($tsd2)+10)) . substr($right_tsd, 2+length($tsd2), 10);
            push @TSD_info, [$seq_name, $insertion_site];
            push @putative_TSD, $tsd2;
            push @put_TSD_names, [$seq_name, $tsd2];
            my $left_flank = substr($starting_left_flank, 0, (-2-length($tsd2)));
            my $right_flank = substr($starting_right_flank, (2+length($tsd2)));
            my $flanks = $left_flank . $right_flank;
            print $flanks_out ">$fname_start\n$flanks\n";
            $tsd2_count++;
        }
        else {
            #store or print to report that TSD was not found for this copy use $seq_name, etc
            #use the end position of tirs to grab the flanking sequences
            my $left_flank = substr($starting_left_flank, 0, -2);
            my $right_flank = substr($starting_right_flank, 2);
            my $flanks = $left_flank . $right_flank;
            push @no_TSD_found, $seq_name;
            print $flanks_out ">$fname_start\n$flanks\n";
        }
    }
}
close($flanks_out);
print "Finished code to find TSDs\n";

my $final_align_len = $final_aln_obj->num_sequences();
$element_info{"copy_num"} = $final_align_len;
my $final_aln_out = File::Spec->catpath($volume, $out_path, $filename . ".final");
$out = Bio::AlignIO->new(-file => ">$final_aln_out", -format => 'fasta');
$out->write_aln($final_aln_obj);

my $no_TSD_found_len = @no_TSD_found;
if ($no_TSD_found_len >= 1) {
    print "Copies without TSDs removed and printed to file\n";
    my $no_TSD_found_out_path = File::Spec->catpath($volume, $out_path, $filename . ".no_TSD_found.info");
    open(my $no_TSD_found_out, ">", $no_TSD_found_out_path);
    foreach my $item (@no_TSD_found) {
        print $no_TSD_found_out "$item\n";
    }
    close($no_TSD_found_out);
}

#if necessary, change TIRS to remove the 2bp TSD that are included in them & calculate the overall percent identity of each TIR
if ($tsd1_count > $tsd2_count) {
    print "Adjusting TIRs by 2bp\n";
    $left_TIR_aln_obj = $final_aln_obj->slice($left_tir_start+2, $left_tir_end,1);
    $right_TIR_aln_obj = $final_aln_obj->slice($right_tir_end, $right_tir_start-2,1);
    $element_aln_obj = $final_aln_obj->slice($left_tir_start+2,$right_tir_start-2,1);
    foreach my $seq_obj ($final_aln_obj->each_seq()) {
        my $seq_name = $seq_obj->id();
        if (!defined $tir_positions{$seq_name}{'left_tir_start'}) {
            my $left_pos_obj = $seq_obj->location_from_column($left_tir_start);
            my $right_pos_obj = $seq_obj->location_from_column($right_tir_start);
            my $left_pos = $left_pos_obj->start();
            my $right_pos = $right_pos_obj->start();
            $tir_positions{$seq_name}{'left_tir_start'} = $left_pos;
            $tir_positions{$seq_name}{'right_tir_start'} = $right_pos;
        }
        $element_info{$seq_name}{"left_tir_start"} = $tir_positions{$seq_name}{'left_tir_start'} + 2;
        $element_info{$seq_name}{"right_tir_start"} = $tir_positions{$seq_name}{'right_tir_start'} - 2;
    }
}
else {
    foreach my $seq_obj ($final_aln_obj->each_seq()) {
        my $seq_name = $seq_obj->id();
        if (!defined $tir_positions{$seq_name}{'left_tir_start'}) {
            my $left_pos_obj = $seq_obj->location_from_column($left_tir_start);
            my $right_pos_obj = $seq_obj->location_from_column($right_tir_start);
            my $left_pos = $left_pos_obj->start();
            my $right_pos = $right_pos_obj->start();
            $tir_positions{$seq_name}{'left_tir_start'} = $left_pos;
            $tir_positions{$seq_name}{'right_tir_start'} = $right_pos;
        }
        $element_info{$seq_name}{"left_tir_start"} = $tir_positions{$seq_name}{'left_tir_start'};
        $element_info{$seq_name}{"right_tir_start"} = $tir_positions{$seq_name}{'right_tir_start'};
    }
}
my $left_tir_id = $left_TIR_aln_obj->percentage_identity;
my $right_tir_id = $right_TIR_aln_obj->percentage_identity;
my $element_id = $element_aln_obj->percentage_identity;
my $left_tir_seq = $left_TIR_aln_obj->consensus_string();
my $right_tir_seq = $right_TIR_aln_obj->consensus_string();

$element_info{"element_id"} = $element_id;
$element_info{"left_tir_seq"} = $left_tir_seq;
$element_info{"left_tir_id"} = $left_tir_id;
$element_info{"right_tir_seq"} = $right_tir_seq;
$element_info{"right_tir_id"} = $right_tir_id;

#Determine the most common TSD by sequence first if possible or by length. If >80% of the TSDs are the same, then that sequence is stored as the TSD for output. Otherwise, look at the lengths of the TSDs and store if >80% of the TSDs have the same length.
my %TSD_counts;
my $TSD_array_length = @put_TSD_names;
if ($TSD_array_length == 0) {
    print "No TSDs found, cleaning up files then exiting program\n";
    print $bad_out "$filename\tNo TSDs found for any copy";
    close($bad_out);
    if (!defined $all) {
        clean_files($out_path);
    } 
    exit 0;
}

my $insertion_site_file = $filename . ".insertion-site.info";
my $insertion_site_out_path = File::Spec->catpath($volume, $out_path, $insertion_site_file);
open(my $insertion_site_out, ">", $insertion_site_out_path) or die "Error creating $insertion_site_out_path. $!\n";
my $insertion_num = @TSD_info;
foreach my $row_ref (@TSD_info) {
    print "Printing insertion site info\n";
    my @pos = @{$row_ref};
    print $insertion_site_out ">$pos[0]\n$pos[1]\n";
}
close($insertion_site_out);

#count the occurances of a TSD sequence and sort with the highest value first in an array
foreach my $row_ref (@put_TSD_names) {
    my @pos = @{$row_ref};
    $TSD_counts{$pos[1]}++;
    $element_info{$pos[0]}{"TSD"} = $pos[1];
}
my @sorted_TSD_keys =  sort { $TSD_counts{$b} <=> $TSD_counts{$a} } keys (%TSD_counts);

#check whether the same TSD sequence was found in >80% of the copies. If not, count the occurances of TSD length, sort by highest occurances, and check if the length of the TSD is the same in >80% of the copies
my $final_TSD_length;
my $final_TSD_seq;
my $TSD_fraction = $TSD_counts{$sorted_TSD_keys[0]} / $TSD_array_length;
my $need_consensus = 0;

if (($TSD_fraction) > 0.8) {
    $element_info{"TSD_fraction"} = $TSD_fraction;
    $final_TSD_length = length($sorted_TSD_keys[0]);
    $final_TSD_seq = $sorted_TSD_keys[0];
    $element_info{"TSD_seq"} = $final_TSD_seq;
    $element_info{"TSD_len"} = $final_TSD_length;
    $element_info{"TSD_fraction"} = $TSD_fraction;
}
else {
    $need_consensus = 1;
    my %TSD_length_counts;
    foreach my $row (@putative_TSD) {
        $TSD_length_counts{length($row)}++;
    }
    my @sorted_TSD_length_keys =  sort { $TSD_length_counts{$b} <=> $TSD_length_counts{$a} } keys (%TSD_length_counts);
    $final_TSD_length = $sorted_TSD_length_keys[0];
    $TSD_fraction = $TSD_array_length/$final_align_len;
    
    $element_info{"TSD_len"} = $final_TSD_length;
    $element_info{"TSD_fraction"} = $TSD_fraction;
}

my $tsd_out_path = File::Spec->catpath($volume, $out_path, $filename . "_tsd.fa");
open(my $tsd_info_out, ">", $tsd_out_path) or die "Error creating $tsd_out_path. $!\n";

foreach my $row_ref (@put_TSD_names) {
    my @pos = @{$row_ref};
    if (length($pos[1]) == $final_TSD_length) {
        print $tsd_info_out ">$pos[0]\n$pos[1]\n";
    }
}
close($tsd_info_out);

if ($need_consensus == 1) {
    #read in tsd file as MSA to generate IUPAC consensus
    my $tsd_in_obj = Bio::AlignIO->new(-file => $tsd_out_path, -format => 'fasta');
    my $tsd_aln_obj = $tsd_in_obj->next_aln();
    my $tsd_consensus = $tsd_aln_obj->consensus_iupac();
    $element_info{"TSD_seq"} = $tsd_consensus;
}

#import TE characteristics from table and generate regular expressions to classify element by TIR and TSD if possible
print "Processing TE charateristics table\n";
my $TE_char_path = $FindBin::Bin . "/DNA_TE_TIR-TSD.table";
open(my $TE_in, "<", $TE_char_path) or die "Error reading $TE_char_path . $!";
my %element_char_hash;

while(my $line = <$TE_in>) {
    chomp $line;
    my @split = split("\t", $line);
    my $ele_name = $split[0];
    if ($split[1] ne ".") {
        my $seq = $split[1];
        my $seq_obj = Bio::Seq->new(-seq => $seq, -alphabet => 'dna');
        my $iupac = Bio::Tools::IUPAC->new(-seq => $seq_obj);
        my $regexp = $iupac->regexp();
        $element_char_hash{$ele_name}{'tir_con'} = $regexp;
    }
    else {
        $element_char_hash{$ele_name}{'tir_con'} = '';
    }
    if ($split[2] =~ m/,/) {
        my @tsd_lengths = split(',', $split[2]);
        my $tsd_array_len = @tsd_lengths;
        my $count = 1;
        my $regexp = '';
        foreach my $tsd_len (@tsd_lengths) {
            if ($count < $tsd_array_len) {
                $regexp = $regexp . $tsd_len . "|";
                $count++;
            }
            else {
                $regexp = $regexp . $tsd_len;
            }
        }
        $element_char_hash{$ele_name}{'tsd_length'} = $regexp;
    }
    else {
        $element_char_hash{$ele_name}{'tsd_length'} = $split[2];
    }
    if ($split[3] ne ".") {
        if ($split[3] =~ m/,/) {
            my @tsd_cons = split(',', $split[3]);
            my $tsd_cons_len = @tsd_cons;
            my $count = 1;
            my $regexp = '';
            foreach my $tsd_consensus (@tsd_cons) {
                my $seq_obj = Bio::Seq->new(-seq => $tsd_consensus, -alphabet => 'dna');
                my $iupac = Bio::Tools::IUPAC->new(-seq => $seq_obj);
                my $regexp_ori = $iupac->regexp();
                if ($count < $tsd_cons_len) {
                    $regexp = $regexp . $regexp_ori . "|";
                    $count++;
                }
                else {
                    $regexp = $regexp . $regexp_ori;
                }
            }
            $element_char_hash{$ele_name}{'tsd_con'} = $regexp;
        }
        else {
            my $seq_obj = Bio::Seq->new(-seq => $split[3], -alphabet => 'dna');
            my $iupac = Bio::Tools::IUPAC->new(-seq => $seq_obj);
            my $regexp = $iupac->regexp();
            $element_char_hash{$ele_name}{'tsd_con'} = $regexp;
        }
    }
    else {
        $element_char_hash{$ele_name}{'tsd_con'} = '';
    }
}
print "Checking element for TE charateristics\n";
my %element_hits;
foreach my $ele_name (keys %element_char_hash) {
    if ($element_info{'TSD_len'} =~ m/$element_char_hash{$ele_name}{"tsd_length"}/) {
        $element_hits{$ele_name}++;
        if ($element_char_hash{$ele_name}{"tsd_con"} ne '' and $element_info{'TSD_seq'} =~ m/$element_char_hash{$ele_name}{"tsd_con"}/i) {
            $element_hits{$ele_name}++;
        }
        if ($element_char_hash{$ele_name}{"tir_con"} ne '') {
            if ($element_info{'left_tir_seq'} =~ m/^$element_char_hash{$ele_name}{"tir_con"}/i) {
                $element_hits{$ele_name}++;
            }
            else {
                delete $element_hits{$ele_name};
            }
        }
    }
}
my @sorted;
foreach my $key ( sort { $element_hits{$b} <=> $element_hits{$a} } keys (%element_hits)){
    my $count = $element_hits{$key};
    push @sorted, [$key, $count];
} 

my $classification = '';
foreach my $row_ref (@sorted) {
    my @info = @{$row_ref};
    if ($info[1] == ${$sorted[0]}[1]) {
        if ($classification eq '') {
            $classification = $classification . $info[0] . "_" . $info[1];
        }
        else {
            $classification = $classification . ", " . $info[0] . "_" . $info[1];
        }
    }
    else {
        last;
    }
}
if ($classification eq '') {
    $classification = "Unknown";
}
$element_info{"classification"} = $classification;
print "Element classification finished\n";
my $element_info_out_path = File::Spec->catpath($volume, $out_path, $filename . ".element_info");
open(my $element_info_out, ">", $element_info_out_path) or die "Error creating $element_info_out_path. $!\n";
print $element_info_out "$fname_fin\t$element_info{'copy_num'}\t$element_id\t$element_info{'left_tir_seq'}\t$element_info{'left_tir_id'}\t$element_info{'right_tir_seq'}\t$element_info{'left_tir_id'}\t$element_info{'TSD_len'}\t$element_info{'TSD_seq'}\t$element_info{'TSD_fraction'}\t$classification";
close($element_info_out);

my $out_fix = $out_path . "/";
my $Blogo_config_path = $FindBin::Bin . "/blogo/Blogo.conf";
system("Blogo_batch.pl file_path=$out_fix file_names=$insertion_site_file img_abs_dir=$out_fix conf_file=$Blogo_config_path");

print "Setting up gff path\n";
my $gff_path = File::Spec->catpath($volume, $out_path, $filename . ".gff");
print "Call gff printer\n";
generate_gff($final_aln_obj, $gff_path, "final", \%element_info);
print "Exited gff printer\n";

if (!defined $all) {
    print "Cleaning up files\n";
    clean_files($out_path);
}
exit 0;

#--------------------------Subroutines---------------------------------#

sub match_tirs {
    my $self = shift;
    my $input_path = shift;
    my $round = shift;
    
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

                #initialize variables
                my $match_len = 0;
                my $start_pos = '';
                my $end_pos = '';
                my $match_query = '';
                my $match_hit = '';
                my $match_mis_aln = 0;
                my $total_mis_aln = 0;
                my $last_good = 0;
                my $hit_pos;
                my $query_pos;
                
                #parse homology string, keeping track of match length and mismatches or gaps
                for (my $count = 0; $count < length($homo_string); $count++) {
                    my $homo_char =  substr($homo_string, $count, 1);
                    my $query_char =  substr($query_str, $count, 1);
                    my $hit_char =  substr($hit_str, $count, 1);
                    if ($count == 6 and $total_mis_aln >= 3) {
                        $match_len = 0;
                                $start_pos = '';
                                $match_query = '';
                                $match_hit = '';
                                $end_pos = '';
                                last;
                    }
                    if ($round == 2) {
                        if ($count == 4 and $total_mis_aln >=2) {
                            $match_len = 0;
                            $start_pos = '';
                            $match_query = '';
                            $match_hit = '';
                            $end_pos = '';
                            last;
                        }
                    }
                    if ($match_len == 0){
                        #if match length equals 0 and position is not a match, continue to next position
                        if ($homo_char eq " ") {
                            $total_mis_aln++;
                            next;
                        }
                        #if position is a match, store info, continue to next position
                        elsif ($homo_char eq ":") {
                            $start_pos = $count;
                            $last_good = $count;
                            $match_len++;
                            $match_query .= $query_char;
                            $match_hit .= $hit_char;
                            next;
                        }
                    }
                    elsif ($match_len >= 1 and $match_len < 4) {
                        #if match length is 1-3 and position is not a match, increment mismatch counter and check if more than one mismatch has occurred
                        if ($homo_char eq " ") {
                            $match_mis_aln++;
                            $total_mis_aln++;
                            #allow one mismatch, store info and continue
                            if ($match_mis_aln <= 1) {
                                $match_len++;
                                $last_good = $count;
                                $match_query .= $query_char;
                                $match_hit .= $hit_char;
                                next;
                            }
                            #more than one mismatch, reset counters and other info, continue
                            elsif ($match_mis_aln > 1){ 
                                $match_len = 0;
                                $start_pos = '';
                                $match_query = '';
                                $match_hit = '';
                                $match_mis_aln = 0;
                                next;
                            }
                            elsif ($total_mis_aln >= 3) {
                                $match_len = 0;
                                $start_pos = '';
                                $match_query = '';
                                $match_hit = '';
                                $end_pos = '';
                                last;
                            }
                        }
                        #position is a match, store info and continue
                        elsif ($homo_char eq ":") {
                            $last_good = $count;
                            $match_len++;
                            $match_query .= $query_char;
                            $match_hit .= $hit_char;
                            next;
                        }
                    }
                    elsif ($match_len >= 4) {
                        #match length is 5 or higher. If position is not a match, increment mismatch counter and check if more than 2 mismatch has occurred. If a match, continue.
                        if ($homo_char eq " ") {
                            $match_mis_aln++;
                            $total_mis_aln++;
                            #mismatches under 3, store info and continue
                            if ($match_mis_aln < 3) {
                                $match_len++;
                                $last_good = $count;
                                $match_query .= $query_char;
                                $match_hit .= $hit_char;
                                next;
                            }
                            elsif($total_mis_aln >= 4) {
                                $match_len = 0;
                                $start_pos = '';
                                $match_query = '';
                                $match_hit = '';
                                $end_pos = '';
                                last;
                            }
                            #mismatches 3 or more, store final info for alignment match and end parsing
                            elsif($match_mis_aln >= 3){
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
                            next;
                        }
                    }
                }
                #add in check to see if TIRs were found. 
                if ($end_pos eq '') {
                    push @result, (0, [0]);
                }
            }
        }
    }
    return (@result);
}

sub generate_gff {
    my $self = shift; #whatever state the msa object is in at time of call
    my $path = shift;
    my $round = shift;
    print "Entered .gff printer\n";
    $self->throw("Need Bio::Align::AlignI argument")
        unless ref $self && $self->isa( 'Bio::Align::AlignI');
    
    if ($round eq 'final') {
        open(my $out, ">", $path) or die "Error creating $path. $!\n";
        print "Round = final\n";
        my $ele_info_ref = shift;
        foreach my $seq_obj ($self->each_seq()) {
            my $seq_name = $seq_obj->id();
            my $seq = $seq_obj->seq();
            $seq =~ s/-//g;
            my $seq_len = length($seq);
            my $ori_end = $seq_len - $flank;
            my $left_comp = $ele_info_ref->{$seq_name}{"left_tir_start"} - 101;
            my $right_comp = $ele_info_ref->{$seq_name}{"right_tir_start"} - $ori_end;
            my $copy_num;
            my $eleid;
            my $seqid;
            my $type = 'terminal_inverted_repeat_element'; 
            my $start;
            my $end;
            my $strand;
            if ($seq_name =~ /^([0-9]*).+[ |_]Query:(.*)[ |_]Sbjct:(.*)[ |_]Length.+Location:\(([0-9]*)[ |_]\-[ |_]([0-9]*)\)[ |_]Direction:(.+)/) {
                $copy_num = $1;
                $eleid = $2;
                $seqid = $3;
                $start = $4+$left_comp;
                $end = $5+$right_comp;
                $strand = $6;
            }
            elsif ($seq_name =~ /^([0-9]*)_(.+)_(.+)-(.+)_(.+)/) {
                $copy_num = $1;
                $seqid = $2;
                $start = $3+$left_comp;
                $end = $4+$right_comp;
                $eleid = $5;
                $strand = "?";
            }
            else {
                print "Header doesn't match TARGeT or RSPB:  $seq_name\n";
                next;
            }
            if ($eleid =~ m/(.+)_TSD/ or $eleid =~ m/(.+)_Unknow/) {
                $eleid = $1;
            }
            my $ltir_end = $start+length($ele_info_ref->{"left_tir_seq"})-1;
            my $rtir_start = $end-(length($ele_info_ref->{"right_tir_seq"})-1);
            my $ele_id = $ele_info_ref->{"element_id"};
            my $ele_class = $ele_info_ref->{"classification"};
            my $tir_id = $ele_info_ref->{"left_tir_id"};
            my $tsd_frac = $ele_info_ref->{"TSD_fraction"};
            my $tsd_con = $ele_info_ref->{"TSD_seq"};
            print $out "$seqid\tactiveTE\t$type\t$start\t$end\t.\t$strand\t.\tID=$eleid-$copy_num;Name=$eleid Copy$copy_num;element_id=$ele_id;element_classification=$ele_class;tir_id=$tir_id;tsd_fraction=$tsd_frac;tsd_consensus=$tsd_con\n";
            print $out "$seqid\tactiveTE\tfive_prime_terminal_inverted_repeat\t$start\t$ltir_end\t.\t.\t.\tParent=$eleid-$copy_num\n";
            print $out "$seqid\tactiveTE\tthree_prime_terminal_inverted_repeat\t$rtir_start\t$end\t.\t.\t.\tParent=$eleid-$copy_num\n";
            if (defined $ele_info_ref->{$seq_name}{"TSD"}) {
                my $ltsd_start = $start-length($ele_info_ref->{$seq_name}{"TSD"});
                my $ltsd_end = $start-1;
                my $rtsd_start = $end+1;
                my $rtsd_end = $end+length($ele_info_ref->{$seq_name}{"TSD"});
                print $out "$seqid\tactiveTE\ttarget_site_duplication\t$ltsd_start\t$ltsd_end\t.\t.\t.\tDerives_from=$eleid-$copy_num\n";
                print $out "$seqid\tactiveTE\ttarget_site_duplication\t$rtsd_start\t$rtsd_end\t.\t.\t.\tDerives_from=$eleid-$copy_num\n";
            }
        }
        close($out);
    }
}

sub clean_files {
    my $out_path = shift;
    opendir(my $in_DIR, $out_path) or die "Cannot open directory: $!";
    while (my $file = readdir($in_DIR)) {
        next if ($file =~ m/^\./);
        if ($file =~ m/\.(final|info|fa|element_info|bad|gff|tif|jpg)$/) {
            next ;
        }
        else{
            my $remove_path = File::Spec->catpath($volume, $out_path, $file);
            system("rm $remove_path");
        }
    }
    closedir($in_DIR);
}
