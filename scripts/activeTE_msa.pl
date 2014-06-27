#!/usr/bin/env perl

=head1 NAME

activeTE_msa - evaluate MSAs for characteristics of active TEs

=head1 SYNOPSIS

activeTE_msa family.aln

=head1 DESCRIPTION

=head1 Command line arguments

 -flank     [integer] Length of flanking sequence (default 100)
 -all       Process all
 -protein   input MSA is from protein searches

=head1 AUTHORS - Brad Cavinder

Brad Cavinder - bradc-AT-ucr.edu

Jason Stajich - jason.stajich-AT-ucr.edu [contributor]

=cut

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
use Data::Dumper;
use Cwd;
use FileHandle;

autoflush STDOUT 1;

#set default flank length and reset if option is given
my $flank = 100;
my $trimal = 'trimal';   # specify the full path to this if it can't be found automatically
my $all;
my $protein;
my $PROGRAM_NAME = "activeTE";

GetOptions(
  'flank:i'    => \$flank,
  'all'        => \$all,
  'protein'    => \$protein,
  'trimal:s'   => \$trimal,
);

# make sure only one valid input file is specified
if (@ARGV != 1 and @ARGV != 2) {
  print "Please specify one input file and, optionally, one protein match file\n";
  help();
}

my $infile = shift @ARGV;
my $protein_match = shift @ARGV;

if (!-f $infile) {
  print "Supplied MSA filepath is not valid: $!\n";
  help();
}

if ($protein_match and !-f $protein_match) {
  print "Supplied protein match filepath is not valid: $!\n";
  help();
}

if ($flank <= 25) {
  print "Please use alignments with flanks >25 nucleotides long\n";
  help();
}

sub help {
  print "

usage:

activeTE_msa.pl -a -f -p <int> <multiple alignment file>

-a keep all intermediate files  DEFAULT = remove
-f length of sequence flanks  DEFAULT = 100
-p MSA is from protein search


The test alignments all have f = 100.

";
  exit 1;
}

my $lflank_len;
my $rflank_len;
if (defined $protein) {
    if ($infile =~ /_trimmed-([0-9]+)/) {
        $lflank_len = $1;
        $rflank_len = $1;
    }
    elsif ($infile =~ /_retrimmed_([0-9]+)-([0-9]+)/) {
        $lflank_len = $1;
        $rflank_len = 2;
    }
}

## cleaning up MSA
my @in_array;
open(my $in, "<", $infile);
my $seq='';
while (my $line = <$in>) {
  if ($line =~ /^>/){
   ## replacing white space with underscores in MSA
   $line =~ s/ /_/g;
   $line =~ s/\//-/g;
   if ($seq){
     ## if there is a seq already stored store it
     my ($header,$justseq) = split /\n/ , $seq;
     if ($justseq !~ /^-+$/){
       push @in_array, $seq;
     }
     $seq = '';
   }
   $seq.=$line;
  }else{
    ##not header, just seq
    chomp $line;
    $seq.=$line;
  }
}
## for last seq
my ($header,$justseq) = split /\n/ , $seq;
if ($justseq !~ /^-+$/){
  push @in_array, $seq;
}
close($in);
open(my $fix_out, ">", $infile);
print $fix_out join ("\n",@in_array),"\n";
close($fix_out);
@in_array = ();

#breakdown directory and filename, create output directory
my ($volume, $in_dir, $filename) = File::Spec->splitpath($infile);

my @fname_fin =  split('\\.', $filename);    # first bit after first '.' is the name
pop @fname_fin;
my $fname_fin    = join('.', @fname_fin);
my $out_dir_name = "aTE_" . $fname_fin;
my $out_path     = File::Spec->catdir($in_dir, $out_dir_name);
my $p_type;

if (!-d $out_path) {
  mkdir($out_path) or die "Can't create dir:$out_path $!\n";
}
else {
    my $glob_path = File::Spec->catpath($volume, $out_path, "*");
    print "Glob path: $glob_path\n";
    my @files = glob $glob_path;
    foreach my $path (@files) {
        unlink $path or warn "Failed to unlink $path: $!";
    }
}
my $log_path = File::Spec->catpath($volume, $out_path, $filename . ".log");
open(my $log_out, '>', $log_path) or die "Can't ropen log file: $!";


if (defined $all) {
  print "Intermediate files will be kept\n";
  print $log_out "Intermediate files will be kept\n";
}
else {
  print "Intermediate files will be removed\n";
  print $log_out "Intermediate files will be removed\n";
}
if (defined $protein) {
    print "Protein tirs flagged\n";
    print $log_out "Protein tirs flagged\n";
    my @fname_type =  split('_', $fname_fin);
    $p_type = pop @fname_type;
    if ($p_type eq 'mariner') {
        $p_type = "Tc_mariner"; 
    }
    elsif ($p_type eq 'harbinger') {
        $p_type = "PIF_harbinger";
    }
}
else {
    print "Protein tirs not flagged\n";
    print $log_out "Protein tirs not flagged\n";
}

my %element_info;    # this hash is for storing element information
my %tir_positions;

#generate trimal alignment gap summary file, using sed to remove the header info
open(my $trimal_run => "$trimal -in $infile -sgc | sed -n '4,\$p' |") || die "Cannot run trimal!";
print "Running Trimal\n";
print $log_out "Running Trimal\n";

# parse the trimal output file with the gap %  of each position in the full alignment
my @gap_id_array;
while (my $line = <$trimal_run>) {
  chomp $line;       # probably not necessary
  if ($line =~ /\S+\s+(\S+)\s+\S+/) {
    push @gap_id_array, $1;
  }
}
close($trimal_run);
warn "Trimal run finished\n";

my $full_aln_obj = get_org_aln($infile);
my $full_aln_len = $full_aln_obj->length();
my $gap_cols = $full_aln_obj->gap_col_matrix();
my $full_aln_num_seqs = $full_aln_obj->num_sequences;
print "Full alignment has $full_aln_num_seqs sequences\n";
print $log_out "Full alignment has $full_aln_num_seqs sequences\n";

if ($full_aln_num_seqs <= 1 and !defined $protein){
  my $abort_out_path = File::Spec->catpath($volume, $out_path, $filename . ".abort");
  error_out($abort_out_path, "$filename\thas only $full_aln_num_seqs sequence(s)\n", 0);
}

#calculate the % nucleotide identity and the fraction of copies with sequence at each position the full alignment and print to file
my $full_id_out_path =  File::Spec->catpath($volume, $out_path, $filename . ".full_id");
my @full_id_array;
my $first_col_80 = 1;
my $last_col_80  = $full_aln_len;

print "Starting Full ID calculation\n";
warn "Starting Full ID calculation - warn\n";
print $log_out "Starting Full ID calculation\n";
@full_id_array = get_percentID_perCol($infile, $full_id_out_path);
my $try = 0;

Cleaning_MSA:
$try++;
if ($try == 2) {
    $full_aln_obj = get_org_aln($infile);
}
elsif ($try > 2) {
    #open a file to store info on why analysis of an element was aborted
    my $abort_out_path = File::Spec->catpath($volume, $out_path, $filename . ".abort");
    error_out($abort_out_path, "$filename\tTIRs not found in MSA. Multiple tries performed.", 0);
}
## remove as many seqs as possible for the cleanest aln
my ($left_tir_start1, $right_tir_start1, $tmp_aln_obj, $ref2tp, $ref2gsr, $ref2gspr) = remove_most($full_aln_obj, \%tir_positions, \@full_id_array, $try);
%tir_positions = %$ref2tp;
my %gap_seq_remove;
my @gap_seq_pos_remove;
if ($ref2gspr and $ref2gsr) {
    %gap_seq_remove     = %$ref2gsr;
    @gap_seq_pos_remove = @$ref2gspr;
}
my $remove_most = 1;
my $current_num_seq = $tmp_aln_obj->num_sequences;
print "Current number of sequences in temp_aln_obj = $current_num_seq.\n";
print $log_out "Current number of sequences in temp_aln_obj = $current_num_seq.\n";

## if this removes too many, remove as few as possible
if ($left_tir_start1 == 0 or $right_tir_start1 == 0 or $current_num_seq <= 5) {
  my $aln_obj = get_org_aln($infile);
  print "run less stringent intial filtering\n";
  warn "run less stringent intial filtering - warn\n";
  print $log_out "run less stringent intial filtering\n";
  ($left_tir_start1, $right_tir_start1, $tmp_aln_obj, $ref2tp, $ref2gsr, $ref2gspr) = remove_least($aln_obj, \%tir_positions, \@full_id_array, $try);
  %tir_positions      = %$ref2tp;
  if (defined $ref2gsr and defined $ref2gspr) {
    %gap_seq_remove = %$ref2gsr;
    @gap_seq_pos_remove = @$ref2gspr;
  }
  $remove_most = 0;
}
##overwrite full_aln with tmp_aln
$full_aln_obj = $tmp_aln_obj;

## printing out the new cleaned up alignments 
my $trim_aln_out = File::Spec->catpath($volume, $out_path, $filename . ".trim");
my $out = Bio::AlignIO->new(
  -file             => ">$trim_aln_out",
  -format           => 'fasta',
  -displayname_flat => 0
);
$out->write_aln($full_aln_obj);

## removes any columns that are now all gaps
$full_aln_obj = $full_aln_obj->remove_gaps('-', 1);
($left_tir_start1, $right_tir_start1) = get_columns($full_aln_obj, \%tir_positions, 1);
my $trimmed_aln_obj = $full_aln_obj;

if (!defined $left_tir_start1 or !defined $right_tir_start1 or $left_tir_start1 == 0 or $right_tir_start1 == 0 or $left_tir_start1 >= $right_tir_start1) {
    if ($try == 1) {
        goto Cleaning_MSA;
    }
    #open a file to store info on why analysis of an element was aborted
    my $abort_out_path = File::Spec->catpath($volume, $out_path, $filename . ".abort");
    error_out($abort_out_path, "$filename\tTIRs not found in MSA by initial search", 0);
}

my ($org_right_tir_start1, $org_left_tir_start1) = ($right_tir_start1, $left_tir_start1);
($left_tir_start1, $right_tir_start1) = adjust_tir_starts($trimmed_aln_obj, $left_tir_start1, $right_tir_start1);
my $left_tir_adjusted = $org_left_tir_start1 - $left_tir_start1;
my $right_tir_adjusted = $right_tir_start1 - $org_right_tir_start1;

#check whether any sequences have gaps within 25bp of start of putative TIRSs and remove them 
my ($left_tir_start, $right_tir_start, $ref2array, $ref2hash);
($left_tir_start, $right_tir_start, $ref2array, $trimmed_aln_obj, $ref2hash) = consensus_filter(\@gap_seq_pos_remove, $trimmed_aln_obj, $left_tir_start1, $right_tir_start1, \%tir_positions, "other", $try);
@gap_seq_pos_remove = @$ref2array;
%tir_positions      = %$ref2hash;
print "new tir starts after consensus filter: $left_tir_start, $right_tir_start\n";

print $log_out "new tir starts after consensus filter: $left_tir_start, $right_tir_start\n";
$current_num_seq = $trimmed_aln_obj->num_sequences;

if ($current_num_seq > 9) {
    ($ref2array, $trimmed_aln_obj) = gap_filter(\@gap_seq_pos_remove, $trimmed_aln_obj, $left_tir_start, $right_tir_start);
    if ($ref2array) {
        @gap_seq_pos_remove = @$ref2array;
    }
}
$trimmed_aln_obj = $trimmed_aln_obj->remove_gaps('-', 1);

print "before get_Col: Trim2 Left TIR start column: $left_tir_start1\n";
print $log_out "before get_Col: Trim2 Left TIR start column: $left_tir_start1\n";


($left_tir_start1, $right_tir_start1) =  get_columns($trimmed_aln_obj, \%tir_positions, 1);

my $test_len = $trimmed_aln_obj->length();
if ($test_len == 0) {
    if ($try ==1) {
        goto Cleaning_MSA;
    }
    #open a file to store info on why analysis of an element was aborted
    my $abort_out_path = File::Spec->catpath($volume, $out_path, $filename . ".abort");
    error_out($abort_out_path, "$filename\tTIRs not found in MSA. All sequences filtered out.", 0);
}

my $trim_aln_out2 = File::Spec->catpath($volume, $out_path, $filename . ".trim2");
my $out2 = Bio::AlignIO->new(-file => ">$trim_aln_out2", -format => 'fasta', -displayname_flat => 0);
$out2->write_aln($trimmed_aln_obj);
print "after get_Col: Trim2 Left TIR start column: $left_tir_start1\n";
print $log_out "after get_Col: Trim2 Left TIR start column: $left_tir_start1\n";
print "Trim2 Right TIR start column: $right_tir_start1 (filename.trim2)\n";
print $log_out "Trim2 Right TIR start column: $right_tir_start1 (filename.trim2)\n";

#check for too many mismatches within 25bp of start of putative TIRSs against the consensus sequence
my @consensus_remove = tir_mismatch_filter($trimmed_aln_obj, $left_tir_start1, 0, $right_tir_start1, 0, 0);

#open file to keep track of why copies were removed
my $removed_out_path = File::Spec->catpath($volume, $out_path, $filename . ".removed_sequences");
open(my $removed_out, ">", $removed_out_path);
foreach my $seq_name (@consensus_remove) {
  my $seq_obj = $trimmed_aln_obj->get_seq_by_id($seq_name);
  if (exists $tir_positions{$seq_name}) {
    delete $tir_positions{$seq_name};
  }
  if (defined $seq_obj) {
    print $removed_out "$seq_name\tToo many mismatches in TIRs\n";
    $trimmed_aln_obj->remove_seq($seq_obj);
  }
}
$trimmed_aln_obj = $trimmed_aln_obj->remove_gaps('-', 1);
($left_tir_start1, $right_tir_start1) =  get_columns($trimmed_aln_obj, \%tir_positions, 1);
print "Last TIRs before .trim3: Left tir start1: $left_tir_start1  Right tir start1: $right_tir_start1\n";

print $log_out "Last TIRs before .trim3: Left tir start1: $left_tir_start1  Right tir start1: $right_tir_start1\n";

$test_len = $trimmed_aln_obj->length();
if ($test_len == 0) {
    if ($try == 1) {
        goto Cleaning_MSA;
    }
    #open a file to store info on why analysis of an element was aborted
    my $abort_out_path = File::Spec->catpath($volume, $out_path, $filename . ".abort");
    error_out($abort_out_path, "$filename\tTIRs not found in MSA. All sequences filtered out.", 0);
}

#write new alignment to file
my $trim_aln_out3 = File::Spec->catpath($volume, $out_path, $filename . ".trim3");
$out = Bio::AlignIO->new(-file => ">$trim_aln_out3", -format => 'fasta' , -displayname_flat => 0);
$out->write_aln($trimmed_aln_obj);


#Store the column positions of the potential TIRs in the trimmed alignment and get the columns of them in the original alignment to remove sequences that cause gaps in the tirs
my %trim_left_pos_hash;
my %trim_right_pos_hash;
my %left_tir_start_check_counts;
my %right_tir_start_check_counts;

foreach my $seq_obj ($trimmed_aln_obj->each_seq()) {
  my $left_res_pos_obj  = $seq_obj->location_from_column($left_tir_start1);
  my $right_res_pos_obj = $seq_obj->location_from_column($right_tir_start1);
  my $seq               = $seq_obj->seq();
  my $seq_len           = length($seq);
  my $seq_name          = $seq_obj->id();
  my $left_res_pos      = $left_res_pos_obj->start();
  my $left_res          = substr($seq, $left_res_pos, 1);
  my $right_res_pos     = $right_res_pos_obj->start();
  my $right_res         = substr($seq, $right_res_pos, 1);

  #my $right_flank_pos   = $seq_len - $right_res_pos;
  $trim_left_pos_hash{$seq_name}  = $left_res_pos;
  $trim_right_pos_hash{$seq_name} = $right_res_pos;
  $left_tir_start_check_counts{$left_res_pos}++;
  $right_tir_start_check_counts{$right_res_pos}++;
}

#sort the left and right TIR starts by largest count to smallest
my @sorted_left_tir_start_keys = sort { $left_tir_start_check_counts{$b} <=> $left_tir_start_check_counts{$a}} keys(%left_tir_start_check_counts);
my @sorted_right_tir_start_keys = sort {$right_tir_start_check_counts{$b} <=> $right_tir_start_check_counts{$a}} keys(%right_tir_start_check_counts);

#check if TIRs start in flanks, indicating the flanks are highly similar
my $left_flank_catch  = 0;
my $right_flank_catch = 0;

print "sorted_left_tir_start_keys[0] = $sorted_left_tir_start_keys[0] ...\nsorted_right_tir_start_keys[0] = $sorted_right_tir_start_keys[0] ...\n";

print $log_out "sorted_left_tir_start_keys[0] = $sorted_left_tir_start_keys[0] ...\nsorted_right_tir_start_keys[0] = $sorted_right_tir_start_keys[0] ...\n";


if (!defined $sorted_left_tir_start_keys[0]) {
    if ($try == 1) {
        goto Cleaning_MSA;
    }
    if (!defined $sorted_right_tir_start_keys[0]) {
    #open a file to store info on why analysis of an element was aborted
    my $abort_out_path = File::Spec->catpath($volume, $out_path, $filename . ".abort");
    error_out($abort_out_path, "$filename\tTIRs not found for both flanks.", 0);
    }
    $left_flank_catch++;
    #open a file to store info on why analysis of an element was aborted
    my $abort_out_path = File::Spec->catpath($volume, $out_path, $filename . ".abort");
    error_out($abort_out_path, "$filename\tLeft flank TIR not found.", 0);
}
elsif (!defined $sorted_right_tir_start_keys[0]) {
    if ($try == 1) {
        goto Cleaning_MSA;
    }
    #open a file to store info on why analysis of an element was aborted
    my $abort_out_path = File::Spec->catpath($volume, $out_path, $filename . ".abort");
    error_out($abort_out_path, "$filename\tRight flank TIR not found.", 0);
}
  print "After removing copies with problems in TIRs these are the most common nucleotide positions of the TIRs:\n\tLeft TIR  start: $sorted_left_tir_start_keys[0]\tRight TIR start: $sorted_right_tir_start_keys[0]\n";
  print $log_out "After removing copies with problems in TIRs these are the most common nucleotide positions of the TIRs:\n\tLeft TIR  start: $sorted_left_tir_start_keys[0]\tRight TIR start: $sorted_right_tir_start_keys[0]\n";
  print "This is after the .trim2 to .trim3 transition\n";
  print $log_out "This is after the .trim2 to .trim3 transition\n";

if ($sorted_left_tir_start_keys[0] <= $flank - 25) {
  $left_flank_catch++;
}
if ($sorted_right_tir_start_keys[0] <= $flank - 25) {
  $right_flank_catch++;
}

if ($left_flank_catch != 0) {
    if ($try == 1) {
        goto Cleaning_MSA;
    }
    if ($right_flank_catch != 0) {
    #open a file to store info on why analysis of an element was aborted
    my $abort_out_path = File::Spec->catpath($volume, $out_path, $filename . ".conserved_flanks");
    error_out($abort_out_path, "$filename\tBoth flanks similar", 0);
    }
    else {
    #open a file to store info on why analysis of an element was aborted
    my $abort_out_path = File::Spec->catpath($volume, $out_path, $filename . ".conserved_left_flank");
    error_out($abort_out_path, "$filename\tLeft flank similar", 0);
    }
}
elsif ($right_flank_catch != 0) {
    if ($try == 1) {
        goto Cleaning_MSA;
    }
    #open a file to store info on why analysis of an element was aborted
    my $abort_out_path = File::Spec->catpath($volume, $out_path, $filename . ".conserved_right_flank");
    error_out($abort_out_path, "$filename\tRight flank similar", 0);
}

TIR_Finding:

my @good_aln;
my @bad_remove;
my @bad_aln;

my $element_len;
my $ele_half_len;
$test_len = $trimmed_aln_obj->length();
my $orig_cwd = cwd;

#grab all sequences from trimmed alignment
foreach my $seq_obj ($trimmed_aln_obj->each_seq()) {

  #grab sequence name and shorten it for use in some output filenames
  my $seq_name    = $seq_obj->id();
  print "Seq name first tir search: $seq_name\n";
  print $log_out "Seq name first tir search: $seq_name\n";
  my $fname_start = $seq_name;
  if ($seq_name =~ /\:/) {
    my @seqn_part   = split(":", $seq_name);
    $fname_start = $seqn_part[0];
    if ($fname_start =~ m/\|/) {
        my @seqn_part2   = split(/\|/, $fname_start);
        $fname_start = $seqn_part2[0];
    }
  }
  else {
    if ($fname_start =~ m/\|/) {
        my @seqn_part2   = split(/\|/, $fname_start);
        $fname_start = $seqn_part2[0];
    }
  }
  print "fname_start first tir search: $fname_start\n";
  print $log_out "fname_start first tir search: $fname_start\n";
  
  #grab sequence & strip hyphens
  my $seq = $seq_obj->seq();
  $seq =~ s/-//g;

  #initialize variables
  my $sub_seq;
  my $first;
  my $last;
  my $first_path = File::Spec->catpath($volume, $out_path, "first_half.txt");
  my $last_path  = File::Spec->catpath($volume, $out_path, "last_half.txt");

  #get end  sequences
  $element_len = $trim_right_pos_hash{$seq_name} - $trim_left_pos_hash{$seq_name}+1;
  $ele_half_len = int($element_len / 2);

  $first = substr($seq, $trim_left_pos_hash{$seq_name} - 1, $ele_half_len);
  $last = substr($seq, $trim_right_pos_hash{$seq_name} - $ele_half_len, $ele_half_len);

  #save the two ends as files to use as inputs for a ggsearch search
  open(my $first_out, ">", $first_path) or die "Error creating $first_path. $!\n";
  open(my $last_out, ">", $last_path) or die "Error creating $last_path. $!\n";
  print $first_out ">first\n$first\n";
  print $last_out ">last\n$last\n";
  close($first_out);
  close($last_out);
  
  #change directory for ggsearch, adjust input filepaths, create output filepath, call ggsearch, change dir back to ori
  chdir $out_path;
  $first_path = "first_half.txt";
  $last_path = "last_half.txt";
  my $out_opt = $fname_start . ".ggsearch.out";
  
  system("ggsearch36 -n -i -T 1 -d 1 $last_path $first_path > $out_opt");
  chdir $orig_cwd;
  $out_opt = File::Spec->catpath($volume, $out_path, $fname_start . ".ggsearch.out");

  my @tir_match_result = match_tirs($seq_obj, $out_opt, 1);

  if ($tir_match_result[0] == 1) {
    push @good_aln, $tir_match_result[1];
  }
  else {
    push @bad_aln,    $seq_name;
    push @bad_remove, $seq_obj;
  }
}

#find the TIRs in the trimmed alignment using the index positions of each TIR stored above
#initialize variables
my %hit_column_counts;
my %hit_match_len;
my %query_column_counts;
my %query_match_len;
my @sorted_hitcolumn_keys;
my @sorted_querycolumn_keys;
my @sorted_hit_len_keys;
my @sorted_query_len_keys;
my @entry;
my $good_aln_len = @good_aln;

if ($good_aln_len == 0) {
  print "No TIRs found from MSA: First round\n";
  print $log_out "No TIRs found from MSA: First round\n";
}
else {
  print "Found TIRs in MSA: 1st round\n";
  print $log_out "Found TIRs in MSA: 1st round\n";

  #go through each entry of @good_aln array
  foreach my $row_ref (@good_aln) {
    my $query_aln_pos;
    my $hit_aln_pos;
    @entry = @{$row_ref};

    #get column position in alignment of left (hit) TIR, using sequence name and index position of TIR in full sequence. Increment the count of that column position in %hit_column_counts. Also, increment the count of TIR length in %hit_match_len
    $hit_aln_pos = $trimmed_aln_obj->column_from_residue_number($entry[0], ${ $entry[1]{"hit"} }[0]);
    $hit_column_counts{$hit_aln_pos}++;
    $hit_match_len{${ $entry[1]{"hit"}}[1]}++;

    #do the same for the query TIR
    $query_aln_pos = $trimmed_aln_obj->column_from_residue_number($entry[0], ${ entry [1]{"query"} }[0]);
    $query_column_counts{$query_aln_pos}++;
    $query_match_len{ ${ $entry[1]{"query"} }[1] }++;

    #Storing column info for $entry[0] in first ggsearch round\nLeft TIR start: \${entry[1]{'hit'}}[0]\tRight TIR start: \${entry[1]{'query'}}[0]\n";
  }

  #sort the hit and query column and match length hashes by largest count to smallest
  @sorted_hitcolumn_keys = sort { $hit_column_counts{$b} <=> $hit_column_counts{$a} } keys(%hit_column_counts);
  @sorted_querycolumn_keys = sort { $query_column_counts{$b} <=> $query_column_counts{$a} } keys(%query_column_counts);
  @sorted_hit_len_keys = sort { $hit_match_len{$b} <=> $hit_match_len{$a} } keys(%hit_match_len);
  @sorted_query_len_keys = sort { $query_match_len{$b} <=> $query_match_len{$a} } keys(%query_match_len);

  print "\@sorted_hitcolumn_keys: $sorted_hitcolumn_keys[0]  \@sorted_hit_len_keys: $sorted_hit_len_keys[0]\n\@sorted_querycolumn_keys: $sorted_querycolumn_keys[0]  \@sorted_hit_len_keys: $sorted_hit_len_keys[0]\n";
  print $log_out "\@sorted_hitcolumn_keys: $sorted_hitcolumn_keys[0]  \@sorted_hit_len_keys: $sorted_hit_len_keys[0]\n\@sorted_querycolumn_keys: $sorted_querycolumn_keys[0]  \@sorted_hit_len_keys: $sorted_hit_len_keys[0]\n";
  
}

#repeat above but shifted out 15bp on each end to catch cases where found TIR is slightly off
my $adjustment = 15;
my %trim_left_pos_hash2;
my %trim_right_pos_hash2;
my @good_aln2;

foreach my $seq_obj ($trimmed_aln_obj->each_seq()) {
  my $left_res_pos_obj  = $seq_obj->location_from_column($left_tir_start1);
  my $right_res_pos_obj = $seq_obj->location_from_column($right_tir_start1);
  my $seq               = $seq_obj->seq();
  my $seq_len           = length($seq);
  my $left_res_pos      = ($left_res_pos_obj->start()) - $adjustment;
  $left_res_pos = $left_res_pos >= 0 ? $left_res_pos : 0;
  my $left_res = substr($seq, $left_res_pos, 1);
  my $right_res_pos = ($right_res_pos_obj->start()) + $adjustment;
  $right_res_pos = $right_res_pos < $seq_len ? $right_res_pos : $seq_len;
  my $right_res       = substr($seq, $right_res_pos, 1);
  my $seq_name        = $seq_obj->id();
  my $right_flank_pos = $seq_len - $right_res_pos;
  $trim_left_pos_hash2{$seq_name}  = $left_res_pos;
  $trim_right_pos_hash2{$seq_name} = $right_res_pos;
}

my @bad_remove2;
my @bad_aln2;

#grab all sequences from trimmed alignment
foreach my $seq_obj ($trimmed_aln_obj->each_seq()) {

  #grab sequence name and shorten it for use in some output filenames
  my $seq_name    = $seq_obj->id();
  print "Seq name 2nd tir search: $seq_name\n";
  print $log_out "Seq name 2nd tir search: $seq_name\n";
  my $fname_start = $seq_name;
  if ($seq_name =~ /\:/) {
    my @seqn_part   = split(":", $seq_name);
    $fname_start = $seqn_part[0];
    if ($fname_start =~ m/\|/) {
        my @seqn_part2   = split(/\|/, $fname_start);
        $fname_start = $seqn_part2[0];
    }
  }
  else {
    if ($fname_start =~ m/\|/) {
        my @seqn_part2   = split(/\|/, $fname_start);
        $fname_start = $seqn_part2[0];
    }
  }
  print "fname_start 2nd tir search: $fname_start\n";
  print $log_out "fname_start 2nd tir search: $fname_start\n";

  #grab sequence, make a copy, strip leading and trailing hyphens
  my $seq = $seq_obj->seq();
  $seq =~ s/-//g;

  #initialize variables
  my $sub_seq;
  my $first;
  my $last;
  my $first_path = File::Spec->catpath($volume, $out_path, "first_half.txt");
  my $last_path  = File::Spec->catpath($volume, $out_path, "last_half.txt");

  #get end  sequences
  $first = substr($seq, $trim_left_pos_hash2{$seq_name}, $ele_half_len);
  $last = substr($seq, ($trim_right_pos_hash2{$seq_name} - $ele_half_len), $ele_half_len);

  #save the two ends as files to use as inputs for a ggsearch search
  open(my $first_out, ">", $first_path) or die "Error creating $first_path. $!\n";
  open(my $last_out, ">", $last_path) or die "Error creating $last_path. $!\n";
  print $first_out ">first\n$first\n";
  print $last_out ">last\n$last\n";
  close($first_out);
  close($last_out);

  #change directory for ggsearch, adjust input filepaths, create output filepath, call ggsearch, change dir back to ori
  chdir $out_path;
  $first_path = "first_half.txt";
  $last_path = "last_half.txt";
  my $out_opt = $fname_start . ".ggsearch.out2";
  
  system("ggsearch36 -n -i -T 1 -d 1 $last_path $first_path > $out_opt");
  chdir $orig_cwd;
  $out_opt = File::Spec->catpath($volume, $out_path, $fname_start . ".ggsearch.out2");
  my @tir_match_result = match_tirs($seq_obj, $out_opt, 2);

  if ($tir_match_result[0] == 1) {
    push @good_aln2, $tir_match_result[1];
  }
  else {
    push @bad_aln2,    $seq_name;
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
my @sorted_hitcolumn_keys2;
my @sorted_querycolumn_keys2;
my @sorted_hit_len_keys2;
my @sorted_query_len_keys2;
my @entry2;

my $good_aln2_len = @good_aln2;

if ($good_aln2_len == 0) {
  print "No TIRs found from MSA: 2nd round\n";
  print $log_out "No TIRs found from MSA: 2nd round\n";
}
else {
  print "Found TIRs in MSA: 2nd round\n";
  print $log_out "Found TIRs in MSA: 2nd round\n";

  #go through each entry of @good_aln2 array
  foreach my $row_ref (@good_aln2) {
    @entry2 = @{$row_ref};

    #get column positions in of TIRs, using sequence name and index position of TIR in full sequence.
    $hit_aln_pos2 = $trimmed_aln_obj->column_from_residue_number($entry2[0], ${$entry2[1]{"hit"}}[0]);
    $hit_column_counts2{$hit_aln_pos2}++;
    $hit_match_len2{ ${ $entry2[1]{"hit"} }[1] }++;

    #do the same for the query (left) TIR
    $query_aln_pos2 = $trimmed_aln_obj->column_from_residue_number($entry2[0], ${$entry2[1]{"query"}}[0]);
    $query_column_counts2{$query_aln_pos2}++;
    $query_match_len2{ ${ $entry2[1]{"query"} }[1] }++;
  }

  #sort the hit and query column and match length hashes by largest count to smallest
  @sorted_hitcolumn_keys2 = sort {$hit_column_counts2{$b} <=> $hit_column_counts2{$a}} keys(%hit_column_counts2);
  @sorted_querycolumn_keys2 = sort {$query_column_counts2{$b} <=> $query_column_counts2{$a}} keys(%query_column_counts2);
  @sorted_hit_len_keys2 = sort {$hit_match_len2{$b} <=> $hit_match_len2{$a}} keys(%hit_match_len2);
  @sorted_query_len_keys2 = sort {$query_match_len2{$b} <=> $query_match_len2{$a}} keys(%query_match_len2);
}

my $good_aln_len2 = @good_aln2;
my $bad_aln_len   = @bad_aln;
my $bad_aln_len2  = @bad_aln2;

print "Good aln len:  $good_aln_len  Good aln len2: $good_aln_len2\n";
print $log_out "Good aln len:  $good_aln_len  Good aln len2: $good_aln_len2\n";
print "Bad aln len:  $bad_aln_len  Bad aln len2: $bad_aln_len2\n";
print $log_out "Bad aln len:  $bad_aln_len  Bad aln len2: $bad_aln_len2\n";

#check which run generated the longer TIRs
my $protein_catch = 0;
my $found = 0;
if ($good_aln_len2 == 0 and $good_aln_len == 0) {
    undef @good_aln;
    #grab all sequences from trimmed alignment
    foreach my $seq_obj ($trimmed_aln_obj->each_seq()) {
    
      #grab sequence name and shorten it for use in some output filenames
      my $seq_name    = $seq_obj->id();
      print "Seq name 3rd tir search: $seq_name\n";
      print $log_out "Seq name 3rd tir search: $seq_name\n";
      my $fname_start = $seq_name;
      if ($seq_name =~ /\:/) {
        my @seqn_part   = split(":", $seq_name);
        $fname_start = $seqn_part[0];
        if ($fname_start =~ m/\|/) {
            my @seqn_part2   = split(/\|/, $fname_start);
            $fname_start = $seqn_part2[0];
        }
      }
      else {
        if ($fname_start =~ m/\|/) {
            my @seqn_part2   = split(/\|/, $fname_start);
            $fname_start = $seqn_part2[0];
        }
      }
      print "fname_start 3rd tir search: $fname_start\n";
      print $log_out "fname_start 3rd tir search: $fname_start\n";
      
      my $out_opt = File::Spec->catpath($volume, $out_path, $fname_start . ".ggsearch.out");
    
      my @tir_match_result = match_tirs2($seq_obj, $out_opt);
    
      if ($tir_match_result[0] == 1) {
        push @good_aln, $tir_match_result[1];
      }
      else {
        push @bad_aln,    $seq_name;
        push @bad_remove, $seq_obj;
      }
    }
    $good_aln_len = @good_aln;
    
    #find the TIRs in the trimmed alignment using the index positions of each TIR stored above
    #initialize variables
    undef %hit_column_counts;
    undef %hit_match_len;
    undef %query_column_counts;
    undef %query_match_len;
    undef @sorted_hitcolumn_keys;
    undef @sorted_querycolumn_keys;
    undef @sorted_hit_len_keys;
    undef @sorted_query_len_keys;
    undef @entry;
    
    if ($good_aln_len == 0) {
      print "No TIRs found from MSA: Third round\n";
      print $log_out "No TIRs found from MSA: Third round\n";
      my $abort_out_path = File::Spec->catpath($volume, $out_path, $filename . ".abort");
        error_out($abort_out_path, "$filename\tTIRs not found by three TIR searches", 0);
      
    }
    else {
      print "Found TIRs in MSA: 3rd round\n";
      print $log_out "Found TIRs in MSA: 3rd round\n";
    
      #go through each entry of @good_aln array
      foreach my $row_ref (@good_aln) {
        my $query_aln_pos;
        my $hit_aln_pos;
        @entry = @{$row_ref};
    
        #get column position in alignment of left (hit) TIR, using sequence name and index position of TIR in full sequence. Increment the count of that column position in %hit_column_counts. Also, increment the count of TIR length in %hit_match_len
        $hit_aln_pos = $trimmed_aln_obj->column_from_residue_number($entry[0], ${ $entry[1]{"hit"} }[0]);
        $hit_column_counts{$hit_aln_pos}++;
        $hit_match_len{${ $entry[1]{"hit"}}[1]}++;
    
        #do the same for the query TIR
        $query_aln_pos = $trimmed_aln_obj->column_from_residue_number($entry[0], ${ entry [1]{"query"} }[0]);
        $query_column_counts{$query_aln_pos}++;
        $query_match_len{ ${ $entry[1]{"query"} }[1] }++;
    
        #Storing column info for $entry[0] in first ggsearch round\nLeft TIR start: \${entry[1]{'hit'}}[0]\tRight TIR start: \${entry[1]{'query'}}[0]\n";
      }
    
      #sort the hit and query column and match length hashes by largest count to smallest
      @sorted_hitcolumn_keys = sort { $hit_column_counts{$b} <=> $hit_column_counts{$a} } keys(%hit_column_counts);
      @sorted_querycolumn_keys = sort { $query_column_counts{$b} <=> $query_column_counts{$a} } keys(%query_column_counts);
      @sorted_hit_len_keys = sort { $hit_match_len{$b} <=> $hit_match_len{$a} } keys(%hit_match_len);
      @sorted_query_len_keys = sort { $query_match_len{$b} <=> $query_match_len{$a} } keys(%query_match_len);
    
      print "\@sorted_hitcolumn_keys: $sorted_hitcolumn_keys[0]  \@sorted_hit_len_keys: $sorted_hit_len_keys[0]\n\@sorted_querycolumn_keys: $sorted_querycolumn_keys[0]  \@sorted_hit_len_keys: $sorted_hit_len_keys[0]\n";
      print $log_out "\@sorted_hitcolumn_keys: $sorted_hitcolumn_keys[0]  \@sorted_hit_len_keys: $sorted_hit_len_keys[0]\n\@sorted_querycolumn_keys: $sorted_querycolumn_keys[0]  \@sorted_hit_len_keys: $sorted_hit_len_keys[0]\n";
      
    }
}
elsif ($good_aln_len2 == 0 and $good_aln_len != 0) {
  print "Alignment set 1 better\n\n";
  print $log_out "Alignment set 1 better\n\n";
  
}
elsif ($good_aln_len2 != 0 and $good_aln_len == 0) {
  print "Alignment set 2 better\n\n";
  print $log_out "Alignment set 2 better\n\n";
  @sorted_hitcolumn_keys   = @sorted_hitcolumn_keys2;
  @sorted_querycolumn_keys = @sorted_querycolumn_keys2;
  @sorted_hit_len_keys     = @sorted_hit_len_keys2;
  @sorted_query_len_keys   = @sorted_query_len_keys2;
  @bad_remove              = @bad_remove2;
  
}
elsif ($good_aln_len2 > $good_aln_len) {
  print "Alignment set 2 better\n";
  print $log_out "Alignment set 2 better\n";
  @sorted_hitcolumn_keys   = @sorted_hitcolumn_keys2;
  @sorted_querycolumn_keys = @sorted_querycolumn_keys2;
  @sorted_hit_len_keys     = @sorted_hit_len_keys2;
  @sorted_query_len_keys   = @sorted_query_len_keys2;
  @bad_remove              = @bad_remove2;
  
}
elsif ($good_aln_len2 == $good_aln_len) {
  if (($sorted_hit_len_keys2[0] > $sorted_hit_len_keys[0] and $sorted_query_len_keys2[0] > $sorted_query_len_keys[0])
    and ($sorted_hitcolumn_keys2[0] >= $sorted_hitcolumn_keys[0] and $sorted_querycolumn_keys2[0] >= $sorted_querycolumn_keys[0])) {
    print "Alignment set 2 better\n";
    print $log_out "Alignment set 2 better\n";
    @sorted_hitcolumn_keys   = @sorted_hitcolumn_keys2;
    @sorted_querycolumn_keys = @sorted_querycolumn_keys2;
    @sorted_hit_len_keys     = @sorted_hit_len_keys2;
    @sorted_query_len_keys   = @sorted_query_len_keys2;
    @bad_remove              = @bad_remove2;
    
  }
  else {
    print "Alignment set 1 better\n";
    print $log_out "Alignment set 1 better\n";
    
  }
}
else {
  print "Alignment set 1 better\n\n";
  print $log_out "Alignment set 1 better\n\n";
  
}

#store info on why analysis of a copy was aborted
foreach my $seq_obj (@bad_remove) {
  my $seq_name = $seq_obj->id();
  my $seq_obj2 = $trimmed_aln_obj->get_seq_by_id($seq_name);
  if (defined $seq_obj2) {
    print $removed_out "$seq_name\tNo TIRs found by first two ggsearch runs\n";
    $trimmed_aln_obj->remove_seq($seq_obj2);
    if (exists $tir_positions{$seq_name}) {
        delete $tir_positions{$seq_name};
    }
  }
}

print "before get_tir_nt_positions: ($sorted_hitcolumn_keys[0],$sorted_hit_len_keys[0],$sorted_querycolumn_keys[0],$sorted_query_len_keys[0])\n";
print $log_out "before get_tir_nt_positions: ($sorted_hitcolumn_keys[0],$sorted_hit_len_keys[0],$sorted_querycolumn_keys[0],$sorted_query_len_keys[0])\n";

my ($ref2_tir_positions, $ref2remove_these) = get_tir_nt_positions($trimmed_aln_obj, \%tir_positions, $sorted_hitcolumn_keys[0], $sorted_hit_len_keys[0], $sorted_querycolumn_keys[0], $sorted_query_len_keys[0]);
%tir_positions = %$ref2_tir_positions;

foreach my $key (keys %$ref2remove_these) {
    my $seq_obj = $$ref2remove_these{$key};
    if (defined $seq_obj) {
        $trimmed_aln_obj->remove_seq($seq_obj);
        if (defined $tir_positions{$key}) {
            delete $tir_positions{$key};
        }
        my @info = ($key, 0, $seq_obj);
        push @gap_seq_pos_remove, [@info];
    }
}

Reimport:
#reread the original alignment file
my $in_obj2     = Bio::AlignIO->new(-file => $infile, -format => 'fasta');
my $ori_aln_obj = $in_obj2->next_aln();
my $ori_aln_len = $ori_aln_obj->length();
my %trimmed_track;
foreach my $seq_obj($trimmed_aln_obj->each_seq) {
    my $seq_id = $seq_obj->id;
    $trimmed_track{$seq_id} = 1;
}

## because we adjusted the left tir the tir_pos is off by the number of gaps removed. but gaps were only removed from a part of the whole.
## need to go back to org, best tir, before adjust, substract the number of nt offset by adjust and by gaps, then use that col
if ($left_tir_adjusted or $right_tir_adjusted) {
  my %col;
  my $tir;
  foreach my $seq_obj($trimmed_aln_obj->each_seq) {
    my $seq_id = $seq_obj->id;
    my $seq = $seq_obj->seq;
    my $left_nt_pos_obj = $seq_obj->location_from_column($left_tir_start1 + $left_tir_adjusted); # plus ntcount adjusted
    my $right_nt_pos_obj = $seq_obj->location_from_column($right_tir_start1 - $right_tir_adjusted); # minus ntcount adjusted
    my $left_nt_pos = $left_nt_pos_obj->start;
    my $right_nt_pos = $right_nt_pos_obj->start;
    my $left_col_num = $ori_aln_obj->column_from_residue_number($seq_id,$left_nt_pos);
    my $right_col_num = $ori_aln_obj->column_from_residue_number($seq_id,$right_nt_pos);
    $col{left}{$left_col_num}++ ;
    $col{right}{$right_col_num}++ ;
   }
  my $left_col = (sort {$col{left}{$b} <=> $col{left}{$a} }keys %{$col{left}})[0];
  my $right_col = (sort {$col{right}{$b} <=> $col{right}{$a}} keys %{$col{right}})[0];
  foreach my $seq_obj($ori_aln_obj->each_seq){
    my $seq_id = $seq_obj->id;
    my $seq = $seq_obj->seq;
    my $before_left_tir = substr($seq,0,$left_col-1);
    #my $after_left_tir = substr($seq,$left_col-1);
    my $element_len = $right_col - $left_col + 1;
    my $element = substr($seq,$left_col-1,$element_len);
    my $after_right_tir = substr($seq,$right_col);
    my $left_gaps = '';
    if ($before_left_tir =~ s/(-+)$//){
      $left_gaps = $1;
    }
    my $right_gaps = '';
    if ($after_right_tir =~ s/^(-+)//){
      $right_gaps = $1;
    }
    #print substr($seq_id,0,4),":$left_gaps...$before_left_tir...$element...$after_right_tir...$right_gaps\n";
    my $new_seq = $left_gaps.$before_left_tir.$element.$after_right_tir.$right_gaps;
    $seq_obj->seq($new_seq);
  }

}

foreach my $key (keys %tir_positions) {
    #print "tir position key: $key\n";
    if (!exists $trimmed_track{$key}) {
        print "tir position key not in trimmed alignment: $key\n";
        delete $tir_positions{$key};
    }
}
#get the column positions of tirs in the original alignment
print "before get_col: lts:$left_tir_start, rts:$right_tir_start,\n";
print $log_out "before get_col: lts:$left_tir_start, rts:$right_tir_start,\n";
my ($left_tir_end, $right_tir_end);
($left_tir_start, $right_tir_start, $left_tir_end, $right_tir_end) = get_columns($ori_aln_obj, \%tir_positions, 2);
print "after get_col: $left_tir_start, $right_tir_start,$left_tir_end, $right_tir_end \n";
print $log_out "after get_col: $left_tir_start, $right_tir_start,$left_tir_end, $right_tir_end \n";

my $reimpport_num_seqs = $ori_aln_obj->num_sequences;

undef %tir_positions;

foreach my $seq_obj ($ori_aln_obj->each_seq()) {
    my $seq_name = $seq_obj->id();
    if (defined $seq_obj->location_from_column($left_tir_start)){
        my $left_tir_start_obj = $seq_obj->location_from_column($left_tir_start);
        my $left_tir_start_pos = $left_tir_start_obj->start();
        $tir_positions{$seq_name}{'left_tir_start'}  = $left_tir_start_pos;
        if (defined $seq_obj->location_from_column($left_tir_end)){
            my $left_tir_end_obj = $seq_obj->location_from_column($left_tir_end);
            my $left_tir_end_pos = $left_tir_end_obj->start();
            $tir_positions{$seq_name}{'left_tir_end'}  = $left_tir_end_pos;
        }
        else {
            delete  $tir_positions{$seq_name};
            next;
        }
    }
    else{
        next;
    }
    if (defined $seq_obj->location_from_column($right_tir_start)){
        my $right_tir_start_obj = $seq_obj->location_from_column($right_tir_start);
        my $right_tir_start_pos = $right_tir_start_obj->start();
        $tir_positions{$seq_name}{'right_tir_start'}  = $right_tir_start_pos;
        if (defined $seq_obj->location_from_column($right_tir_end)){
            my $right_tir_end_obj = $seq_obj->location_from_column($right_tir_end);
            my $right_tir_end_pos = $right_tir_end_obj->start();
            $tir_positions{$seq_name}{'right_tir_end'}  = $right_tir_end_pos;
        }
        else {
            delete  $tir_positions{$seq_name};
            next;
        }
    }
    else {
        delete  $tir_positions{$seq_name};
        next;
    }
}

($left_tir_start, $right_tir_start, $left_tir_end, $right_tir_end) = get_columns($ori_aln_obj, \%tir_positions, 2);

print "after tir_positions rebuild: $left_tir_start, $right_tir_start,$left_tir_end, $right_tir_end \n";
print $log_out "after tir_positions rebuild: $left_tir_start, $right_tir_start,$left_tir_end, $right_tir_end \n";

my $tir_length = $left_tir_end - $left_tir_start;

if ($reimpport_num_seqs > 4) {
    #go back through the full alignment and remove sequences that create gaps in the TIRs
    my %gap_seq_remove2;
    my %search_tirs;
    print "Num Seqs before first filtering after reimporting ori: ",$ori_aln_obj->num_sequences(),"\n";
    print $log_out "Num Seqs before first filtering after reimporting ori: ",$ori_aln_obj->num_sequences(),"\n";
    
    #print "remove_most is $remove_most\n";
    #print $log_out "remove_most is $remove_most\n";
    my $gap_cutoff = 50;
    my %recheck;
    
    for (my $i = 0 ; $i < $ori_aln_len ; $i++) {
      my $id_row_ref      = $full_id_array[$i];
      my @id_info         = @{$id_row_ref};
      my @gap_col_array   = @{$gap_cols};
      my $gap_col_hashref = $gap_col_array[$i];
      my %gap_col_hash    = %{$gap_col_hashref};
      my $base_count      = 0;
      my $total_count     = 0;
      if (($i >= $left_tir_start and $i <= $left_tir_end) or ($i >= $right_tir_end and $i <= $right_tir_start)) {
        foreach my $key (keys %gap_col_hash) {
          #print "\ngap col key: $key\n";
          #print $log_out "\ngap col key: $key\n";
          if ($gap_col_hash{$key} != 1) {
            $base_count++;
          }
          $total_count++;
        }
        my $pos_present = $base_count / $total_count;
        #print "i: $i  gap_id_array value: $gap_id_array[$i]  id_info[1]: $id_info[1]  pos_present: $pos_present\n";
        #print $log_out "i: $i  gap_id_array value: $gap_id_array[$i]  id_info[1]: $id_info[1]  pos_present: $pos_present\n";
        if ($gap_id_array[$i] >= $gap_cutoff and ($id_info[1] < $gap_cutoff or $pos_present < .5)) {
          foreach my $key (keys %gap_col_hash) {
            if ($gap_col_hash{$key} != 1) {
                print "Is this bad?  $key\n";
                print $log_out "Is this bad?  $key\n";
                $gap_seq_remove2{$key}++;
            }
          }
        }
      }
      elsif (($i >= $left_tir_start-10 and $i < $left_tir_start) or ($i > $right_tir_start and $i <= $right_tir_start + 10)) {
        foreach my $key (keys %gap_col_hash) {
          #print "gap col key: $key\n";
          #print $log_out "gap col key: $key\n";
          if ($gap_col_hash{$key} != 1) {
            $base_count++;
          }
          $total_count++;
        }
        my $pos_present = $base_count / $total_count;
        if ($gap_id_array[$i] < $gap_cutoff) {
          foreach my $key (keys %gap_col_hash) {
            if ($gap_col_hash{$key} == 1) {
                $recheck{$key}++;
                $gap_seq_remove2{$key}++;
            }
          }
        }
      }
    }
    
    #keep track of removed sequences to print to file
    my %removed_seq_hash;
    foreach my $key (keys %gap_seq_remove2) {
      my $remove = $ori_aln_obj->get_seq_by_id($key);
      if (exists $recheck{$key}) {
          #print "says key exists!\n";
          #print $log_out "says key exists!\n";
          $recheck{$key} = $remove;
      }
      my $seq_id = $remove->id();
      print $removed_out "$seq_id\tStep2:Sequence caused or contained gaps in at least one of the TIRs\n";
      $removed_seq_hash{$seq_id}++;
      if (defined $remove) {
        $ori_aln_obj->remove_seq($remove);
        if (exists $tir_positions{$key}) {
            delete $tir_positions{$key};
        }
      }
      if (exists $recheck{$key}) {
          #print "key = $key\n";
          #print $log_out "key = $key\n";
          if (exists $gap_seq_remove2{$key}) {
              #print "Deleting from gap_seq_remove2\n";
              #print $log_out "Deleting from gap_seq_remove2\n";
            delete $gap_seq_remove2{$key};
          }
      }
    }
    #remove tir mismatch and no tir copies already found in trimmed alignment from the reimported original
    #foreach my $seq_name (@consensus_remove) {
    #    my $remove = $ori_aln_obj->get_seq_by_id($seq_name);
    #    if (defined $remove) {
    #        $ori_aln_obj->remove_seq($remove);
    #    }
    #}
    foreach my $seq_obj (@bad_remove) {
       my $seq_name = $seq_obj->id();
       my $remove = $ori_aln_obj->get_seq_by_id($seq_name);
       if (defined $remove) {
            $ori_aln_obj->remove_seq($remove);
            if (exists $tir_positions{$seq_name}) {
                delete $tir_positions{$seq_name};
            }
       }
    }
    #add the TIR mismatch sequences removed from the trimmed alignment to gap_seq_remove2 hash so subsequent steps are performed only once
    #foreach my $seq_name (@consensus_remove) {
    #    $gap_seq_remove2{$seq_name}++;
    #}
    print "Before tir mismatch filter in first filtering after reimporting\nlts: $left_tir_start  lte: $left_tir_end  rte: $right_tir_end  rts: $right_tir_start\n";
    my @consensus_remove2 = tir_mismatch_filter($ori_aln_obj, $left_tir_start, $left_tir_end, $right_tir_start, $right_tir_end, 1, \%recheck);
    foreach my $seq_name (@consensus_remove2) {
        my $remove = $ori_aln_obj->get_seq_by_id($seq_name);
        if (defined $remove) {
            $ori_aln_obj->remove_seq($remove);
            if (exists $tir_positions{$seq_name}) {
                delete $tir_positions{$seq_name};
            }
        }
    }
    print "After tir mismatch filter in first filtering after reimporting\nlts: $left_tir_start  lte: $left_tir_end  rte: $right_tir_end  rts: $right_tir_start\n\n";
    print "Num Seqs after first filtering after reimporting ori: ",$ori_aln_obj->num_sequences(),"\n";
    print $log_out "Num Seqs after first filtering after reimporting ori: ",$ori_aln_obj->num_sequences(),"\n";
}

undef @sorted_hitcolumn_keys;
undef @sorted_querycolumn_keys;
undef @sorted_hit_len_keys;
undef @sorted_query_len_keys;
undef @bad_remove;
undef @sorted_hitcolumn_keys2;
undef @sorted_querycolumn_keys2;
undef @sorted_hit_len_keys2;
undef @sorted_query_len_keys2;
undef @bad_aln;
undef @bad_aln2;
undef @good_aln;

my $final_aln_obj = $ori_aln_obj->remove_gaps('-', 1);
my $int_aln_out =
  File::Spec->catpath($volume, $out_path, $filename . ".intermediate");
$out = Bio::AlignIO->new(
  -file             => ">$int_aln_out",
  -format           => 'fasta',
  -displayname_flat => 0
);
$out->write_aln($final_aln_obj);
my $final_len = $final_aln_obj->length();
my $final_num_seqs = $final_aln_obj->num_sequences();
print "MSA length is $final_len now\n";
print $log_out "MSA length is $final_len now\n";
print "Current number of sequences: $final_num_seqs\n";
print $log_out "Current number of sequences: $final_num_seqs\n";

if ($left_tir_start == 0 or $right_tir_start == 0){
  #open a file to store info on why analysis of an element was aborted
  my $abort_out_path = File::Spec->catpath($volume, $out_path, $filename . ".abort");
  error_out($abort_out_path, "$filename\tTIRs not found in MSA after reimporting and first round filtering of original alignment", 0);
}

#get the column positions of tirs in the intermediate alignment
($left_tir_start, $right_tir_start, $left_tir_end, $right_tir_end) = get_columns($final_aln_obj, \%tir_positions, 2);
print "2nd column grab after removing some TIR disrupting copies and removing gap only columns\nLeft Tir Start: $left_tir_start  Left Tir End: $left_tir_end\nRight Tir End: $right_tir_end  Right TIR Start: $right_tir_start (filename.intermediate)\n";
print $log_out "2nd column grab after removing some TIR disrupting copies and removing gap only columns\nLeft Tir Start: $left_tir_start  Left Tir End: $left_tir_end\nRight Tir End: $right_tir_end  Right TIR Start: $right_tir_start (filename.intermediate)\n";

if ($left_tir_start == 0 or $right_tir_start == 0){
    #open a file to store info on why analysis of an element was aborted
    my $abort_out_path = File::Spec->catpath($volume, $out_path, $filename . ".abort");
    error_out($abort_out_path, "$filename\tTIRs not found in MSA after reimporting and second round filtering of original alignment", 0);
}

if ($final_num_seqs > 4) {
    my $new_gap_cols = $final_aln_obj->gap_col_matrix();
    my %new_gap_seq_remove;
    
    for (my $i = 0 ; $i < $final_len ; $i++) {
      if (($i >= $left_tir_start-1  and $i <= $left_tir_end-1 ) or ($i >= $right_tir_end-1 and $i <= $right_tir_start-1)) {
        my $trim_gap_col_hashref = $new_gap_cols->[$i];
        my $base_count           = 0;
        my $total_count          = 0;
        foreach my $key (keys %{$trim_gap_col_hashref}) {
          if ($trim_gap_col_hashref->{$key} != 1) {
            $base_count++;
          }
          $total_count++;
        }
        my $pos_present = $base_count / $total_count;
        if ($pos_present < 0.5) {
          foreach my $key (keys %{$trim_gap_col_hashref}) {
            if ($trim_gap_col_hashref->{$key} != 1) {
              my $seq_obj = $final_aln_obj->get_seq_by_id($key);
              $new_gap_seq_remove{$key} = $seq_obj;
            }
          }
        }
        #else {
        #  foreach my $key (keys %{$trim_gap_col_hashref}) {
        #    if ($trim_gap_col_hashref->{$key} == 1) {
        #      my $seq_obj = $final_aln_obj->get_seq_by_id($key);
        #      $new_gap_seq_remove{$key} = $seq_obj;
        #    }
        #  }
        #}
      }
    }
    
    #if ($ori_aln_obj->num_sequences() - (keys %new_gap_seq_remove) > 10){
    foreach my $key (keys %new_gap_seq_remove) {
      my $seq_obj = $new_gap_seq_remove{$key};
      my $seq_id  = $seq_obj->id();
      if (defined $seq_obj) {
        print $removed_out "$seq_id\tSequence caused or contained gaps in at least one of the TIRs\n";
        print "$seq_id\tSequence caused or contained gaps in at least one of the TIRs\n";
        $final_aln_obj->remove_seq($seq_obj);
      }
    }
    #}
    $final_aln_obj = $final_aln_obj->remove_gaps('-', 1);
    my $check_len = $final_aln_obj->length();
    $final_num_seqs = $final_aln_obj->num_sequences();
    my $int_aln_out2 = File::Spec->catpath($volume, $out_path, $filename . ".intermediate2");
    $out = Bio::AlignIO->new(-file => ">$int_aln_out2", -format => 'fasta', -displayname_flat => 0);
    $out->write_aln($final_aln_obj);
    print "Current number of sequences - check2: $final_num_seqs\n";
    print $log_out "Current number of sequences - check2: $final_num_seqs\n";
    print "MSA length now: $check_len\n";
    print $log_out "MSA length now: $check_len\n";
    
}

#get the column positions of tirs in the intermediate alignment
($left_tir_start, $right_tir_start, $left_tir_end, $right_tir_end) = get_columns($final_aln_obj, \%tir_positions, 2);
print "3rd column grab after removing more copies with TIR isssues\nLeft Tir Start: $left_tir_start  Left Tir End: $left_tir_end\nRight Tir End: $right_tir_end  Right TIR Start: $right_tir_start ($filename.intermediate)\n";
print $log_out "3rd column grab after removing more copies with TIR isssues\nLeft Tir Start: $left_tir_start  Left Tir End: $left_tir_end\nRight Tir End: $right_tir_end  Right TIR Start: $right_tir_start ($filename.intermediate)\n";

if ($left_tir_start == 0 or $right_tir_start == 0){
    #open a file to store info on why analysis of an element was aborted
    my $abort_out_path = File::Spec->catpath($volume, $out_path, $filename . ".abort");
    error_out($abort_out_path, "$filename\tTIRs not found in MSA after reimporting and third round filtering of original alignment", 0);
}

foreach my $seq_obj ($final_aln_obj->each_seq()) {
  my $seq_name = $seq_obj->id();
  next if !defined $seq_obj;
  print "Seq name final tir search: $seq_name\n";
  print $log_out "Seq name final tir search: $seq_name\n";
  my $fname_start = $seq_name;
  if ($seq_name =~ /\:/) {
    my @seqn_part   = split(":", $seq_name);
    $fname_start = $seqn_part[0];
    if ($fname_start =~ m/\|/) {
        my @seqn_part2   = split(/\|/, $fname_start);
        $fname_start = $seqn_part2[0];
    }
  }
  else {
    if ($fname_start =~ m/\|/) {
        my @seqn_part2   = split(/\|/, $fname_start);
        $fname_start = $seqn_part2[0];
    }
  }
  print "fname_start final tir search: $fname_start\n";
  print $log_out "fname_start final tir search: $fname_start\n";
  
  my $seq         = $seq_obj->seq();
  my $seq_len     = length($seq);
  print "Error $seq_name:seq_len is 0\n" if $seq_len == 0;
  print $log_out "Error $seq_name:seq_len is 0\n" if $seq_len == 0;

  my $left_test_pos  = substr($seq, $left_tir_start - 1,  1);
  my $right_test_pos = substr($seq, $right_tir_start - 1, 1);
  if ($left_test_pos eq "-" or $right_test_pos eq "-") {
    print $removed_out "$seq_name\tContains a gap at the start of one or both TIRs\n";
    $final_aln_obj->remove_seq($seq_obj);
    next;
  }

  my $left_tir_start_obj  = $seq_obj->location_from_column($left_tir_start);
  my $left_tir_start_pos  = $left_tir_start_obj->start();
  my $right_tir_start_obj = $seq_obj->location_from_column($right_tir_start);
  my $right_tir_start_pos = $right_tir_start_obj->start();

  $seq =~ s/-//g;
  my $first;
  my $last;
  my $first_path = File::Spec->catpath($volume, $out_path, "first_half.txt");
  my $last_path  = File::Spec->catpath($volume, $out_path, "last_half.txt");

  $first = substr($seq, ($left_tir_start_pos - 1), $ele_half_len);
  $last = substr($seq, ($right_tir_start_pos - $ele_half_len), $ele_half_len);

  #continue TIR search
  open(my $first_out, ">", $first_path) or die "Error creating $first_path. $!\n";
  open(my $last_out, ">", $last_path) or die "Error creating $last_path. $!\n";
  print $first_out ">first\n$first\n";
  print $last_out ">last\n$last\n";
  close($first_out);
  close($last_out);
  
  
  #change directory for ggsearch, adjust input filepaths, create output filepath, call ggsearch, change dir back to ori
  chdir $out_path;
  $first_path = "first_half.txt";
  $last_path = "last_half.txt";
  my $out_opt = $fname_start . ".ggsearch3.out";
  
  system("ggsearch36 -n -i -T 1 -d 1 $last_path $first_path > $out_opt");
  chdir $orig_cwd;
  $out_opt = File::Spec->catpath($volume, $out_path, $fname_start . ".ggsearch3.out");
  
  my @tir_match_result = match_tirs($seq_obj, $out_opt, 3);

  if ($tir_match_result[0] == 0 or !@tir_match_result) {
    my $seq_id = $seq_obj->id();
    $final_aln_obj->remove_seq($seq_obj);
    delete  $tir_positions{$seq_name} if exists $tir_positions{$seq_name};
    print $removed_out "$seq_id\tNo TIR matches found by last ggsearch run\n";
  }
  else {
      push @good_aln, $tir_match_result[1];
  }
}

undef @entry;
print "Found TIRs in MSA: Final round\n";
print $log_out "Found TIRs in MSA: Final round\n";

#go through each entry of @good_aln array
foreach my $row_ref (@good_aln) {
    my $query_aln_pos;
    my $hit_aln_pos;
    @entry = @{$row_ref};

    #get column position in alignment of left (hit) TIR, using sequence name and index position of TIR in full sequence. Increment the count of that column position in %hit_column_counts. Also, increment the count of TIR length in %hit_match_len
    $hit_aln_pos = $final_aln_obj->column_from_residue_number($entry[0], ${ $entry[1]{"hit"} }[0]);
    $hit_column_counts{$hit_aln_pos}++;
    $hit_match_len{${ $entry[1]{"hit"}}[1]}++;

    #do the same for the query TIR
    $query_aln_pos = $final_aln_obj->column_from_residue_number($entry[0], ${ entry [1]{"query"} }[0]);
    $query_column_counts{$query_aln_pos}++;
    $query_match_len{ ${ $entry[1]{"query"} }[1] }++;

    #Storing column info for $entry[0] in first ggsearch round\nLeft TIR start: \${entry[1]{'hit'}}[0]\tRight TIR start: \${entry[1]{'query'}}[0]\n";
}

#sort the hit and query column and match length hashes by largest count to smallest
@sorted_hitcolumn_keys = sort { $hit_column_counts{$b} <=> $hit_column_counts{$a} } keys(%hit_column_counts);
@sorted_querycolumn_keys = sort { $query_column_counts{$b} <=> $query_column_counts{$a} } keys(%query_column_counts);
@sorted_hit_len_keys = sort { $hit_match_len{$b} <=> $hit_match_len{$a} } keys(%hit_match_len);
@sorted_query_len_keys = sort { $query_match_len{$b} <=> $query_match_len{$a} } keys(%query_match_len);

print
"\@sorted_hitcolumn_keys: $sorted_hitcolumn_keys[0]  \@sorted_hit_len_keys: $sorted_hit_len_keys[0]\n\@sorted_querycolumn_keys: $sorted_querycolumn_keys[0]  \@sorted_query_len_keys: $sorted_query_len_keys[0]\n";

print $log_out
"\@sorted_hitcolumn_keys: $sorted_hitcolumn_keys[0]  \@sorted_hit_len_keys: $sorted_hit_len_keys[0]\n\@sorted_querycolumn_keys: $sorted_querycolumn_keys[0]  \@sorted_query_len_keys: $sorted_query_len_keys[0]\n";


$left_tir_start = $sorted_hitcolumn_keys[0];
$left_tir_end = $sorted_hitcolumn_keys[0] + $sorted_hit_len_keys[0] - 1;
$right_tir_end = $sorted_querycolumn_keys[0] - $sorted_hit_len_keys[0] -1;
$right_tir_start = $sorted_querycolumn_keys[0];


if ($left_tir_start == 0 or $right_tir_start == 0 or $left_tir_end <= 0 or $right_tir_end == 0){
  #open a file to store info on why analysis of an element was aborted
  my $abort_out_path = File::Spec->catpath($volume, $out_path, $filename . ".abort");
  error_out($abort_out_path, "$filename\tTIRs not found in Final MSA", 0);
}
print "Before get columns\nlt_start: $left_tir_start  lt_end: $left_tir_end  rt_end: $right_tir_end  rt_start: $right_tir_start\n";
print $log_out "Before get columns\nlt_start: $left_tir_start  lt_end: $left_tir_end  rt_end: $right_tir_end  rt_start: $right_tir_start\n";


$final_aln_obj = $final_aln_obj->remove_gaps('-', 1);
my $int_aln_out3 = File::Spec->catpath($volume, $out_path, $filename . ".intermediate3");
$out = Bio::AlignIO->new(-file => ">$int_aln_out3", -format => 'fasta', -displayname_flat => 0);
$out->write_aln($final_aln_obj);
my $last_len = $final_aln_obj->length();

#get the column positions of tirs in the final alignment
($left_tir_start, $right_tir_start, $left_tir_end, $right_tir_end) = get_columns($final_aln_obj, \%tir_positions, 2);
print "After - lt_start: $left_tir_start  lt_end: $left_tir_end  rt_end: $right_tir_end  rt_start: $right_tir_start\n";
print $log_out "After - lt_start: $left_tir_start  lt_end: $left_tir_end  rt_end: $right_tir_end  rt_start: $right_tir_start\n";

#check for too many mismatches in putative TIRSs against the consensus sequence
my @consensus_remove3 = tir_mismatch_filter($final_aln_obj, $left_tir_start, $left_tir_end, $right_tir_start, $right_tir_end, 2);

#print to file why copies were removed
foreach my $seq_name (@consensus_remove3) {
  my $seq_obj = $final_aln_obj->get_seq_by_id($seq_name);
  print $removed_out "$seq_name\tToo many mismatches in TIRs\n";
  if (defined $seq_obj) {
    $final_aln_obj->remove_seq($seq_obj);
  }
}
$final_aln_obj = $final_aln_obj->remove_gaps('-', 1);
my $int_aln_out4 = File::Spec->catpath($volume, $out_path, $filename . ".intermediate4");
$out = Bio::AlignIO->new(-file => ">$int_aln_out4", -format => 'fasta', -displayname_flat => 0);
$out->write_aln($final_aln_obj);
$last_len = $final_aln_obj->length();

#get the column positions of tirs in the final alignment
($left_tir_start, $right_tir_start, $left_tir_end, $right_tir_end) = get_columns($final_aln_obj, \%tir_positions, 2);

print "4th column grab after removing copies with mismatches in TIRs\nLeft Tir Start: $left_tir_start  Left Tir End: $left_tir_end\nRight Tir End: $right_tir_end  Right TIR Start: $right_tir_start ($filename.intermediate3)\n";
print $log_out "4th column grab after removing copies with mismatches in TIRs\nLeft Tir Start: $left_tir_start  Left Tir End: $left_tir_end\nRight Tir End: $right_tir_end  Right TIR Start: $right_tir_start ($filename.intermediate3)\n";

if ($left_tir_start == 0 or $right_tir_start == 0 or $left_tir_end ==0 or $right_tir_end==0){
  #open a file to store info on why analysis of an element was aborted
  my $abort_out_path = File::Spec->catpath($volume, $out_path, $filename . ".abort");
  error_out($abort_out_path, "$filename\tTIRs not found after removing sequences with mismatches in TIRs", 0);
}

TSD_Finding:

#Extract the left and right TIRs as new alignment objects
my $left_TIR_aln_obj = $final_aln_obj->slice($left_tir_start, $left_tir_end, 1);
my $right_TIR_aln_obj = $final_aln_obj->slice($right_tir_end, $right_tir_start, 1);
my $element_aln_obj = $final_aln_obj->slice($left_tir_start, $right_tir_start, 1);

#find TSDs based off column # of TIRs. Must check for 2-bp and 4-bp TSD included in TIRS and then a longer TSD. Longer TSDs that start and finish with same 2- or 4-bp as the 2- or 4-bp TSDs supersede the shorter TSDs. Also, grab flanks and save them to file.
my @putative_TSD;
my @TSD_info;
my @put_TSD_names;
my @putative_TSD1;
my @TSD_info1;
my @put_TSD_names1;
my @putative_TSD2;
my @TSD_info2;
my @put_TSD_names2;
my @putative_TSD3;
my @TSD_info3;
my @put_TSD_names3;
my @putative_TSD4;
my @TSD_info4;
my @put_TSD_names4;
my @putative_TSD5;
my @TSD_info5;
my @put_TSD_names5;
my @putative_TSD6;
my @TSD_info6;
my @put_TSD_names6;

my @no_TSD_found;

#my $flanks_out_path = File::Spec->catpath($volume, $out_path, $filename . "_flanks.fa");
#open(my $flanks_out, ">", $flanks_out_path) or die "Error creating $flanks_out_path. $!\n";
my $tsd1_precount = 0;
my $tsd2_precount = 0;
my $tsd3_precount = 0;
my $tsd4_precount = 0;
my $tsd5_precount = 0;
my $tsd6_precount = 0;
my $left_short_count = 0;
my $right_short_count = 0;
my $aln_count_check = $final_aln_obj->num_sequences();
my $no_TSD_found_out_path = File::Spec->catpath($volume, $out_path, $filename . ".TSD_issues.info");
open(my $no_TSD_found_out, ">", $no_TSD_found_out_path);

foreach my $seq_obj ($final_aln_obj->each_seq()) {
  my $seq_name    = $seq_obj->id();
  my @seqn_part   = split(":", $seq_name);
  my $fname_start = $seqn_part[0];
  my $seq_ori     = $seq_obj->seq();
  my $seq         = $seq_ori;
  $seq =~ s/-//g;
  my $left_tsd;
  my $left_tsd_substr;
  my $alt_left_tsd_substr;
  my $starting_left_flank;
  my $right_tsd;
  my $right_tsd_substr;
  my $alt_right_tsd_substr;
  my $starting_right_flank;
  my $last = 0;

  #print "Seq name: $seq_name\n";
  #print $log_out "Seq name: $seq_name\n";

  #get seq_pos of the start of TIRs from column # in alignment. Adjust to include the 4-bp at outside ends of TIRs and then split into the 2-, 4-, and 10-bp sequences to search for TSDs
  my $left_tsd_loc_obj = $seq_obj->location_from_column($left_tir_start);
  my $left_tsd_end_pos = $left_tsd_loc_obj->start();
  $starting_left_flank = substr($seq, 0, $left_tsd_end_pos + 4);
  $left_tsd = substr($seq, $left_tsd_end_pos - 31, 34);
  if (length $left_tsd < 34){
    my $message = "Left flank too short to look at TSDs.\n";
    $left_short_count++;
    print $no_TSD_found_out $message;
    print "$message\n";
    print $log_out "$message\n";
    next;
  }
  $left_tsd_substr = substr($left_tsd, 10, 20);
  $alt_left_tsd_substr = substr($left_tsd, 11, 20);

  my $right_tsd_loc_obj   = $seq_obj->location_from_column($right_tir_start);
  my $right_tsd_start_pos = $right_tsd_loc_obj->start();
  $starting_right_flank = substr($seq, $right_tsd_start_pos - 4);
  $right_tsd            = substr($seq, $right_tsd_start_pos - 4, 34);
  if (length $right_tsd < 34) {
    my $message = "Right flank too short to look at TSDs.";
    $right_short_count++;
    print $no_TSD_found_out $message;
    print "$message\n";
    print $log_out "$message\n";
    next;
  }
  $right_tsd_substr = substr($right_tsd, 4, 20);
  $alt_right_tsd_substr = substr($right_tsd, 3, 20);
  
  my $left_2bp  = substr($left_tsd,  -4, 2);
  my $right_2bp = substr($right_tsd, 2,  2);
  my $left_4bp  = substr($left_tsd,  -4, 4);
  my $right_4bp = substr($right_tsd, 0,  4);
  
  if ($left_2bp eq $right_2bp) {
    $tsd1_precount++;
  }
  if ($left_4bp eq $right_4bp) {
    $tsd2_precount++;
  }
  for (my $i = 0 ; $i < 19 ; $i++) {
    if (substr($left_tsd_substr, $i) eq substr($right_tsd_substr, 0, -($i))) {
      $tsd3_precount++;
      last;
    }
  }
  for (my $i = 0 ; $i < 19 ; $i++) {
    if (substr($alt_left_tsd_substr, $i) eq substr($alt_right_tsd_substr, 0, -($i))) {
      $tsd4_precount++;
      last;
    }
  }
  for (my $i = 0 ; $i < 19 ; $i++) {
    if (substr($left_tsd_substr, $i) eq substr($alt_right_tsd_substr, 0, -($i))) {
      $tsd5_precount++;
      last;
    }
  }
  for (my $i = 0 ; $i < 19 ; $i++) {
    if (substr($alt_left_tsd_substr, $i) eq substr($right_tsd_substr, 0, -($i))) {
      $tsd6_precount++;
      last;
    }
  }
}

if ($aln_count_check - $left_short_count < 1) {
    if ($aln_count_check - $right_short_count < 1) {
        my $abort_out_path = File::Spec->catpath($volume, $out_path, $filename . ".conserved_flanks");
        error_out($abort_out_path, "$filename\tBoth flanks similar, TIRs too close to ends for TSD finding.", 0);
    }
    else {
        my $abort_out_path = File::Spec->catpath($volume, $out_path, $filename . ".conserved_left_flank");
        error_out($abort_out_path, "$filename\tLeft flank similar, TIR start too close to end for TSD finding.", 0);
    }
}
elsif ($aln_count_check - $right_short_count < 1) {
    my $abort_out_path = File::Spec->catpath($volume, $out_path, $filename . ".conserved_right_flank");
    error_out($abort_out_path, "$filename\tRight flank similar, TIR start too close to end for TSD finding.", 0);
}


print "Starting code to find TSDs\n";
print $log_out "Starting code to find TSDs\n";
my $tsd1_count = 0;
my $tsd2_count = 0;
my $tsd3_count = 0;
my $tsd4_count = 0;
my $tsd5_count = 0;
my $tsd6_count = 0;

my $counter = 0;
foreach my $seq_obj ($final_aln_obj->each_seq()) {
  my $seq_name    = $seq_obj->id();
  my @seqn_part   = split(":", $seq_name);
  my $fname_start = $seqn_part[0];
  my $seq_ori     = $seq_obj->seq();
  my $seq         = $seq_ori;
  $seq =~ s/-//g;
  my $left_tsd;
  my $left_tsd_substr;
  my $alt_left_tsd_substr;
  my $starting_left_flank;
  my $right_tsd;
  my $right_tsd_substr;
  my $alt_right_tsd_substr;
  my $starting_right_flank;
  my $last = 0;

  #print "Seq name: $seq_name\n";
  #print $log_out "Seq name: $seq_name\n";

  #get seq_pos of the start of TIRs from column # in alignment. Adjust to include the 4-bp at outside ends of TIRs and then split into the 2-, 4-, and 10-bp sequences to search for TSDs
  my $left_tsd_loc_obj = $seq_obj->location_from_column($left_tir_start);
  my $left_tsd_end_pos = $left_tsd_loc_obj->start();
  $starting_left_flank = substr($seq, 0, $left_tsd_end_pos + 4);
  $left_tsd = substr($seq, $left_tsd_end_pos - 31, 34);
  if (length $left_tsd < 34){
    #my $message = "Left flank too short to look at TSDs.\n";
    #$left_short_count++;
    #print $no_TSD_found_out $message;
    #print "$message\n";
    #print $log_out "$message\n";
    next;
  }
  $left_tsd_substr = substr($left_tsd, 10, 20);
  $alt_left_tsd_substr = substr($left_tsd, 11, 20);

  my $right_tsd_loc_obj   = $seq_obj->location_from_column($right_tir_start);
  my $right_tsd_start_pos = $right_tsd_loc_obj->start();
  $starting_right_flank = substr($seq, $right_tsd_start_pos - 4);
  $right_tsd            = substr($seq, $right_tsd_start_pos - 4, 34);
  if (length $right_tsd < 34) {
    #my $message = "Right flank too short to look at TSDs.";
    #$right_short_count++;
    #print $no_TSD_found_out $message;
    #print "$message\n";
    #print $log_out "$message\n";
    next;
  }
  $right_tsd_substr = substr($right_tsd, 4, 20);
  $alt_right_tsd_substr = substr($right_tsd, 3, 20);

  print "$seq_name\nleft tsd: $left_tsd_substr  alt left tsd: $alt_left_tsd_substr\nright tsd: $right_tsd_substr  alt right tsd: $alt_right_tsd_substr\n";
  print $log_out "$seq_name\nleft tsd: $left_tsd_substr  alt left tsd: $alt_left_tsd_substr\nright tsd: $right_tsd_substr  alt right tsd: $alt_right_tsd_substr\n";

  my $left_2bp  = substr($left_tsd,  -4, 2);
  my $right_2bp = substr($right_tsd, 2,  2);
  my $left_4bp  = substr($left_tsd,  -4, 4);
  my $right_4bp = substr($right_tsd, 0,  4);
  my $tsd1;
  my $tsd2;
  my $tsd3;
  my $tsd4;
  my $tsd5;
  my $tsd6;
  
  my $tsd4_catch = 0;

  if ($left_2bp eq $right_2bp) {
    print "Round 1 TSD search - Yes! $fname_start\n";
    print $log_out "Round 1 TSD search - Yes! $fname_start\n";
    $tsd1 = $left_2bp;
  }
  if ($left_4bp eq $right_4bp) {
    print "Round 2 TSD search - Yes! $fname_start\n";
    print $log_out "Round 2 TSD search - Yes! $fname_start\n";
    $tsd2 = $left_4bp;
  }
  for (my $i = 0 ; $i < 19 ; $i++) {
    if (substr($left_tsd_substr, $i) eq substr($right_tsd_substr, 0, -($i))) {
      print "Round 3 TSD search - Yes! $i $fname_start\n";
      print $log_out "Round 3 TSD search - Yes! $i $fname_start\n";
      $tsd3 = substr($left_tsd_substr, $i);
      last;
    }
  }
  for (my $i = 0 ; $i < 19 ; $i++) {
    if (substr($alt_left_tsd_substr, $i) eq substr($alt_right_tsd_substr, 0, -($i))) {
      print "Round 4 TSD search - Yes! $i $fname_start\n";
      print $log_out "Round 4 TSD search - Yes! $i $fname_start\n";
      $tsd4 = substr($alt_left_tsd_substr, $i);
      last;
    }
  }
  for (my $i = 0 ; $i < 19 ; $i++) {
    if (substr($left_tsd_substr, $i) eq substr($alt_right_tsd_substr, 0, -($i))) {
      print "Round 5 TSD search - Yes! $i $fname_start\n";
      print $log_out "Round 5 TSD search - Yes! $i $fname_start\n";
      $tsd5 = substr($alt_left_tsd_substr, $i);
      last;
    }
  }
  for (my $i = 0 ; $i < 19 ; $i++) {
    if (substr($alt_left_tsd_substr, $i) eq substr($right_tsd_substr, 0, -($i))) {
      print "Round 6 TSD search - Yes! $i $fname_start\n";
      print $log_out "Round 6 TSD search - Yes! $i $fname_start\n";
      $tsd6 = substr($alt_left_tsd_substr, $i);
      last;
    }
  }

  #Save found TSD to an array or report that a TSD wasn't found
  if ($tsd4 or $tsd5 or $tsd6) {
        $tsd4_catch = 1;
        if ($tsd4) {
            if ($tsd5 and !$tsd6) {
                if ((length($tsd4) >= length($tsd5) or $tsd4_count > $tsd5_count) and ($tsd5_count <= $tsd4_count or $tsd5_precount < $tsd4_precount)) {
                    if (!$tsd1 and !$tsd2 and !$tsd3) {
                        my $insertion_site = substr($left_tsd, (-3-length($tsd4)-10), (length($tsd4)+10)) . substr($right_tsd, 3+length($tsd4), 10);
                        push @TSD_info4, [ $seq_name, $insertion_site, $tsd4 ];
                        push @putative_TSD4, $tsd4;
                        push @put_TSD_names4, [ $seq_name, $tsd4 ];
                        print "tsd4: $fname_start\n\n";
                        print $log_out "tsd4: $fname_start\n\n";
                        $tsd4_count++;
                        $tsd4_catch = 1;
                        undef $tsd1;
                        undef $tsd2;
                        undef $tsd3;
                        undef $tsd5;
                        undef $tsd6;
                    }
                    elsif (!$tsd1 and !$tsd2 and $tsd3) {
                        if ((length($tsd4) >= length($tsd3) or $tsd4_count > $tsd3_count) and ($tsd3_count <= $tsd4_count or $tsd3_precount < $tsd4_precount)) {
                            my $insertion_site = substr($left_tsd, (-3-length($tsd4)-10), (length($tsd4)+10)) . substr($right_tsd, 3+length($tsd4), 10);
                            push @TSD_info4, [ $seq_name, $insertion_site, $tsd4 ];
                            push @putative_TSD4, $tsd4;
                            push @put_TSD_names4, [ $seq_name, $tsd4 ];
                            print "tsd4: $fname_start\n\n";
                            print $log_out "tsd4: $fname_start\n\n";
                            $tsd4_count++;
                            $tsd4_catch = 1;
                            undef $tsd1;
                            undef $tsd2;
                            undef $tsd3;
                            undef $tsd5;
                            undef $tsd6;
                        }
                        else {
                            my $insertion_site = substr($left_tsd, (-4 - length($tsd3) - 10),(length($tsd3) + 10)) . substr($right_tsd, 4 + length($tsd3), 10);
                            push @TSD_info3, [ $seq_name, $insertion_site, $tsd3 ];
                            push @putative_TSD3, $tsd3;
                            push @put_TSD_names3, [ $seq_name, $tsd3 ];
                            print "tsd3: $fname_start\n\n";
                            print $log_out "tsd3: $fname_start\n\n";
                            $tsd3_count++;
                            undef $tsd1;
                            undef $tsd2;
                            undef $tsd4;
                            undef $tsd5;
                            undef $tsd6;
                        }
                    }
                    elsif (!$tsd1_count and !$tsd3_count and $tsd2) {
                        if ((length($tsd4) >= length($tsd2) or $tsd4_count > $tsd2_count) and ($tsd2_count <= $tsd4_count or $tsd2_precount < $tsd4_precount)) {
                            my $insertion_site = substr($left_tsd, (-3-length($tsd4)-10), (length($tsd4)+10)) . substr($right_tsd, 3+length($tsd4), 10);
                            push @TSD_info4, [ $seq_name, $insertion_site, $tsd4 ];
                            push @putative_TSD4, $tsd4;
                            push @put_TSD_names4, [ $seq_name, $tsd4 ];
                            print "tsd4: $fname_start\n\n";
                            print $log_out "tsd4: $fname_start\n\n";
                            $tsd4_count++;
                            $tsd4_catch = 1;
                            undef $tsd1;
                            undef $tsd2;
                            undef $tsd3;
                            undef $tsd5;
                            undef $tsd6;
                        }
                        else {
                            my $insertion_site = substr($left_tsd, (-14)) . substr($right_tsd, 4, 10);
                            push @TSD_info2, [ $seq_name, $insertion_site, $tsd2 ];
                            push @putative_TSD2, $tsd2;
                            push @put_TSD_names2, [ $seq_name, $tsd2 ];
                            print "tsd2: $fname_start\n\n";
                            print $log_out "tsd2: $fname_start\n\n";
                            $tsd2_count++;
                            undef $tsd1;
                            undef $tsd4;
                            undef $tsd3;
                            undef $tsd5;
                            undef $tsd6;
                        }
                    }
                    elsif (!$tsd2 and !$tsd3 and $tsd1) {
                        if ((length($tsd4) >= length($tsd1) or $tsd4_count > $tsd1_count) and ($tsd1_count <= $tsd4_count or $tsd1_precount < $tsd4_precount)) {
                            my $insertion_site = substr($left_tsd, (-3-length($tsd4)-10), (length($tsd4)+10)) . substr($right_tsd, 3+length($tsd4), 10);
                            push @TSD_info4, [ $seq_name, $insertion_site, $tsd4 ];
                            push @putative_TSD4, $tsd4;
                            push @put_TSD_names4, [ $seq_name, $tsd4 ];
                            print "tsd4: $fname_start\n\n";
                            print $log_out "tsd4: $fname_start\n\n";
                            $tsd4_count++;
                            $tsd4_catch = 1;
                            undef $tsd1;
                            undef $tsd2;
                            undef $tsd3;
                            undef $tsd5;
                            undef $tsd6;
                        }
                        else {
                            my $insertion_site = substr($left_tsd, -14, 12) . substr($right_tsd, 4, 10);
                            push @TSD_info1, [ $seq_name, $insertion_site, $tsd1 ];
                            push @putative_TSD1, $tsd1;
                            push @put_TSD_names1, [ $seq_name, $tsd1 ];
                            print "tsd1: $fname_start\n\n";
                            print $log_out "tsd1: $fname_start\n\n";
                            $tsd1_count++;
                            undef $tsd4;
                            undef $tsd2;
                            undef $tsd3;
                            undef $tsd5;
                            undef $tsd6;
                        }
                    }
                }
                else {
                    if (!$tsd1 and !$tsd2 and !$tsd3) {
                        my $insertion_site = substr($left_tsd, (-4-length($tsd5)-10), (length($tsd5)+10)) . substr($right_tsd, 3+length($tsd5), 10);
                        push @TSD_info5, [ $seq_name, $insertion_site, $tsd5 ];
                        push @putative_TSD5, $tsd5;
                        push @put_TSD_names5, [ $seq_name, $tsd5 ];
                        print "tsd5: $fname_start\n\n";
                        print $log_out "tsd5: $fname_start\n\n";
                        $tsd5_count++;
                        $tsd4_catch = 1;
                        undef $tsd1;
                        undef $tsd2;
                        undef $tsd3;
                        undef $tsd4;
                        undef $tsd6;
                    }
                    elsif (!$tsd1 and !$tsd2 and $tsd3) {
                        if ((length($tsd5) >= length($tsd3) or $tsd5_count > $tsd3_count) and ($tsd3_count <= $tsd5_count or $tsd3_precount < $tsd5_precount)) {
                            my $insertion_site = substr($left_tsd, (-4-length($tsd5)-10), (length($tsd5)+10)) . substr($right_tsd, 3+length($tsd5), 10);
                            push @TSD_info5, [ $seq_name, $insertion_site, $tsd5 ];
                            push @putative_TSD5, $tsd5;
                            push @put_TSD_names5, [ $seq_name, $tsd5 ];
                            print "tsd5: $fname_start\n\n";
                            print $log_out "tsd5: $fname_start\n\n";
                            $tsd5_count++;
                            $tsd4_catch = 1;
                            undef $tsd1;
                            undef $tsd2;
                            undef $tsd3;
                            undef $tsd4;
                            undef $tsd6;
                        }
                        else {
                            my $insertion_site = substr($left_tsd, (-4 - length($tsd3) - 10),(length($tsd3) + 10)) . substr($right_tsd, 4 + length($tsd3), 10);
                            push @TSD_info3, [ $seq_name, $insertion_site, $tsd3 ];
                            push @putative_TSD3, $tsd3;
                            push @put_TSD_names3, [ $seq_name, $tsd3 ];
                            print "tsd3: $fname_start\n\n";
                            print $log_out "tsd3: $fname_start\n\n";
                            $tsd3_count++;
                            undef $tsd1;
                            undef $tsd2;
                            undef $tsd4;
                            undef $tsd6;
                            undef $tsd5;
                        }
                    }
                    elsif (!$tsd1_count and !$tsd3_count and $tsd2) {
                        if ((length($tsd5) >= length($tsd2) or $tsd5_count > $tsd2_count) and ($tsd2_count <= $tsd5_count or $tsd2_precount < $tsd5_precount)) {
                            my $insertion_site = substr($left_tsd, (-4-length($tsd5)-10), (length($tsd5)+10)) . substr($right_tsd, 3+length($tsd5), 10);
                            push @TSD_info5, [ $seq_name, $insertion_site, $tsd5 ];
                            push @putative_TSD5, $tsd5;
                            push @put_TSD_names5, [ $seq_name, $tsd5 ];
                            print "tsd5: $fname_start\n\n";
                            print $log_out "tsd5: $fname_start\n\n";
                            $tsd5_count++;
                            $tsd4_catch = 1;
                            undef $tsd1;
                            undef $tsd2;
                            undef $tsd3;
                            undef $tsd4;
                            undef $tsd6;
                        }
                        else {
                            my $insertion_site = substr($left_tsd, (-14)) . substr($right_tsd, 4, 10);
                            push @TSD_info2, [ $seq_name, $insertion_site, $tsd2 ];
                            push @putative_TSD2, $tsd2;
                            push @put_TSD_names2, [ $seq_name, $tsd2 ];
                            print "tsd2: $fname_start\n\n";
                            print $log_out "tsd2: $fname_start\n\n";
                            $tsd2_count++;
                            undef $tsd1;
                            undef $tsd4;
                            undef $tsd3;
                            undef $tsd6;
                            undef $tsd5;
                        }
                    }
                    elsif (!$tsd2 and !$tsd3 and $tsd1) {
                        if ((length($tsd5) >= length($tsd1) or $tsd5_count > $tsd1_count) and ($tsd1_count <= $tsd5_count or $tsd1_precount < $tsd5_precount)) {
                            my $insertion_site = substr($left_tsd, (-4-length($tsd5)-10), (length($tsd5)+10)) . substr($right_tsd, 3+length($tsd5), 10);
                            push @TSD_info5, [ $seq_name, $insertion_site, $tsd5 ];
                            push @putative_TSD5, $tsd5;
                            push @put_TSD_names5, [ $seq_name, $tsd5 ];
                            print "tsd5: $fname_start\n\n";
                            print $log_out "tsd5: $fname_start\n\n";
                            $tsd5_count++;
                            $tsd4_catch = 1;
                            undef $tsd1;
                            undef $tsd2;
                            undef $tsd3;
                            undef $tsd4;
                            undef $tsd6;
                        }
                        else {
                            my $insertion_site = substr($left_tsd, -14, 12) . substr($right_tsd, 4, 10);
                            push @TSD_info1, [ $seq_name, $insertion_site, $tsd1 ];
                            push @putative_TSD1, $tsd1;
                            push @put_TSD_names1, [ $seq_name, $tsd1 ];
                            print "tsd1: $fname_start\n\n";
                            print $log_out "tsd1: $fname_start\n\n";
                            $tsd1_count++;
                            undef $tsd4;
                            undef $tsd2;
                            undef $tsd3;
                            undef $tsd6;
                            undef $tsd5;
                        }
                    }
                }
            }
            elsif ($tsd6 and !$tsd5) {
                if ((length($tsd4) >= length($tsd6) or $tsd4_count > $tsd6_count) and ($tsd6_count <= $tsd4_count or $tsd6_precount < $tsd4_precount)) {
                    if (!$tsd1 and !$tsd2 and !$tsd3) {
                        my $insertion_site = substr($left_tsd, (-3-length($tsd4)-10), (length($tsd4)+10)) . substr($right_tsd, 3+length($tsd4), 10);
                        push @TSD_info4, [ $seq_name, $insertion_site, $tsd4 ];
                        push @putative_TSD4, $tsd4;
                        push @put_TSD_names4, [ $seq_name, $tsd4 ];
                        print "tsd4: $fname_start\n\n";
                        print $log_out "tsd4: $fname_start\n\n";
                        $tsd4_count++;
                        $tsd4_catch = 1;
                        undef $tsd1;
                        undef $tsd2;
                        undef $tsd3;
                        undef $tsd5;
                        undef $tsd6;
                    }
                    elsif (!$tsd1 and !$tsd2 and $tsd3) {
                        if ((length($tsd4) >= length($tsd3) or $tsd4_count > $tsd3_count) and ($tsd3_count <= $tsd4_count or $tsd3_precount < $tsd4_precount)) {
                            my $insertion_site = substr($left_tsd, (-3-length($tsd4)-10), (length($tsd4)+10)) . substr($right_tsd, 3+length($tsd4), 10);
                            push @TSD_info4, [ $seq_name, $insertion_site, $tsd4 ];
                            push @putative_TSD4, $tsd4;
                            push @put_TSD_names4, [ $seq_name, $tsd4 ];
                            print "tsd4: $fname_start\n\n";
                            print $log_out "tsd4: $fname_start\n\n";
                            $tsd4_count++;
                            $tsd4_catch = 1;
                            undef $tsd1;
                            undef $tsd2;
                            undef $tsd3;
                            undef $tsd6;
                            undef $tsd5;
                        }
                        else {
                            my $insertion_site = substr($left_tsd, (-4 - length($tsd3) - 10),(length($tsd3) + 10)) . substr($right_tsd, 4 + length($tsd3), 10);
                            push @TSD_info3, [ $seq_name, $insertion_site, $tsd3 ];
                            push @putative_TSD3, $tsd3;
                            push @put_TSD_names3, [ $seq_name, $tsd3 ];
                            print "tsd3: $fname_start\n\n";
                            print $log_out "tsd3: $fname_start\n\n";
                            $tsd3_count++;
                            undef $tsd1;
                            undef $tsd2;
                            undef $tsd5;
                            undef $tsd4;
                            undef $tsd6;
                        }
                    }
                    elsif (!$tsd1_count and !$tsd3_count and $tsd2) {
                        if ((length($tsd4) >= length($tsd2) or $tsd4_count > $tsd2_count) and ($tsd2_count <= $tsd4_count or $tsd2_precount < $tsd4_precount)) {
                            my $insertion_site = substr($left_tsd, (-3-length($tsd4)-10), (length($tsd4)+10)) . substr($right_tsd, 3+length($tsd4), 10);
                            push @TSD_info4, [ $seq_name, $insertion_site, $tsd4 ];
                            push @putative_TSD4, $tsd4;
                            push @put_TSD_names4, [ $seq_name, $tsd4 ];
                            print "tsd4: $fname_start\n\n";
                            print $log_out "tsd4: $fname_start\n\n";
                            $tsd4_count++;
                            $tsd4_catch = 1;
                            undef $tsd1;
                            undef $tsd2;
                            undef $tsd3;
                            undef $tsd5;
                            undef $tsd6;
                        }
                        else {
                            my $insertion_site = substr($left_tsd, (-14)) . substr($right_tsd, 4, 10);
                            push @TSD_info2, [ $seq_name, $insertion_site, $tsd2 ];
                            push @putative_TSD2, $tsd2;
                            push @put_TSD_names2, [ $seq_name, $tsd2 ];
                            print "tsd2: $fname_start\n\n";
                            print $log_out "tsd2: $fname_start\n\n";
                            $tsd2_count++;
                            undef $tsd1;
                            undef $tsd4;
                            undef $tsd3;
                            undef $tsd5;
                            undef $tsd6;
                        }
                    }
                    elsif (!$tsd2 and !$tsd3 and $tsd1) {
                        if ((length($tsd4) >= length($tsd1) or $tsd4_count > $tsd1_count) and ($tsd1_count <= $tsd4_count or $tsd1_precount < $tsd4_precount)) {
                            my $insertion_site = substr($left_tsd, (-3-length($tsd4)-10), (length($tsd4)+10)) . substr($right_tsd, 3+length($tsd4), 10);
                            push @TSD_info4, [ $seq_name, $insertion_site, $tsd4 ];
                            push @putative_TSD4, $tsd4;
                            push @put_TSD_names4, [ $seq_name, $tsd4 ];
                            print "tsd4: $fname_start\n\n";
                            print $log_out "tsd4: $fname_start\n\n";
                            $tsd4_count++;
                            $tsd4_catch = 1;
                            undef $tsd1;
                            undef $tsd2;
                            undef $tsd3;
                            undef $tsd5;
                            undef $tsd6;
                        }
                        else {
                            my $insertion_site = substr($left_tsd, -14, 12) . substr($right_tsd, 4, 10);
                            push @TSD_info1, [ $seq_name, $insertion_site, $tsd1 ];
                            push @putative_TSD1, $tsd1;
                            push @put_TSD_names1, [ $seq_name, $tsd1 ];
                            print "tsd1: $fname_start\n\n";
                            print $log_out "tsd1: $fname_start\n\n";
                            $tsd1_count++;
                            undef $tsd4;
                            undef $tsd2;
                            undef $tsd3;
                            undef $tsd5;
                            undef $tsd6;
                        }
                    }
                }
                else {
                    if (!$tsd1 and !$tsd2 and !$tsd3) {
                        my $insertion_site = substr($left_tsd, (-3-length($tsd6)-10), (length($tsd6)+10)) . substr($right_tsd, 4+length($tsd6), 10);
                        push @TSD_info6, [ $seq_name, $insertion_site, $tsd6 ];
                        push @putative_TSD6, $tsd6;
                        push @put_TSD_names6, [ $seq_name, $tsd6 ];
                        print "tsd6: $fname_start\n\n";
                        print $log_out "tsd6: $fname_start\n\n";
                        $tsd6_count++;
                        $tsd4_catch = 1;
                        undef $tsd1;
                        undef $tsd2;
                        undef $tsd3;
                        undef $tsd4;
                        undef $tsd5;
                    }
                    elsif (!$tsd1 and !$tsd2 and $tsd3) {
                        if ((length($tsd6) >= length($tsd3) or $tsd6_count > $tsd3_count) and ($tsd3_count <= $tsd6_count or $tsd3_precount < $tsd6_precount)) {
                            my $insertion_site = substr($left_tsd, (-3-length($tsd6)-10), (length($tsd6)+10)) . substr($right_tsd, 4+length($tsd6), 10);
                            push @TSD_info6, [ $seq_name, $insertion_site, $tsd6 ];
                            push @putative_TSD6, $tsd6;
                            push @put_TSD_names6, [ $seq_name, $tsd6 ];
                            print "tsd6: $fname_start\n\n";
                            print $log_out "tsd6: $fname_start\n\n";
                            $tsd6_count++;
                            $tsd4_catch = 1;
                            undef $tsd1;
                            undef $tsd2;
                            undef $tsd3;
                            undef $tsd4;
                            undef $tsd5;
                        }
                        else {
                            my $insertion_site = substr($left_tsd, (-4 - length($tsd3) - 10),(length($tsd3) + 10)) . substr($right_tsd, 4 + length($tsd3), 10);
                            push @TSD_info3, [ $seq_name, $insertion_site, $tsd3 ];
                            push @putative_TSD3, $tsd3;
                            push @put_TSD_names3, [ $seq_name, $tsd3 ];
                            print "tsd3: $fname_start\n\n";
                            print $log_out "tsd3: $fname_start\n\n";
                            $tsd3_count++;
                            undef $tsd1;
                            undef $tsd2;
                            undef $tsd4;
                            undef $tsd6;
                            undef $tsd5;
                        }
                    }
                    elsif (!$tsd1_count and !$tsd3_count and $tsd2) {
                        if ((length($tsd6) >= length($tsd2) or $tsd6_count > $tsd2_count) and ($tsd2_count <= $tsd6_count or $tsd2_precount < $tsd6_precount)) {
                            my $insertion_site = substr($left_tsd, (-3-length($tsd6)-10), (length($tsd6)+10)) . substr($right_tsd, 4+length($tsd6), 10);
                            push @TSD_info6, [ $seq_name, $insertion_site, $tsd6 ];
                            push @putative_TSD6, $tsd6;
                            push @put_TSD_names6, [ $seq_name, $tsd6 ];
                            print "tsd6: $fname_start\n\n";
                            print $log_out "tsd6: $fname_start\n\n";
                            $tsd6_count++;
                            $tsd4_catch = 1;
                            undef $tsd1;
                            undef $tsd2;
                            undef $tsd3;
                            undef $tsd4;
                            undef $tsd5;
                        }
                        else {
                            my $insertion_site = substr($left_tsd, (-14)) . substr($right_tsd, 4, 10);
                            push @TSD_info2, [ $seq_name, $insertion_site, $tsd2 ];
                            push @putative_TSD2, $tsd2;
                            push @put_TSD_names2, [ $seq_name, $tsd2 ];
                            print "tsd2: $fname_start\n\n";
                            print $log_out "tsd2: $fname_start\n\n";
                            $tsd2_count++;
                            undef $tsd1;
                            undef $tsd4;
                            undef $tsd3;
                            undef $tsd6;
                            undef $tsd5;
                        }
                    }
                    elsif (!$tsd2 and !$tsd3 and $tsd1) {
                        if ((length($tsd6) >= length($tsd1) or $tsd6_count > $tsd1_count) and ($tsd1_count <= $tsd6_count or $tsd1_precount < $tsd6_precount)) {
                            my $insertion_site = substr($left_tsd, (-3-length($tsd6)-10), (length($tsd6)+10)) . substr($right_tsd, 4+length($tsd6), 10);
                            push @TSD_info6, [ $seq_name, $insertion_site, $tsd6 ];
                            push @putative_TSD6, $tsd6;
                            push @put_TSD_names6, [ $seq_name, $tsd6 ];
                            print "tsd6: $fname_start\n\n";
                            print $log_out "tsd6: $fname_start\n\n";
                            $tsd6_count++;
                            $tsd4_catch = 1;
                            undef $tsd1;
                            undef $tsd2;
                            undef $tsd3;
                            undef $tsd4;
                            undef $tsd5;
                        }
                        else {
                            my $insertion_site = substr($left_tsd, -14, 12) . substr($right_tsd, 4, 10);
                            push @TSD_info1, [ $seq_name, $insertion_site, $tsd1 ];
                            push @putative_TSD1, $tsd1;
                            push @put_TSD_names1, [ $seq_name, $tsd1 ];
                            print "tsd1: $fname_start\n\n";
                            print $log_out "tsd1: $fname_start\n\n";
                            $tsd1_count++;
                            undef $tsd4;
                            undef $tsd2;
                            undef $tsd3;
                            undef $tsd6;
                            undef $tsd5;
                        }
                    }
                }
            }
            elsif ($tsd5 and $tsd6) {
                if ((length($tsd4) >= length($tsd5) or $tsd4_count > $tsd5_count) and ($tsd5_count <= $tsd4_count or $tsd5_precount < $tsd4_precount)) {
                    if ((length($tsd4) >= length($tsd6) or $tsd4_count > $tsd6_count) and ($tsd6_count <= $tsd4_count or $tsd6_precount < $tsd4_precount)) {
                        if (!$tsd1 and !$tsd2 and !$tsd3) {
                            my $insertion_site = substr($left_tsd, (-3-length($tsd4)-10), (length($tsd4)+10)) . substr($right_tsd, 3+length($tsd4), 10);
                            push @TSD_info4, [ $seq_name, $insertion_site, $tsd4 ];
                            push @putative_TSD4, $tsd4;
                            push @put_TSD_names4, [ $seq_name, $tsd4 ];
                            print "tsd4: $fname_start\n\n";
                            print $log_out "tsd4: $fname_start\n\n";
                            $tsd4_count++;
                            $tsd4_catch = 1;
                            undef $tsd1;
                            undef $tsd2;
                            undef $tsd3;
                            undef $tsd5;
                            undef $tsd6;
                        }
                        elsif (!$tsd1 and !$tsd2 and $tsd3) {
                            if ((length($tsd4) >= length($tsd3) or $tsd4_count > $tsd3_count) and ($tsd3_count <= $tsd4_count or $tsd3_precount < $tsd4_precount)) {
                                my $insertion_site = substr($left_tsd, (-3-length($tsd4)-10), (length($tsd4)+10)) . substr($right_tsd, 3+length($tsd4), 10);
                                push @TSD_info4, [ $seq_name, $insertion_site, $tsd4 ];
                                push @putative_TSD4, $tsd4;
                                push @put_TSD_names4, [ $seq_name, $tsd4 ];
                                print "tsd4: $fname_start\n\n";
                                print $log_out "tsd4: $fname_start\n\n";
                                $tsd4_count++;
                                $tsd4_catch = 1;
                                undef $tsd1;
                                undef $tsd2;
                                undef $tsd3;
                                undef $tsd5;
                                undef $tsd6;
                            }
                            else {
                                my $insertion_site = substr($left_tsd, (-4 - length($tsd3) - 10),(length($tsd3) + 10)) . substr($right_tsd, 4 + length($tsd3), 10);
                                push @TSD_info3, [ $seq_name, $insertion_site, $tsd3 ];
                                push @putative_TSD3, $tsd3;
                                push @put_TSD_names3, [ $seq_name, $tsd3 ];
                                print "tsd3: $fname_start\n\n";
                                print $log_out "tsd3: $fname_start\n\n";
                                $tsd3_count++;
                                undef $tsd1;
                                undef $tsd2;
                                undef $tsd4;
                                undef $tsd5;
                                undef $tsd6;
                            }
                        }
                        elsif (!$tsd1_count and !$tsd3_count and $tsd2) {
                            if ((length($tsd4) >= length($tsd2) or $tsd4_count > $tsd2_count) and ($tsd2_count <= $tsd4_count or $tsd2_precount < $tsd4_precount)) {
                                my $insertion_site = substr($left_tsd, (-3-length($tsd4)-10), (length($tsd4)+10)) . substr($right_tsd, 3+length($tsd4), 10);
                                push @TSD_info4, [ $seq_name, $insertion_site, $tsd4 ];
                                push @putative_TSD4, $tsd4;
                                push @put_TSD_names4, [ $seq_name, $tsd4 ];
                                print "tsd4: $fname_start\n\n";
                                print $log_out "tsd4: $fname_start\n\n";
                                $tsd4_count++;
                                $tsd4_catch = 1;
                                undef $tsd1;
                                undef $tsd2;
                                undef $tsd3;
                                undef $tsd5;
                                undef $tsd6;
                            }
                            else {
                                my $insertion_site = substr($left_tsd, (-14)) . substr($right_tsd, 4, 10);
                                push @TSD_info2, [ $seq_name, $insertion_site, $tsd2 ];
                                push @putative_TSD2, $tsd2;
                                push @put_TSD_names2, [ $seq_name, $tsd2 ];
                                print "tsd2: $fname_start\n\n";
                                print $log_out "tsd2: $fname_start\n\n";
                                $tsd2_count++;
                                undef $tsd1;
                                undef $tsd4;
                                undef $tsd3;
                                undef $tsd5;
                                undef $tsd6;
                            }
                        }
                        elsif (!$tsd2 and !$tsd3 and $tsd1) {
                            if ((length($tsd4) >= length($tsd1) or $tsd4_count > $tsd1_count) and ($tsd1_count <= $tsd4_count or $tsd1_precount < $tsd4_precount)) {
                                my $insertion_site = substr($left_tsd, (-3-length($tsd4)-10), (length($tsd4)+10)) . substr($right_tsd, 3+length($tsd4), 10);
                                push @TSD_info4, [ $seq_name, $insertion_site, $tsd4 ];
                                push @putative_TSD4, $tsd4;
                                push @put_TSD_names4, [ $seq_name, $tsd4 ];
                                print "tsd4: $fname_start\n\n";
                                print $log_out "tsd4: $fname_start\n\n";
                                $tsd4_count++;
                                $tsd4_catch = 1;
                                undef $tsd1;
                                undef $tsd2;
                                undef $tsd3;
                                undef $tsd5;
                                undef $tsd6;
                            }
                            else {
                                my $insertion_site = substr($left_tsd, -14, 12) . substr($right_tsd, 4, 10);
                                push @TSD_info1, [ $seq_name, $insertion_site, $tsd1 ];
                                push @putative_TSD1, $tsd1;
                                push @put_TSD_names1, [ $seq_name, $tsd1 ];
                                print "tsd1: $fname_start\n\n";
                                print $log_out "tsd1: $fname_start\n\n";
                                $tsd1_count++;
                                undef $tsd4;
                                undef $tsd2;
                                undef $tsd3;
                                undef $tsd5;
                                undef $tsd6;
                            }
                        }
                    }
                    else {
                        if (!$tsd1 and !$tsd2 and !$tsd3) {
                            my $insertion_site = substr($left_tsd, (-3-length($tsd6)-10), (length($tsd6)+10)) . substr($right_tsd, 4+length($tsd6), 10);
                            push @TSD_info6, [ $seq_name, $insertion_site, $tsd6 ];
                            push @putative_TSD6, $tsd6;
                            push @put_TSD_names6, [ $seq_name, $tsd6 ];
                            print "tsd6: $fname_start\n\n";
                            print $log_out "tsd6: $fname_start\n\n";
                            $tsd6_count++;
                            $tsd4_catch = 1;
                            undef $tsd1;
                            undef $tsd2;
                            undef $tsd3;
                            undef $tsd4;
                            undef $tsd5;
                        }
                        elsif (!$tsd1 and !$tsd2 and $tsd3) {
                            if ((length($tsd6) >= length($tsd3) or $tsd6_count > $tsd3_count) and ($tsd3_count <= $tsd6_count or $tsd3_precount < $tsd6_precount)) {
                                my $insertion_site = substr($left_tsd, (-3-length($tsd6)-10), (length($tsd6)+10)) . substr($right_tsd, 4+length($tsd6), 10);
                                push @TSD_info6, [ $seq_name, $insertion_site, $tsd6 ];
                                push @putative_TSD6, $tsd6;
                                push @put_TSD_names6, [ $seq_name, $tsd6 ];
                                print "tsd6: $fname_start\n\n";
                                print $log_out "tsd6: $fname_start\n\n";
                                $tsd6_count++;
                                $tsd4_catch = 1;
                                undef $tsd1;
                                undef $tsd2;
                                undef $tsd3;
                                undef $tsd4;
                                undef $tsd5;
                            }
                            else {
                                my $insertion_site = substr($left_tsd, (-4 - length($tsd3) - 10),(length($tsd3) + 10)) . substr($right_tsd, 4 + length($tsd3), 10);
                                push @TSD_info3, [ $seq_name, $insertion_site, $tsd3 ];
                                push @putative_TSD3, $tsd3;
                                push @put_TSD_names3, [ $seq_name, $tsd3 ];
                                print "tsd3: $fname_start\n\n";
                                print $log_out "tsd3: $fname_start\n\n";
                                $tsd3_count++;
                                undef $tsd1;
                                undef $tsd2;
                                undef $tsd4;
                                undef $tsd6;
                                undef $tsd5;
                            }
                        }
                        elsif (!$tsd1_count and !$tsd3_count and $tsd2) {
                            if ((length($tsd6) >= length($tsd2) or $tsd6_count > $tsd2_count) and ($tsd2_count <= $tsd6_count or $tsd2_precount < $tsd6_precount)) {
                                my $insertion_site = substr($left_tsd, (-3-length($tsd6)-10), (length($tsd6)+10)) . substr($right_tsd, 4+length($tsd6), 10);
                                push @TSD_info6, [ $seq_name, $insertion_site, $tsd6 ];
                                push @putative_TSD6, $tsd6;
                                push @put_TSD_names6, [ $seq_name, $tsd6 ];
                                print "tsd6: $fname_start\n\n";
                                print $log_out "tsd6: $fname_start\n\n";
                                $tsd6_count++;
                                $tsd4_catch = 1;
                                undef $tsd1;
                                undef $tsd2;
                                undef $tsd3;
                                undef $tsd4;
                                undef $tsd5;
                            }
                            else {
                                my $insertion_site = substr($left_tsd, (-14)) . substr($right_tsd, 4, 10);
                                push @TSD_info2, [ $seq_name, $insertion_site, $tsd2 ];
                                push @putative_TSD2, $tsd2;
                                push @put_TSD_names2, [ $seq_name, $tsd2 ];
                                print "tsd2: $fname_start\n\n";
                                print $log_out "tsd2: $fname_start\n\n";
                                $tsd2_count++;
                                undef $tsd1;
                                undef $tsd4;
                                undef $tsd3;
                                undef $tsd6;
                                undef $tsd5;
                            }
                        }
                        elsif (!$tsd2 and !$tsd3 and $tsd1) {
                            if ((length($tsd6) >= length($tsd1) or $tsd6_count > $tsd1_count) and ($tsd1_count <= $tsd6_count or $tsd1_precount < $tsd6_precount)) {
                                my $insertion_site = substr($left_tsd, (-3-length($tsd6)-10), (length($tsd6)+10)) . substr($right_tsd, 4+length($tsd6), 10);
                                push @TSD_info6, [ $seq_name, $insertion_site, $tsd6 ];
                                push @putative_TSD6, $tsd6;
                                push @put_TSD_names6, [ $seq_name, $tsd6 ];
                                print "tsd6: $fname_start\n\n";
                                print $log_out "tsd6: $fname_start\n\n";
                                $tsd6_count++;
                                $tsd4_catch = 1;
                                undef $tsd1;
                                undef $tsd2;
                                undef $tsd3;
                                undef $tsd4;
                                undef $tsd5;
                            }
                            else {
                                my $insertion_site = substr($left_tsd, -14, 12) . substr($right_tsd, 4, 10);
                                push @TSD_info1, [ $seq_name, $insertion_site, $tsd1 ];
                                push @putative_TSD1, $tsd1;
                                push @put_TSD_names1, [ $seq_name, $tsd1 ];
                                print "tsd1: $fname_start\n\n";
                                print $log_out "tsd1: $fname_start\n\n";
                                $tsd1_count++;
                                undef $tsd4;
                                undef $tsd2;
                                undef $tsd3;
                                undef $tsd6;
                                undef $tsd5;
                            }
                        }
                    }
                }
                elsif ((length($tsd5) >= length($tsd6) or $tsd5_count > $tsd6_count) and ($tsd6_count <= $tsd5_count or $tsd6_precount < $tsd5_precount)) {
                    if (!$tsd1 and !$tsd2 and !$tsd3) {
                        my $insertion_site = substr($left_tsd, (-4-length($tsd5)-10), (length($tsd5)+10)) . substr($right_tsd, 3+length($tsd5), 10);
                        push @TSD_info5, [ $seq_name, $insertion_site, $tsd5 ];
                        push @putative_TSD5, $tsd5;
                        push @put_TSD_names5, [ $seq_name, $tsd5 ];
                        print "tsd5: $fname_start\n\n";
                        print $log_out "tsd5: $fname_start\n\n";
                        $tsd5_count++;
                        $tsd4_catch = 1;
                        undef $tsd1;
                        undef $tsd2;
                        undef $tsd3;
                        undef $tsd4;
                        undef $tsd6;
                    }
                    elsif (!$tsd1 and !$tsd2 and $tsd3) {
                        if ((length($tsd5) >= length($tsd3) or $tsd5_count > $tsd3_count) and ($tsd3_count <= $tsd5_count or $tsd3_precount < $tsd5_precount)) {
                            my $insertion_site = substr($left_tsd, (-4-length($tsd5)-10), (length($tsd5)+10)) . substr($right_tsd, 3+length($tsd5), 10);
                            push @TSD_info5, [ $seq_name, $insertion_site, $tsd5 ];
                            push @putative_TSD5, $tsd5;
                            push @put_TSD_names5, [ $seq_name, $tsd5 ];
                            print "tsd5: $fname_start\n\n";
                            print $log_out "tsd5: $fname_start\n\n";
                            $tsd5_count++;
                            $tsd4_catch = 1;
                            undef $tsd1;
                            undef $tsd2;
                            undef $tsd3;
                            undef $tsd4;
                            undef $tsd6;
                        }
                        else {
                            my $insertion_site = substr($left_tsd, (-4 - length($tsd3) - 10),(length($tsd3) + 10)) . substr($right_tsd, 4 + length($tsd3), 10);
                            push @TSD_info3, [ $seq_name, $insertion_site, $tsd3 ];
                            push @putative_TSD3, $tsd3;
                            push @put_TSD_names3, [ $seq_name, $tsd3 ];
                            print "tsd3: $fname_start\n\n";
                            print $log_out "tsd3: $fname_start\n\n";
                            $tsd3_count++;
                            undef $tsd1;
                            undef $tsd2;
                            undef $tsd4;
                            undef $tsd6;
                            undef $tsd5;
                        }
                    }
                    elsif (!$tsd1_count and !$tsd3_count and $tsd2) {
                        if ((length($tsd5) >= length($tsd2) or $tsd5_count > $tsd2_count) and ($tsd2_count <= $tsd5_count or $tsd2_precount < $tsd5_precount)) {
                            my $insertion_site = substr($left_tsd, (-4-length($tsd5)-10), (length($tsd5)+10)) . substr($right_tsd, 3+length($tsd5), 10);
                            push @TSD_info5, [ $seq_name, $insertion_site, $tsd5 ];
                            push @putative_TSD5, $tsd5;
                            push @put_TSD_names5, [ $seq_name, $tsd5 ];
                            print "tsd5: $fname_start\n\n";
                            print $log_out "tsd5: $fname_start\n\n";
                            $tsd5_count++;
                            $tsd4_catch = 1;
                            undef $tsd1;
                            undef $tsd2;
                            undef $tsd3;
                            undef $tsd4;
                            undef $tsd6;
                        }
                        else {
                            my $insertion_site = substr($left_tsd, (-14)) . substr($right_tsd, 4, 10);
                            push @TSD_info2, [ $seq_name, $insertion_site, $tsd2 ];
                            push @putative_TSD2, $tsd2;
                            push @put_TSD_names2, [ $seq_name, $tsd2 ];
                            print "tsd2: $fname_start\n\n";
                            print $log_out "tsd2: $fname_start\n\n";
                            $tsd2_count++;
                            undef $tsd1;
                            undef $tsd4;
                            undef $tsd3;
                            undef $tsd6;
                            undef $tsd5;
                        }
                    }
                    elsif (!$tsd2 and !$tsd3 and $tsd1) {
                        if ((length($tsd5) >= length($tsd1) or $tsd5_count > $tsd1_count) and ($tsd1_count <= $tsd5_count or $tsd1_precount < $tsd5_precount)) {
                            my $insertion_site = substr($left_tsd, (-4-length($tsd5)-10), (length($tsd5)+10)) . substr($right_tsd, 3+length($tsd5), 10);
                            push @TSD_info5, [ $seq_name, $insertion_site, $tsd5 ];
                            push @putative_TSD5, $tsd5;
                            push @put_TSD_names5, [ $seq_name, $tsd5 ];
                            print "tsd5: $fname_start\n\n";
                            print $log_out "tsd5: $fname_start\n\n";
                            $tsd5_count++;
                            $tsd4_catch = 1;
                            undef $tsd1;
                            undef $tsd2;
                            undef $tsd3;
                            undef $tsd4;
                            undef $tsd6;
                        }
                        else {
                            my $insertion_site = substr($left_tsd, -14, 12) . substr($right_tsd, 4, 10);
                            push @TSD_info1, [ $seq_name, $insertion_site, $tsd1 ];
                            push @putative_TSD1, $tsd1;
                            push @put_TSD_names1, [ $seq_name, $tsd1 ];
                            print "tsd1: $fname_start\n\n";
                            print $log_out "tsd1: $fname_start\n\n";
                            $tsd1_count++;
                            undef $tsd4;
                            undef $tsd2;
                            undef $tsd3;
                            undef $tsd6;
                            undef $tsd5;
                        }
                    }
                }
                else {
                    if (!$tsd1 and !$tsd2 and !$tsd3) {
                        my $insertion_site = substr($left_tsd, (-3-length($tsd6)-10), (length($tsd6)+10)) . substr($right_tsd, 4+length($tsd6), 10);
                        push @TSD_info6, [ $seq_name, $insertion_site, $tsd6 ];
                        push @putative_TSD6, $tsd6;
                        push @put_TSD_names6, [ $seq_name, $tsd6 ];
                        print "tsd6: $fname_start\n\n";
                        print $log_out "tsd6: $fname_start\n\n";
                        $tsd6_count++;
                        $tsd4_catch = 1;
                        undef $tsd1;
                        undef $tsd2;
                        undef $tsd3;
                        undef $tsd4;
                        undef $tsd5;
                    }
                    elsif (!$tsd1 and !$tsd2 and $tsd3) {
                        if ((length($tsd6) >= length($tsd3) or $tsd6_count > $tsd3_count) and ($tsd3_count <= $tsd6_count or $tsd3_precount < $tsd6_precount)) {
                            my $insertion_site = substr($left_tsd, (-3-length($tsd6)-10), (length($tsd6)+10)) . substr($right_tsd, 4+length($tsd6), 10);
                            push @TSD_info6, [ $seq_name, $insertion_site, $tsd6 ];
                            push @putative_TSD6, $tsd6;
                            push @put_TSD_names6, [ $seq_name, $tsd6 ];
                            print "tsd6: $fname_start\n\n";
                            print $log_out "tsd6: $fname_start\n\n";
                            $tsd6_count++;
                            $tsd4_catch = 1;
                            undef $tsd1;
                            undef $tsd2;
                            undef $tsd3;
                            undef $tsd4;
                            undef $tsd5;
                        }
                        else {
                            my $insertion_site = substr($left_tsd, (-4 - length($tsd3) - 10),(length($tsd3) + 10)) . substr($right_tsd, 4 + length($tsd3), 10);
                            push @TSD_info3, [ $seq_name, $insertion_site, $tsd3 ];
                            push @putative_TSD3, $tsd3;
                            push @put_TSD_names3, [ $seq_name, $tsd3 ];
                            print "tsd3: $fname_start\n\n";
                            print $log_out "tsd3: $fname_start\n\n";
                            $tsd3_count++;
                            undef $tsd1;
                            undef $tsd2;
                            undef $tsd4;
                            undef $tsd6;
                            undef $tsd5;
                        }
                    }
                    elsif (!$tsd1_count and !$tsd3_count and $tsd2) {
                        if ((length($tsd6) >= length($tsd2) or $tsd6_count > $tsd2_count) and ($tsd2_count <= $tsd6_count or $tsd2_precount < $tsd6_precount)) {
                            my $insertion_site = substr($left_tsd, (-3-length($tsd6)-10), (length($tsd6)+10)) . substr($right_tsd, 4+length($tsd6), 10);
                            push @TSD_info6, [ $seq_name, $insertion_site, $tsd6 ];
                            push @putative_TSD6, $tsd6;
                            push @put_TSD_names6, [ $seq_name, $tsd6 ];
                            print "tsd6: $fname_start\n\n";
                            print $log_out "tsd6: $fname_start\n\n";
                            $tsd6_count++;
                            $tsd4_catch = 1;
                            undef $tsd1;
                            undef $tsd2;
                            undef $tsd3;
                            undef $tsd4;
                            undef $tsd5;
                        }
                        else {
                            my $insertion_site = substr($left_tsd, (-14)) . substr($right_tsd, 4, 10);
                            push @TSD_info2, [ $seq_name, $insertion_site, $tsd2 ];
                            push @putative_TSD2, $tsd2;
                            push @put_TSD_names2, [ $seq_name, $tsd2 ];
                            print "tsd2: $fname_start\n\n";
                            print $log_out "tsd2: $fname_start\n\n";
                            $tsd2_count++;
                            undef $tsd1;
                            undef $tsd4;
                            undef $tsd3;
                            undef $tsd6;
                            undef $tsd5;
                        }
                    }
                    elsif (!$tsd2 and !$tsd3 and $tsd1) {
                        if ((length($tsd6) >= length($tsd1) or $tsd6_count > $tsd1_count) and ($tsd1_count <= $tsd6_count or $tsd1_precount < $tsd6_precount)) {
                            my $insertion_site = substr($left_tsd, (-3-length($tsd6)-10), (length($tsd6)+10)) . substr($right_tsd, 4+length($tsd6), 10);
                            push @TSD_info6, [ $seq_name, $insertion_site, $tsd6 ];
                            push @putative_TSD6, $tsd6;
                            push @put_TSD_names6, [ $seq_name, $tsd6 ];
                            print "tsd6: $fname_start\n\n";
                            print $log_out "tsd6: $fname_start\n\n";
                            $tsd6_count++;
                            $tsd4_catch = 1;
                            undef $tsd1;
                            undef $tsd2;
                            undef $tsd3;
                            undef $tsd4;
                            undef $tsd5;
                        }
                        else {
                            my $insertion_site = substr($left_tsd, -14, 12) . substr($right_tsd, 4, 10);
                            push @TSD_info1, [ $seq_name, $insertion_site, $tsd1 ];
                            push @putative_TSD1, $tsd1;
                            push @put_TSD_names1, [ $seq_name, $tsd1 ];
                            print "tsd1: $fname_start\n\n";
                            print $log_out "tsd1: $fname_start\n\n";
                            $tsd1_count++;
                            undef $tsd4;
                            undef $tsd2;
                            undef $tsd3;
                            undef $tsd6;
                            undef $tsd5;
                        }
                    }
                }
            }
            else {
                if (!$tsd1 and !$tsd2 and !$tsd3) {
                    my $insertion_site = substr($left_tsd, (-3-length($tsd4)-10), (length($tsd4)+10)) . substr($right_tsd, 3+length($tsd4), 10);
                    push @TSD_info4, [ $seq_name, $insertion_site, $tsd4 ];
                    push @putative_TSD4, $tsd4;
                    push @put_TSD_names4, [ $seq_name, $tsd4 ];
                    print "tsd4: $fname_start\n\n";
                    print $log_out "tsd4: $fname_start\n\n";
                    $tsd4_count++;
                    $tsd4_catch = 1;
                    undef $tsd1;
                    undef $tsd2;
                    undef $tsd3;
                    undef $tsd5;
                    undef $tsd6;
                }
                elsif (!$tsd1 and !$tsd2 and $tsd3) {
                    if ((length($tsd4) >= length($tsd3) or $tsd4_count > $tsd3_count) and ($tsd3_count <= $tsd4_count or $tsd3_precount < $tsd4_precount)) {
                        my $insertion_site = substr($left_tsd, (-3-length($tsd4)-10), (length($tsd4)+10)) . substr($right_tsd, 3+length($tsd4), 10);
                        push @TSD_info4, [ $seq_name, $insertion_site, $tsd4 ];
                        push @putative_TSD4, $tsd4;
                        push @put_TSD_names4, [ $seq_name, $tsd4 ];
                        print "tsd4: $fname_start\n\n";
                        print $log_out "tsd4: $fname_start\n\n";
                        $tsd4_count++;
                        $tsd4_catch = 1;
                        undef $tsd1;
                        undef $tsd2;
                        undef $tsd3;
                        undef $tsd5;
                        undef $tsd6;
                    }
                    else {
                        my $insertion_site = substr($left_tsd, (-4 - length($tsd3) - 10),(length($tsd3) + 10)) . substr($right_tsd, 4 + length($tsd3), 10);
                        push @TSD_info3, [ $seq_name, $insertion_site, $tsd3 ];
                        push @putative_TSD3, $tsd3;
                        push @put_TSD_names3, [ $seq_name, $tsd3 ];
                        print "tsd3: $fname_start\n\n";
                        print $log_out "tsd3: $fname_start\n\n";
                        $tsd3_count++;
                        undef $tsd1;
                        undef $tsd2;
                        undef $tsd4;
                        undef $tsd5;
                        undef $tsd6;
                    }
                }
                elsif (!$tsd1_count and !$tsd3_count and $tsd2) {
                    if ((length($tsd4) >= length($tsd2) or $tsd4_count > $tsd2_count) and ($tsd2_count <= $tsd4_count or $tsd2_precount < $tsd4_precount)) {
                        my $insertion_site = substr($left_tsd, (-3-length($tsd4)-10), (length($tsd4)+10)) . substr($right_tsd, 3+length($tsd4), 10);
                        push @TSD_info4, [ $seq_name, $insertion_site, $tsd4 ];
                        push @putative_TSD4, $tsd4;
                        push @put_TSD_names4, [ $seq_name, $tsd4 ];
                        print "tsd4: $fname_start\n\n";
                        print $log_out "tsd4: $fname_start\n\n";
                        $tsd4_count++;
                        $tsd4_catch = 1;
                        undef $tsd1;
                        undef $tsd2;
                        undef $tsd3;
                        undef $tsd5;
                        undef $tsd6;
                    }
                    else {
                        my $insertion_site = substr($left_tsd, (-14)) . substr($right_tsd, 4, 10);
                        push @TSD_info2, [ $seq_name, $insertion_site, $tsd2 ];
                        push @putative_TSD2, $tsd2;
                        push @put_TSD_names2, [ $seq_name, $tsd2 ];
                        print "tsd2: $fname_start\n\n";
                        print $log_out "tsd2: $fname_start\n\n";
                        $tsd2_count++;
                        undef $tsd1;
                        undef $tsd4;
                        undef $tsd3;
                        undef $tsd5;
                        undef $tsd6;
                    }
                }
                elsif (!$tsd2 and !$tsd3 and $tsd1) {
                    if ((length($tsd4) >= length($tsd1) or $tsd4_count > $tsd1_count) and ($tsd1_count <= $tsd4_count or $tsd1_precount < $tsd4_precount)) {
                        my $insertion_site = substr($left_tsd, (-3-length($tsd4)-10), (length($tsd4)+10)) . substr($right_tsd, 3+length($tsd4), 10);
                        push @TSD_info4, [ $seq_name, $insertion_site, $tsd4 ];
                        push @putative_TSD4, $tsd4;
                        push @put_TSD_names4, [ $seq_name, $tsd4 ];
                        print "tsd4: $fname_start\n\n";
                        print $log_out "tsd4: $fname_start\n\n";
                        $tsd4_count++;
                        $tsd4_catch = 1;
                        undef $tsd1;
                        undef $tsd2;
                        undef $tsd3;
                        undef $tsd5;
                        undef $tsd6;
                    }
                    else {
                        my $insertion_site = substr($left_tsd, -14, 12) . substr($right_tsd, 4, 10);
                        push @TSD_info1, [ $seq_name, $insertion_site, $tsd1 ];
                        push @putative_TSD1, $tsd1;
                        push @put_TSD_names1, [ $seq_name, $tsd1 ];
                        print "tsd1: $fname_start\n\n";
                        print $log_out "tsd1: $fname_start\n\n";
                        $tsd1_count++;
                        undef $tsd4;
                        undef $tsd2;
                        undef $tsd3;
                        undef $tsd5;
                        undef $tsd6;
                    }
                }
            }
        }
        elsif ($tsd5) {
            if (!$tsd6) {
                if (!$tsd1 and !$tsd2 and !$tsd3) {
                    my $insertion_site = substr($left_tsd, (-4-length($tsd5)-10), (length($tsd5)+10)) . substr($right_tsd, 3+length($tsd5), 10);
                    push @TSD_info5, [ $seq_name, $insertion_site, $tsd5 ];
                    push @putative_TSD5, $tsd5;
                    push @put_TSD_names5, [ $seq_name, $tsd5 ];
                    print "tsd5: $fname_start\n\n";
                    print $log_out "tsd5: $fname_start\n\n";
                    $tsd5_count++;
                    $tsd4_catch = 1;
                    undef $tsd1;
                    undef $tsd2;
                    undef $tsd3;
                    undef $tsd4;
                }
                elsif (!$tsd1 and !$tsd2 and $tsd3) {
                    if ((length($tsd5) >= length($tsd3) or $tsd5_count > $tsd3_count) and ($tsd3_count <= $tsd5_count or $tsd3_precount < $tsd5_precount)) {
                        my $insertion_site = substr($left_tsd, (-4-length($tsd5)-10), (length($tsd5)+10)) . substr($right_tsd, 3+length($tsd5), 10);
                        push @TSD_info5, [ $seq_name, $insertion_site, $tsd5 ];
                        push @putative_TSD5, $tsd5;
                        push @put_TSD_names5, [ $seq_name, $tsd5 ];
                        print "tsd5: $fname_start\n\n";
                        print $log_out "tsd5: $fname_start\n\n";
                        $tsd5_count++;
                        $tsd4_catch = 1;
                        undef $tsd1;
                        undef $tsd2;
                        undef $tsd3;
                        undef $tsd4;
                    }
                    else {
                        my $insertion_site = substr($left_tsd, (-4 - length($tsd3) - 10),(length($tsd3) + 10)) . substr($right_tsd, 4 + length($tsd3), 10);
                        push @TSD_info3, [ $seq_name, $insertion_site, $tsd3 ];
                        push @putative_TSD3, $tsd3;
                        push @put_TSD_names3, [ $seq_name, $tsd3 ];
                        print "tsd3: $fname_start\n\n";
                        print $log_out "tsd3: $fname_start\n\n";
                        $tsd3_count++;
                        undef $tsd1;
                        undef $tsd2;
                        undef $tsd4;
                        undef $tsd5;
                    }
                }
                elsif (!$tsd1_count and !$tsd3_count and $tsd2) {
                    if ((length($tsd5) >= length($tsd2) or $tsd5_count > $tsd2_count) and ($tsd2_count <= $tsd5_count or $tsd2_precount < $tsd5_precount)) {
                        my $insertion_site = substr($left_tsd, (-4-length($tsd5)-10), (length($tsd5)+10)) . substr($right_tsd, 3+length($tsd5), 10);
                        push @TSD_info5, [ $seq_name, $insertion_site, $tsd5 ];
                        push @putative_TSD5, $tsd5;
                        push @put_TSD_names5, [ $seq_name, $tsd5 ];
                        print "tsd5: $fname_start\n\n";
                        print $log_out "tsd5: $fname_start\n\n";
                        $tsd5_count++;
                        $tsd4_catch = 1;
                        undef $tsd1;
                        undef $tsd2;
                        undef $tsd3;
                        undef $tsd4;
                    }
                    else {
                        my $insertion_site = substr($left_tsd, (-14)) . substr($right_tsd, 4, 10);
                        push @TSD_info2, [ $seq_name, $insertion_site, $tsd2 ];
                        push @putative_TSD2, $tsd2;
                        push @put_TSD_names2, [ $seq_name, $tsd2 ];
                        print "tsd2: $fname_start\n\n";
                        print $log_out "tsd2: $fname_start\n\n";
                        $tsd2_count++;
                        undef $tsd1;
                        undef $tsd5;
                        undef $tsd3;
                        undef $tsd4;
                    }
                }
                elsif (!$tsd2 and !$tsd3 and $tsd1) {
                    if ((length($tsd5) >= length($tsd1) or $tsd5_count > $tsd1_count) and ($tsd1_count <= $tsd5_count or $tsd1_precount < $tsd5_precount)) {
                        my $insertion_site = substr($left_tsd, (-4-length($tsd5)-10), (length($tsd5)+10)) . substr($right_tsd, 3+length($tsd5), 10);
                        push @TSD_info5, [ $seq_name, $insertion_site, $tsd5 ];
                        push @putative_TSD5, $tsd5;
                        push @put_TSD_names5, [ $seq_name, $tsd5 ];
                        print "tsd5: $fname_start\n\n";
                        print $log_out "tsd5: $fname_start\n\n";
                        $tsd5_count++;
                        $tsd4_catch = 1;
                        undef $tsd1;
                        undef $tsd2;
                        undef $tsd3;
                        undef $tsd4;
                    }
                    else {
                        my $insertion_site = substr($left_tsd, -14, 12) . substr($right_tsd, 4, 10);
                        push @TSD_info1, [ $seq_name, $insertion_site, $tsd1 ];
                        push @putative_TSD1, $tsd1;
                        push @put_TSD_names1, [ $seq_name, $tsd1 ];
                        print "tsd1: $fname_start\n\n";
                        print $log_out "tsd1: $fname_start\n\n";
                        $tsd1_count++;
                        undef $tsd5;
                        undef $tsd2;
                        undef $tsd3;
                        undef $tsd4;
                    }
                }
            }
            else {
                if ((length($tsd5) >= length($tsd6) or $tsd5_count > $tsd6_count) and ($tsd6_count <= $tsd5_count or $tsd6_precount < $tsd5_precount)) {
                    if (!$tsd1 and !$tsd2 and !$tsd3) {
                        my $insertion_site = substr($left_tsd, (-4-length($tsd5)-10), (length($tsd5)+10)) . substr($right_tsd, 3+length($tsd5), 10);
                        push @TSD_info5, [ $seq_name, $insertion_site, $tsd5 ];
                        push @putative_TSD5, $tsd5;
                        push @put_TSD_names5, [ $seq_name, $tsd5 ];
                        print "tsd5: $fname_start\n\n";
                        print $log_out "tsd5: $fname_start\n\n";
                        $tsd5_count++;
                        $tsd4_catch = 1;
                        undef $tsd1;
                        undef $tsd2;
                        undef $tsd3;
                        undef $tsd4;
                        undef $tsd6;
                    }
                    elsif (!$tsd1 and !$tsd2 and $tsd3) {
                        if ((length($tsd5) >= length($tsd3) or $tsd5_count > $tsd3_count) and ($tsd3_count <= $tsd5_count or $tsd3_precount < $tsd5_precount)) {
                            my $insertion_site = substr($left_tsd, (-4-length($tsd5)-10), (length($tsd5)+10)) . substr($right_tsd, 3+length($tsd5), 10);
                            push @TSD_info5, [ $seq_name, $insertion_site, $tsd5 ];
                            push @putative_TSD5, $tsd5;
                            push @put_TSD_names5, [ $seq_name, $tsd5 ];
                            print "tsd5: $fname_start\n\n";
                            print $log_out "tsd5: $fname_start\n\n";
                            $tsd5_count++;
                            $tsd4_catch = 1;
                            undef $tsd1;
                            undef $tsd2;
                            undef $tsd3;
                            undef $tsd4;
                            undef $tsd6;
                        }
                        else {
                            my $insertion_site = substr($left_tsd, (-4 - length($tsd3) - 10),(length($tsd3) + 10)) . substr($right_tsd, 4 + length($tsd3), 10);
                            push @TSD_info3, [ $seq_name, $insertion_site, $tsd3 ];
                            push @putative_TSD3, $tsd3;
                            push @put_TSD_names3, [ $seq_name, $tsd3 ];
                            print "tsd3: $fname_start\n\n";
                            print $log_out "tsd3: $fname_start\n\n";
                            $tsd3_count++;
                            undef $tsd1;
                            undef $tsd2;
                            undef $tsd4;
                            undef $tsd6;
                            undef $tsd5;
                        }
                    }
                    elsif (!$tsd1 and !$tsd3 and $tsd2) {
                        if ((length($tsd5) >= length($tsd2) or $tsd5_count > $tsd2_count) and ($tsd2_count <= $tsd5_count or $tsd2_precount < $tsd5_precount)) {
                            my $insertion_site = substr($left_tsd, (-4-length($tsd5)-10), (length($tsd5)+10)) . substr($right_tsd, 3+length($tsd5), 10);
                            push @TSD_info5, [ $seq_name, $insertion_site, $tsd5 ];
                            push @putative_TSD5, $tsd5;
                            push @put_TSD_names5, [ $seq_name, $tsd5 ];
                            print "tsd5: $fname_start\n\n";
                            print $log_out "tsd5: $fname_start\n\n";
                            $tsd5_count++;
                            $tsd4_catch = 1;
                            undef $tsd1;
                            undef $tsd2;
                            undef $tsd3;
                            undef $tsd4;
                            undef $tsd6;
                        }
                        else {
                            my $insertion_site = substr($left_tsd, (-14)) . substr($right_tsd, 4, 10);
                            push @TSD_info2, [ $seq_name, $insertion_site, $tsd2 ];
                            push @putative_TSD2, $tsd2;
                            push @put_TSD_names2, [ $seq_name, $tsd2 ];
                            print "tsd2: $fname_start\n\n";
                            print $log_out "tsd2: $fname_start\n\n";
                            $tsd2_count++;
                            undef $tsd1;
                            undef $tsd4;
                            undef $tsd3;
                            undef $tsd6;
                            undef $tsd5;
                        }
                    }
                    elsif (!$tsd2 and !$tsd3 and $tsd1) {
                        if ((length($tsd5) >= length($tsd1) or $tsd5_count > $tsd1_count) and ($tsd1_count <= $tsd5_count or $tsd1_precount < $tsd5_precount)) {
                            my $insertion_site = substr($left_tsd, (-4-length($tsd5)-10), (length($tsd5)+10)) . substr($right_tsd, 3+length($tsd5), 10);
                            push @TSD_info5, [ $seq_name, $insertion_site, $tsd5 ];
                            push @putative_TSD5, $tsd5;
                            push @put_TSD_names5, [ $seq_name, $tsd5 ];
                            print "tsd5: $fname_start\n\n";
                            print $log_out "tsd5: $fname_start\n\n";
                            $tsd5_count++;
                            $tsd4_catch = 1;
                            undef $tsd1;
                            undef $tsd2;
                            undef $tsd3;
                            undef $tsd4;
                            undef $tsd6;
                        }
                        else {
                            my $insertion_site = substr($left_tsd, -14, 12) . substr($right_tsd, 4, 10);
                            push @TSD_info1, [ $seq_name, $insertion_site, $tsd1 ];
                            push @putative_TSD1, $tsd1;
                            push @put_TSD_names1, [ $seq_name, $tsd1 ];
                            print "tsd1: $fname_start\n\n";
                            print $log_out "tsd1: $fname_start\n\n";
                            $tsd1_count++;
                            undef $tsd4;
                            undef $tsd2;
                            undef $tsd3;
                            undef $tsd6;
                            undef $tsd5;
                        }
                    }
                }
                else {
                    if (!$tsd1 and !$tsd2 and !$tsd3) {
                        my $insertion_site = substr($left_tsd, (-3-length($tsd6)-10), (length($tsd6)+10)) . substr($right_tsd, 4+length($tsd6), 10);
                        push @TSD_info6, [ $seq_name, $insertion_site, $tsd6 ];
                        push @putative_TSD6, $tsd6;
                        push @put_TSD_names6, [ $seq_name, $tsd6 ];
                        print "tsd6: $fname_start\n\n";
                        print $log_out "tsd6: $fname_start\n\n";
                        $tsd6_count++;
                        $tsd4_catch = 1;
                        undef $tsd1;
                        undef $tsd2;
                        undef $tsd3;
                        undef $tsd4;
                        undef $tsd5;
                    }
                    elsif (!$tsd1 and !$tsd2 and $tsd3) {
                        if ((length($tsd6) >= length($tsd3) or $tsd6_count > $tsd3_count) and ($tsd3_count <= $tsd6_count or $tsd3_precount < $tsd6_precount)) {
                            my $insertion_site = substr($left_tsd, (-3-length($tsd6)-10), (length($tsd6)+10)) . substr($right_tsd, 4+length($tsd6), 10);
                            push @TSD_info6, [ $seq_name, $insertion_site, $tsd6 ];
                            push @putative_TSD6, $tsd6;
                            push @put_TSD_names6, [ $seq_name, $tsd6 ];
                            print "tsd6: $fname_start\n\n";
                            print $log_out "tsd6: $fname_start\n\n";
                            $tsd6_count++;
                            $tsd4_catch = 1;
                            undef $tsd1;
                            undef $tsd2;
                            undef $tsd3;
                            undef $tsd4;
                            undef $tsd5;
                        }
                        else {
                            my $insertion_site = substr($left_tsd, (-4 - length($tsd3) - 10),(length($tsd3) + 10)) . substr($right_tsd, 4 + length($tsd3), 10);
                            push @TSD_info3, [ $seq_name, $insertion_site, $tsd3 ];
                            push @putative_TSD3, $tsd3;
                            push @put_TSD_names3, [ $seq_name, $tsd3 ];
                            print "tsd3: $fname_start\n\n";
                            print $log_out "tsd3: $fname_start\n\n";
                            $tsd3_count++;
                            undef $tsd1;
                            undef $tsd2;
                            undef $tsd4;
                            undef $tsd6;
                            undef $tsd5;
                        }
                    }
                    elsif (!$tsd1 and !$tsd3 and $tsd2) {
                        if ((length($tsd6) >= length($tsd2) or $tsd6_count > $tsd2_count) and ($tsd2_count <= $tsd6_count or $tsd2_precount < $tsd6_precount)) {
                            my $insertion_site = substr($left_tsd, (-3-length($tsd6)-10), (length($tsd6)+10)) . substr($right_tsd, 4+length($tsd6), 10);
                            push @TSD_info6, [ $seq_name, $insertion_site, $tsd6 ];
                            push @putative_TSD6, $tsd6;
                            push @put_TSD_names6, [ $seq_name, $tsd6 ];
                            print "tsd6: $fname_start\n\n";
                            print $log_out "tsd6: $fname_start\n\n";
                            $tsd6_count++;
                            $tsd4_catch = 1;
                            undef $tsd1;
                            undef $tsd2;
                            undef $tsd3;
                            undef $tsd4;
                            undef $tsd5;
                        }
                        else {
                            my $insertion_site = substr($left_tsd, (-14)) . substr($right_tsd, 4, 10);
                            push @TSD_info2, [ $seq_name, $insertion_site, $tsd2 ];
                            push @putative_TSD2, $tsd2;
                            push @put_TSD_names2, [ $seq_name, $tsd2 ];
                            print "tsd2: $fname_start\n\n";
                            print $log_out "tsd2: $fname_start\n\n";
                            $tsd2_count++;
                            undef $tsd1;
                            undef $tsd4;
                            undef $tsd3;
                            undef $tsd6;
                            undef $tsd5;
                        }
                    }
                    elsif (!$tsd2 and !$tsd3 and $tsd1) {
                        if ((length($tsd6) >= length($tsd1) or $tsd6_count > $tsd1_count) and ($tsd1_count <= $tsd6_count or $tsd1_precount < $tsd6_precount)) {
                            my $insertion_site = substr($left_tsd, (-3-length($tsd6)-10), (length($tsd6)+10)) . substr($right_tsd, 4+length($tsd6), 10);
                            push @TSD_info6, [ $seq_name, $insertion_site, $tsd6 ];
                            push @putative_TSD6, $tsd6;
                            push @put_TSD_names6, [ $seq_name, $tsd6 ];
                            print "tsd6: $fname_start\n\n";
                            print $log_out "tsd6: $fname_start\n\n";
                            $tsd6_count++;
                            $tsd4_catch = 1;
                            undef $tsd1;
                            undef $tsd2;
                            undef $tsd3;
                            undef $tsd4;
                            undef $tsd5;
                        }
                        else {
                            my $insertion_site = substr($left_tsd, -14, 12) . substr($right_tsd, 4, 10);
                            push @TSD_info1, [ $seq_name, $insertion_site, $tsd1 ];
                            push @putative_TSD1, $tsd1;
                            push @put_TSD_names1, [ $seq_name, $tsd1 ];
                            print "tsd1: $fname_start\n\n";
                            print $log_out "tsd1: $fname_start\n\n";
                            $tsd1_count++;
                            undef $tsd4;
                            undef $tsd2;
                            undef $tsd3;
                            undef $tsd6;
                            undef $tsd5;
                        }
                    }
                }
            }
        }
        elsif ($tsd6) {
            if (!$tsd1 and !$tsd2 and !$tsd3) {
                my $insertion_site = substr($left_tsd, (-3-length($tsd6)-10), (length($tsd6)+10)) . substr($right_tsd, 4+length($tsd6), 10);
                push @TSD_info6, [ $seq_name, $insertion_site, $tsd6 ];
                push @putative_TSD6, $tsd6;
                push @put_TSD_names6, [ $seq_name, $tsd6 ];
                print "tsd6: $fname_start\n\n";
                print $log_out "tsd6: $fname_start\n\n";
                $tsd6_count++;
                $tsd4_catch = 1;
                undef $tsd1;
                undef $tsd2;
                undef $tsd3;
                undef $tsd4;
                undef $tsd5;
            }
            elsif (!$tsd1 and !$tsd2 and $tsd3) {
                if ((length($tsd6) >= length($tsd3) or $tsd6_count > $tsd3_count) and ($tsd3_count <= $tsd6_count or $tsd3_precount < $tsd6_precount)) {
                    my $insertion_site = substr($left_tsd, (-3-length($tsd6)-10), (length($tsd6)+10)) . substr($right_tsd, 4+length($tsd6), 10);
                    push @TSD_info6, [ $seq_name, $insertion_site, $tsd6 ];
                    push @putative_TSD6, $tsd6;
                    push @put_TSD_names6, [ $seq_name, $tsd6 ];
                    print "tsd6: $fname_start\n\n";
                    print $log_out "tsd6: $fname_start\n\n";
                    $tsd6_count++;
                    $tsd4_catch = 1;
                    undef $tsd1;
                    undef $tsd2;
                    undef $tsd3;
                    undef $tsd4;
                    undef $tsd5;
                }
                else {
                    my $insertion_site = substr($left_tsd, (-4 - length($tsd3) - 10),(length($tsd3) + 10)) . substr($right_tsd, 4 + length($tsd3), 10);
                    push @TSD_info3, [ $seq_name, $insertion_site, $tsd3 ];
                    push @putative_TSD3, $tsd3;
                    push @put_TSD_names3, [ $seq_name, $tsd3 ];
                    print "tsd3: $fname_start\n\n";
                    print $log_out "tsd3: $fname_start\n\n";
                    $tsd3_count++;
                    undef $tsd1;
                    undef $tsd2;
                    undef $tsd4;
                    undef $tsd6;
                    undef $tsd5;
                }
            }
            elsif (!$tsd1 and !$tsd3 and $tsd2) {
                if ((length($tsd6) >= length($tsd2) or $tsd6_count > $tsd2_count) and ($tsd2_count <= $tsd6_count or $tsd2_precount < $tsd6_precount)) {
                    my $insertion_site = substr($left_tsd, (-3-length($tsd6)-10), (length($tsd6)+10)) . substr($right_tsd, 4+length($tsd6), 10);
                    push @TSD_info6, [ $seq_name, $insertion_site, $tsd6 ];
                    push @putative_TSD6, $tsd6;
                    push @put_TSD_names6, [ $seq_name, $tsd6 ];
                    print "tsd6: $fname_start\n\n";
                    print $log_out "tsd6: $fname_start\n\n";
                    $tsd6_count++;
                    $tsd4_catch = 1;
                    undef $tsd1;
                    undef $tsd2;
                    undef $tsd3;
                    undef $tsd4;
                    undef $tsd5;
                }
                else {
                    my $insertion_site = substr($left_tsd, (-14)) . substr($right_tsd, 4, 10);
                    push @TSD_info2, [ $seq_name, $insertion_site, $tsd2 ];
                    push @putative_TSD2, $tsd2;
                    push @put_TSD_names2, [ $seq_name, $tsd2 ];
                    print "tsd2: $fname_start\n\n";
                    print $log_out "tsd2: $fname_start\n\n";
                    $tsd2_count++;
                    undef $tsd1;
                    undef $tsd4;
                    undef $tsd3;
                    undef $tsd6;
                    undef $tsd5;
                }
            }
            elsif (!$tsd2 and !$tsd3 and $tsd1) {
                if ((length($tsd6) >= length($tsd1) or $tsd6_count > $tsd1_count) and ($tsd1_count <= $tsd6_count or $tsd1_precount < $tsd6_precount)) {
                    my $insertion_site = substr($left_tsd, (-3-length($tsd6)-10), (length($tsd6)+10)) . substr($right_tsd, 4+length($tsd6), 10);
                    push @TSD_info6, [ $seq_name, $insertion_site, $tsd6 ];
                    push @putative_TSD6, $tsd6;
                    push @put_TSD_names6, [ $seq_name, $tsd6 ];
                    print "tsd6: $fname_start\n\n";
                    print $log_out "tsd6: $fname_start\n\n";
                    $tsd6_count++;
                    $tsd4_catch = 1;
                    undef $tsd1;
                    undef $tsd2;
                    undef $tsd3;
                    undef $tsd4;
                    undef $tsd5;
                }
                else {
                    my $insertion_site = substr($left_tsd, -14, 12) . substr($right_tsd, 4, 10);
                    push @TSD_info1, [ $seq_name, $insertion_site, $tsd1 ];
                    push @putative_TSD1, $tsd1;
                    push @put_TSD_names1, [ $seq_name, $tsd1 ];
                    print "tsd1: $fname_start\n\n";
                    print $log_out "tsd1: $fname_start\n\n";
                    $tsd1_count++;
                    undef $tsd4;
                    undef $tsd2;
                    undef $tsd3;
                    undef $tsd6;
                    undef $tsd5;
                }
            }
        }
  }

  if ($tsd1 and $tsd4_catch == 0) {
    if (!$tsd2 and !$tsd3) {
      my $insertion_site = substr($left_tsd, -14, 12) . substr($right_tsd, 4, 10);
      push @TSD_info1, [ $seq_name, $insertion_site, $tsd1 ];
      push @putative_TSD1, $tsd1;
      push @put_TSD_names1, [ $seq_name, $tsd1 ];
      #my $left_flank  = substr($starting_left_flank,  0, -2);
      #my $right_flank = substr($starting_right_flank, 2);
      #my $flanks      = $left_flank . $right_flank;
      #print $flanks_out ">$fname_start\n$flanks\n";
      print "tsd1: $fname_start\n\n";
      print $log_out "tsd1: $fname_start\n\n";
      $tsd1_count++;
    }
    elsif (!$tsd2 and $tsd3) {
      if ((length($tsd3) > length($tsd1)) or ($tsd1 =~ /TA/i and $p_type eq "hAT")) {
        my $insertion_site = substr($left_tsd, (-4 - length($tsd3) - 10),(length($tsd3) + 10)) . substr($right_tsd, 4 + length($tsd3), 10);
        push @TSD_info3, [ $seq_name, $insertion_site, $tsd3 ];
        push @putative_TSD3, $tsd3;
        push @put_TSD_names3, [ $seq_name, $tsd3 ];
        print "tsd3: $fname_start\n\n";
        print $log_out "tsd3: $fname_start\n\n";
        $tsd3_count++;
      }
      else {
        my $insertion_site = substr($left_tsd, -12) . substr($right_tsd, 2, 10);
        push @TSD_info1, [$seq_name, $insertion_site, $tsd1];
        push @putative_TSD1, $tsd1;
        push @put_TSD_names1, [ $seq_name, $tsd1 ];
        print "tsd1: $fname_start\n\n";
        print $log_out "tsd1: $fname_start\n\n";
        $tsd1_count++;
      }
    }
    elsif ($tsd2 and !$tsd3) {
      if (substr($tsd2, 0, 2) eq $tsd1 and substr($tsd2, -2) eq $tsd1) {
        my $insertion_site = substr($left_tsd, (-14)) . substr($right_tsd, 4, 10);
        push @TSD_info2, [ $seq_name, $insertion_site, $tsd2 ];
        push @putative_TSD2, $tsd2;
        push @put_TSD_names2, [ $seq_name, $tsd2 ];
        print "tsd2: $fname_start\n\n";
        print $log_out "tsd2: $fname_start\n\n";
        $tsd2_count++;
      }
      else {
        my $insertion_site = substr($left_tsd, -14, 12) . substr($right_tsd, 2, 10);
        push @TSD_info1, [ $seq_name, $insertion_site, $tsd1 ];
        push @putative_TSD1, $tsd1;
        push @put_TSD_names1, [ $seq_name, $tsd1 ];
        print "tsd1: $fname_start\n\n";
        print $log_out "tsd1: $fname_start\n\n";
        $tsd1_count++;
      }
    }
    elsif ($tsd2 and $tsd3) {
      if ( (substr($tsd3, 0, 2) eq $tsd1) and (substr($tsd3, -2) eq $tsd1)) {
        my $insertion_site = substr($left_tsd, (-4 - length($tsd3) - 10), (length($tsd3) + 10)) . substr($right_tsd, 4 + length($tsd3), 10);
        push @TSD_info3, [ $seq_name, $insertion_site, $tsd3 ];
        push @putative_TSD3, $tsd3;
        push @put_TSD_names3, [ $seq_name, $tsd3 ];
        print "tsd3: $fname_start\n\n";
        print $log_out "tsd3: $fname_start\n\n";
        $tsd3_count++;
      }
      elsif ((substr($tsd2, 0, 2) eq $tsd1) and (substr($tsd2, -2) eq $tsd1)) {
        my $insertion_site = substr($left_tsd,(-2-length($tsd2)-10), (length($tsd2) + 10)) . substr($right_tsd, 2 +length($tsd2), 10);
        push @TSD_info2, [ $seq_name, $insertion_site, $tsd2 ];
        push @putative_TSD2, $tsd2;
        push @put_TSD_names2, [ $seq_name, $tsd2 ];
        print "tsd2: $fname_start\n\n";
        print $log_out "tsd2: $fname_start\n\n";
        $tsd2_count++;
      }
      else {
        my $insertion_site = substr($left_tsd, -14, 12) . substr($right_tsd, 4, 10);
        push @TSD_info1, [ $seq_name, $insertion_site, $tsd1 ];
        push @putative_TSD1, $tsd1;
        push @put_TSD_names1, [ $seq_name, $tsd1 ];
        print "tsd1: $fname_start\n\n";
        print $log_out "tsd1: $fname_start\n\n";
        $tsd1_count++;
      }
    }
  }
  elsif ($tsd3 and $tsd4_catch == 0) {
    if (!$tsd2) {
      my $insertion_site = substr($left_tsd, (-4-length($tsd3)-10), (length($tsd3)+10)) . substr($right_tsd, 4+length($tsd3), 10);
      push @TSD_info3, [ $seq_name, $insertion_site, $tsd3 ];
      push @putative_TSD3, $tsd3;
      push @put_TSD_names3, [ $seq_name, $tsd3 ];
      print "tsd3: $fname_start\n\n";
      print $log_out "tsd3: $fname_start\n\n";
      $tsd3_count++;
    }
    else {
      if ((length($tsd3) > length($tsd2)) and (substr($tsd3, 0, 4) eq $tsd2) and (substr($tsd3, -4) eq $tsd2)) {
        my $insertion_site = substr($left_tsd, (-4-length($tsd3)-10), (length($tsd3)+10)) . substr($right_tsd, 4+length($tsd3), 10);
        push @TSD_info3, [ $seq_name, $insertion_site, $tsd3 ];
        push @putative_TSD3, $tsd3;
        push @put_TSD_names3, [ $seq_name, $tsd3 ];
        print "tsd3: $fname_start\n\n";
        print $log_out "tsd3: $fname_start\n\n";
        $tsd3_count++;
      }
      else {
        my $insertion_site = substr($left_tsd, -14) . substr($right_tsd, 4, 10);
        push @TSD_info2, [ $seq_name, $insertion_site, $tsd2 ];
        push @putative_TSD2, $tsd2;
        push @put_TSD_names2, [ $seq_name, $tsd2 ];
        print "tsd2: $fname_start\n\n";
        print $log_out "tsd2: $fname_start\n\n";
        $tsd2_count++;
      }
    }
  }
  elsif ($tsd2 and $tsd4_catch == 0) {
    my $insertion_site = substr($left_tsd, -14) . substr($right_tsd, 4, 10);
    push @TSD_info2, [ $seq_name, $insertion_site, $tsd2 ];
    push @putative_TSD2, $tsd2;
    push @put_TSD_names2, [ $seq_name, $tsd2 ];
    print "tsd2: $fname_start\n\n";
    print $log_out "tsd2: $fname_start\n\n";
    $tsd2_count++;
  }
  elsif (!$tsd1 and !$tsd2 and !$tsd3 and !$tsd4 and !$tsd5 and !$tsd6) {
    push @no_TSD_found, $seq_name;
    print "no TSD: $fname_start\n\n";
    print $log_out "no TSD: $fname_start\n\n";
  }
}
print "\ntsd4 count: $tsd4_count\ntsd5 count: $tsd5_count\ntsd6 count: $tsd6_count\ntsd1 count: $tsd1_count\ntsd2 count: $tsd2_count\ntsd3 count: $tsd3_count\n\n";
print $log_out "\ntsd4 count: $tsd4_count\ntsd5 count: $tsd5_count\ntsd6 count: $tsd6_count\ntsd1 count: $tsd1_count\ntsd2 count: $tsd2_count\ntsd3 count: $tsd3_count\n\n";

print "Finished code to find TSDs\n\n";
print $log_out "Finished code to find TSDs\n\n";

my $final_align_len = $final_aln_obj->num_sequences();
$element_info{"copy_num"} = $final_align_len;
my $final_aln_out = File::Spec->catpath($volume, $out_path, $filename . ".final");
$out = Bio::AlignIO->new(
  -file             => ">$final_aln_out",
  -format           => 'fasta',
  -displayname_flat => 0
);
$out->write_aln($final_aln_obj);

my $no_TSD_found_len = @no_TSD_found;
if ($no_TSD_found_len >= 1) {
  print "Copies without TSDs printed to file\n";
  print $log_out "Copies without TSDs printed to file\n";
  foreach my $item (@no_TSD_found) {
    print $no_TSD_found_out "$item\tNo TSDs found\n";
  }
}

TIR_Adjustment:
#if necessary, change TIRS to remove the 2bp or 4bp TSD, or 1bp partial TSD, that are included in them
#calculate the overall percent identity of each TIR
if (($tsd4_count > $tsd1_count) and ($tsd4_count > $tsd2_count) and ($tsd4_count > $tsd3_count) and ($tsd4_count > $tsd5_count) and ($tsd4_count > $tsd6_count)) {
    print "Adjusting TIRs by 1bp\n\n";
    print $log_out "Adjusting TIRs by 1bp\n\n";

    $left_tir_start += 1;
    $right_tir_start -= 1;
    my $last_count;
    foreach my $seq_obj ($final_aln_obj->each_seq()) {
        $last_count++;
        my $seq_name = $seq_obj->id();
        if (!defined $tir_positions{$seq_name}{'left_tir_start'}) {
            my $left_pos_obj  = $seq_obj->location_from_column($left_tir_start);
            my $right_pos_obj = $seq_obj->location_from_column($right_tir_start);
            my $left_pos      = $left_pos_obj->start();
            my $right_pos     = $right_pos_obj->start();
            $tir_positions{$seq_name}{'left_tir_start'}  = $left_pos;
            $tir_positions{$seq_name}{'right_tir_start'} = $right_pos;
        }
        $element_info{$seq_name}{"left_tir_start"} = $tir_positions{$seq_name}{'left_tir_start'} + 1;
        $element_info{$seq_name}{"right_tir_start"} = $tir_positions{$seq_name}{'right_tir_start'} - 1;
    }

    $left_TIR_aln_obj = $final_aln_obj->slice($left_tir_start, $left_tir_end, 1);
    $right_TIR_aln_obj = $final_aln_obj->slice($right_tir_end, $right_tir_start, 1);
    $element_aln_obj = $final_aln_obj->slice($left_tir_start, $right_tir_start, 1);
    @putative_TSD = @putative_TSD4;
    @TSD_info = @TSD_info4;
    @put_TSD_names = @put_TSD_names4;

    undef @putative_TSD1;
    undef @TSD_info1;
    undef @put_TSD_names1;
    undef @putative_TSD2;
    undef @TSD_info2;
    undef @put_TSD_names2;
    undef @putative_TSD3;
    undef @TSD_info3;
    undef @put_TSD_names3;
    undef @putative_TSD5;
    undef @TSD_info5;
    undef @put_TSD_names5;
    undef @putative_TSD6;
    undef @TSD_info6;
    undef @put_TSD_names6;
    goto Summation;
    
}
TIR_Adjustment2: if (($tsd5_count > $tsd1_count) and ($tsd5_count > $tsd3_count) and ($tsd5_count > $tsd2_count) and ($tsd5_count > $tsd6_count)) {
    print "Adjusting right TIR by 1bp\n\n";
    print $log_out "Adjusting right TIR by 1bp\n\n";
    my $matches = 0;
    my $mismatches = 0;
    
    my $test_left_tir_start = $left_tir_start;
    my $test_right_tir_start = $right_tir_start - 1;
    for (my $count = 0 ; $count < 2 ; $count++) {
        foreach my $seq_obj ($final_aln_obj->each_seq()) {
            my $seq = $seq_obj->seq();
            $seq =~ s/-//g;
            my $left_pos_obj  = $seq_obj->location_from_column($test_left_tir_start);
            my $right_pos_obj = $seq_obj->location_from_column($test_right_tir_start);
            my $left_pos      = $left_pos_obj->start();
            my $right_pos     = $right_pos_obj->start();
            my $left_nt = substr($seq, $left_pos-1+$count, 1);
            my $right_nt = substr($seq, $right_pos-1-$count, 1);
            $right_nt =~ tr/ATGCatgc/TACGtacg/;
            $right_nt = reverse($right_nt);
            if ($left_nt eq $right_nt) {
                $matches++;
            }
            else {
                $mismatches++;
            }
        }
        if ($mismatches > 0 and $matches/$mismatches < 0.8) {
            goto TIR_Adjustment3;
        }
    }

    $left_tir_start += 0;
    $right_tir_start -= 1;
    my $last_count;
    foreach my $seq_obj ($final_aln_obj->each_seq()) {
        $last_count++;
        my $seq_name = $seq_obj->id();
        if (!defined $tir_positions{$seq_name}{'left_tir_start'}) {
            my $left_pos_obj  = $seq_obj->location_from_column($left_tir_start);
            my $right_pos_obj = $seq_obj->location_from_column($right_tir_start);
            my $left_pos      = $left_pos_obj->start();
            my $right_pos     = $right_pos_obj->start();
            $tir_positions{$seq_name}{'left_tir_start'}  = $left_pos;
            $tir_positions{$seq_name}{'right_tir_start'} = $right_pos;
        }
        $element_info{$seq_name}{"right_tir_start"} = $tir_positions{$seq_name}{'right_tir_start'} - 1;
    }

    $left_TIR_aln_obj = $final_aln_obj->slice($left_tir_start, $left_tir_end, 1);
    $right_TIR_aln_obj = $final_aln_obj->slice($right_tir_end, $right_tir_start, 1);
    $element_aln_obj = $final_aln_obj->slice($left_tir_start, $right_tir_start, 1);
    @putative_TSD = @putative_TSD5;
    @TSD_info = @TSD_info5;
    @put_TSD_names = @put_TSD_names5;

    undef @putative_TSD1;
    undef @TSD_info1;
    undef @put_TSD_names1;
    undef @putative_TSD2;
    undef @TSD_info2;
    undef @put_TSD_names2;
    undef @putative_TSD3;
    undef @TSD_info3;
    undef @put_TSD_names3;
    undef @putative_TSD4;
    undef @TSD_info4;
    undef @put_TSD_names4;
    undef @putative_TSD6;
    undef @TSD_info6;
    undef @put_TSD_names6;
    goto Summation;
    
}
TIR_Adjustment3: if (($tsd6_count > $tsd1_count) and ($tsd6_count > $tsd3_count) and ($tsd6_count > $tsd2_count)) {
    print "Adjusting left TIR by 1bp\n\n";
    print $log_out "Adjusting left TIR by 1bp\n\n";
    
    my $matches = 0;
    my $mismatches = 0;
    
    my $test_left_tir_start = $left_tir_start + 1;
    my $test_right_tir_start = $right_tir_start;
    
    for (my $count = 0 ; $count < 2 ; $count++) {
        foreach my $seq_obj ($final_aln_obj->each_seq()) {
            my $seq = $seq_obj->seq();
            $seq =~ s/-//g;
            my $left_pos_obj  = $seq_obj->location_from_column($test_left_tir_start);
            my $right_pos_obj = $seq_obj->location_from_column($test_right_tir_start);
            my $left_pos      = $left_pos_obj->start();
            my $right_pos     = $right_pos_obj->start();
            my $left_nt = substr($seq, $left_pos-1+$count, 1);
            my $right_nt = substr($seq, $right_pos-1-$count, 1);
            $right_nt =~ tr/ATGCatgc/TACGtacg/;
            $right_nt = reverse($right_nt);
            if ($left_nt eq $right_nt) {
                $matches++;
            }
            else {
                $mismatches++;
            }
        }
        if ($mismatches > 0 and $matches/$mismatches < 0.8) {
            goto TIR_Adjustment4;
        }
    }
    
    $left_tir_start += 1;
    $right_tir_start -= 0;
    my $last_count;
    foreach my $seq_obj ($final_aln_obj->each_seq()) {
        $last_count++;
        my $seq_name = $seq_obj->id();
        if (!defined $tir_positions{$seq_name}{'left_tir_start'}) {
            my $left_pos_obj  = $seq_obj->location_from_column($left_tir_start);
            my $right_pos_obj = $seq_obj->location_from_column($right_tir_start);
            my $left_pos      = $left_pos_obj->start();
            my $right_pos     = $right_pos_obj->start();
            $tir_positions{$seq_name}{'left_tir_start'}  = $left_pos;
            $tir_positions{$seq_name}{'right_tir_start'} = $right_pos;
        }
        $element_info{$seq_name}{"right_tir_start"} = $tir_positions{$seq_name}{'right_tir_start'} - 1;
    }

    $left_TIR_aln_obj = $final_aln_obj->slice($left_tir_start, $left_tir_end, 1);
    $right_TIR_aln_obj = $final_aln_obj->slice($right_tir_end, $right_tir_start, 1);
    $element_aln_obj = $final_aln_obj->slice($left_tir_start, $right_tir_start, 1);
    @putative_TSD = @putative_TSD6;
    @TSD_info = @TSD_info6;
    @put_TSD_names = @put_TSD_names6;

    undef @putative_TSD1;
    undef @TSD_info1;
    undef @put_TSD_names1;
    undef @putative_TSD2;
    undef @TSD_info2;
    undef @put_TSD_names2;
    undef @putative_TSD3;
    undef @TSD_info3;
    undef @put_TSD_names3;
    undef @putative_TSD5;
    undef @TSD_info5;
    undef @put_TSD_names5;
    undef @putative_TSD4;
    undef @TSD_info4;
    undef @put_TSD_names4;
    goto Summation;
    
}
TIR_Adjustment4: if (($tsd1_count > $tsd2_count) and ($tsd1_count > $tsd3_count) and ($tsd1_count > $tsd5_count) and ($tsd1_count > $tsd6_count)) {
  print "Adjusting TIRs by 2bp\n\n";
  print $log_out "Adjusting TIRs by 2bp\n\n";
  $left_tir_start += 2;
  $right_tir_start -= 2;
  my $last_count;
  foreach my $seq_obj ($final_aln_obj->each_seq()) {
    $last_count++;
    my $seq_name = $seq_obj->id();
    if (!defined $tir_positions{$seq_name}{'left_tir_start'}) {
      my $left_pos_obj  = $seq_obj->location_from_column($left_tir_start);
      my $right_pos_obj = $seq_obj->location_from_column($right_tir_start);
      my $left_pos      = $left_pos_obj->start();
      my $right_pos     = $right_pos_obj->start();
      $tir_positions{$seq_name}{'left_tir_start'}  = $left_pos;
      $tir_positions{$seq_name}{'right_tir_start'} = $right_pos;
    }
    $element_info{$seq_name}{"left_tir_start"} = $tir_positions{$seq_name}{'left_tir_start'} + 2;
    $element_info{$seq_name}{"right_tir_start"} = $tir_positions{$seq_name}{'right_tir_start'} - 2;
  }

  $left_TIR_aln_obj = $final_aln_obj->slice($left_tir_start, $left_tir_end, 1);
  $right_TIR_aln_obj = $final_aln_obj->slice($right_tir_end, $right_tir_start, 1);
  $element_aln_obj = $final_aln_obj->slice($left_tir_start, $right_tir_start, 1);
  @putative_TSD = @putative_TSD1;
  @TSD_info = @TSD_info1;
  @put_TSD_names = @put_TSD_names1;
  undef @putative_TSD4;
  undef @TSD_info4;
  undef @put_TSD_names4;
  undef @putative_TSD2;
  undef @TSD_info2;
  undef @put_TSD_names2;
  undef @putative_TSD3;
  undef @TSD_info3;
  undef @put_TSD_names3;
  undef @putative_TSD5;
  undef @TSD_info5;
  undef @put_TSD_names5;
  undef @putative_TSD6;
  undef @TSD_info6;
  undef @put_TSD_names6;
  goto Summation;
}
TIR_Adjustment5: if (($tsd2_count > $tsd1_count) and ($tsd2_count > $tsd3_count) and ($tsd2_count > $tsd5_count) and ($tsd2_count > $tsd6_count)) {
  print "Adjusting TIRs by 4bp\n\n";
  print $log_out "Adjusting TIRs by 4bp\n\n";
  $left_tir_start += 4;
  $right_tir_start -= 4;
  my $last_count;

  foreach my $seq_obj ($final_aln_obj->each_seq()) {
    $last_count++;
    my $seq_name = $seq_obj->id();
    if (!defined $tir_positions{$seq_name}{'left_tir_start'}) {
      my $left_pos_obj  = $seq_obj->location_from_column($left_tir_start);
      my $right_pos_obj = $seq_obj->location_from_column($right_tir_start);
      my $left_pos      = $left_pos_obj->start();
      my $right_pos     = $right_pos_obj->start();
      $tir_positions{$seq_name}{'left_tir_start'}  = $left_pos;
      $tir_positions{$seq_name}{'right_tir_start'} = $right_pos;
    }
    $element_info{$seq_name}{"left_tir_start"} = $tir_positions{$seq_name}{'left_tir_start'} + 4;
    $element_info{$seq_name}{"right_tir_start"} = $tir_positions{$seq_name}{'right_tir_start'} - 4;
  }

  $left_TIR_aln_obj = $final_aln_obj->slice($left_tir_start, $left_tir_end, 1);
  $right_TIR_aln_obj = $final_aln_obj->slice($right_tir_end, $right_tir_start, 1);
  $element_aln_obj = $final_aln_obj->slice($left_tir_start, $right_tir_start, 1);
  @putative_TSD = @putative_TSD2;
  @TSD_info = @TSD_info2;
  @put_TSD_names = @put_TSD_names2;
  undef @putative_TSD1;
  undef @TSD_info1;
  undef @put_TSD_names1;
  undef @putative_TSD4;
  undef @TSD_info4;
  undef @put_TSD_names4;
  undef @putative_TSD3;
  undef @TSD_info3;
  undef @put_TSD_names3;
  undef @putative_TSD5;
  undef @TSD_info5;
  undef @put_TSD_names5;
  undef @putative_TSD6;
  undef @TSD_info6;
  undef @put_TSD_names6;
  goto Summation;
}

TIR_Adjustment6:  foreach my $seq_obj ($final_aln_obj->each_seq()) {
my $seq_name = $seq_obj->id();
if (!defined $tir_positions{$seq_name}{'left_tir_start'}) {
  my $left_pos_obj  = $seq_obj->location_from_column($left_tir_start);
  my $right_pos_obj = $seq_obj->location_from_column($right_tir_start);
  my $left_pos      = $left_pos_obj->start();
  my $right_pos     = $right_pos_obj->start();
  $tir_positions{$seq_name}{'left_tir_start'}  = $left_pos;
  $tir_positions{$seq_name}{'right_tir_start'} = $right_pos;
}
$element_info{$seq_name}{"left_tir_start"} = $tir_positions{$seq_name}{'left_tir_start'};
$element_info{$seq_name}{"right_tir_start"} = $tir_positions{$seq_name}{'right_tir_start'};
}
@putative_TSD = @putative_TSD3;
@TSD_info = @TSD_info3;
@put_TSD_names = @put_TSD_names3;
undef @putative_TSD1;
undef @TSD_info1;
undef @put_TSD_names1;
undef @putative_TSD2;
undef @TSD_info2;
undef @put_TSD_names2;
undef @putative_TSD4;
undef @TSD_info4;
undef @put_TSD_names4;
undef @putative_TSD5;
undef @TSD_info5;
undef @put_TSD_names5;
undef @putative_TSD6;
undef @TSD_info6;
undef @put_TSD_names6;

Summation: my $left_tir_id   = $left_TIR_aln_obj->percentage_identity;
my $right_tir_id  = $right_TIR_aln_obj->percentage_identity;
my $element_id    = $element_aln_obj->percentage_identity;

my $element_aln_out = File::Spec->catpath($volume, $out_path, $filename . "_element-only-msa.fa");
$out = Bio::AlignIO->new(-file => ">$element_aln_out", -format => 'fasta', -displayname_flat => 0);
$out->write_aln($element_aln_obj);
my $trimal_element_aln_out = File::Spec->catpath($volume, $out_path, $filename . "_element-only-msa_trimal.fa");
system("$trimal -in $element_aln_out -out $trimal_element_aln_out -gt 0.5");
my $in_obj = Bio::AlignIO->new(-file => $trimal_element_aln_out, -format => 'fasta' , -alphabet => 'dna');
my $trimmed_element_aln_obj = $in_obj->next_aln();
my $element_consensus = $trimmed_element_aln_obj->consensus_string();
my $left_tir_seq  = $left_TIR_aln_obj->consensus_string();
my $right_tir_seq = $right_TIR_aln_obj->consensus_string();

my $left_tir_iupac = $left_TIR_aln_obj->consensus_iupac();
$left_tir_iupac = uc($left_tir_iupac);
my $left_tir_iupac_seq_obj = Bio::Seq->new(-seq => $left_tir_iupac, -alphabet => 'dna');
my $left_tir_iupac_obj   = Bio::Tools::IUPAC->new(-seq => $left_tir_iupac_seq_obj);
my $left_tir_iupac_regexp  = $left_tir_iupac_obj->regexp();
my $left_tir_regex_out_path = File::Spec->catpath($volume, $out_path, $filename . "_left_tir_regex.info");
open(my $left_tir_regex_out, ">", $left_tir_regex_out_path);
print $left_tir_regex_out "$fname_fin\t$left_tir_iupac_regexp\n";
close($left_tir_regex_out);

my $right_tir_iupac = $right_TIR_aln_obj->consensus_iupac();
$right_tir_iupac = uc($right_tir_iupac);
my $right_tir_iupac_seq_obj = Bio::Seq->new(-seq => $right_tir_iupac, -alphabet => 'dna');
my $right_tir_iupac_obj   = Bio::Tools::IUPAC->new(-seq => $right_tir_iupac_seq_obj);
my $right_tir_iupac_regexp  = $right_tir_iupac_obj->regexp();
my $right_tir_regex_out_path = File::Spec->catpath($volume, $out_path, $filename . "_right_tir_regex.info");
open(my $right_tir_regex_out, ">", $right_tir_regex_out_path);
print $right_tir_regex_out "$fname_fin\t$right_tir_iupac_regexp\n";
close($right_tir_regex_out);

my $both_tirs_regex = $left_tir_iupac_regexp . ".+" . $right_tir_iupac_regexp;
my $both_tirs_regex_out_path = File::Spec->catpath($volume, $out_path, $filename . "_both_tirs_regex.info");
open(my $both_tirs_regex_out, ">", $both_tirs_regex_out_path);
print $both_tirs_regex_out "$fname_fin\t$both_tirs_regex\n";
close($both_tirs_regex_out);

$element_info{"element_id"}    = $element_id;
$element_info{"element_consensus"} = $element_consensus;
$element_info{"left_tir_seq"}  = $left_tir_seq;
$element_info{"left_tir_id"}   = $left_tir_id;
$element_info{"right_tir_seq"} = $right_tir_seq;
$element_info{"right_tir_id"}  = $right_tir_id;

my $left_tir_out_path = File::Spec->catpath($volume, $out_path, $filename . ".left-tir.fa");
my $right_tir_out_path = File::Spec->catpath($volume, $out_path, $filename . ".right-tir.fa");
open(my $left_tir_out, ">", $left_tir_out_path) or die "Error creating $left_tir_out_path. $!\n";
print $left_tir_out join("\n", ">" . $fname_fin . "_left_tir_consensus", $left_tir_seq), "\n";
close($left_tir_out);

open(my $right_tir_out, ">", $right_tir_out_path) or die "Error creating $right_tir_out_path. $!\n";
print $right_tir_out join("\n", ">" . $fname_fin . "_right_tir_consensus", $right_tir_seq), "\n";
close($right_tir_out);

my $element_consensus_out_path = File::Spec->catpath($volume, $out_path, $filename . ".consensus");
open(my $element_consensus_out, ">", $element_consensus_out_path) or die "Error creating $element_consensus_out_path. $!\n";
print $element_consensus_out join("\n", ">" . $fname_fin . "_consensus", $element_info{'element_consensus'}), "\n";
close($element_consensus_out);
print "Printing fasta\n";
print $log_out "Printing fasta\n";
my $element_fasta_out_path = File::Spec->catpath($volume, $out_path, $filename . ".fasta");
print_fasta($element_fasta_out_path, $element_aln_obj);


#Determine the most common TSD by sequence first if possible or by length. If >80% of the TSDs are the same, then that sequence is stored as the TSD for output. Otherwise, look at the lengths of the TSDs and store if >80% of the TSDs have the same length.
my %TSD_counts;
my $TSD_array_length = @putative_TSD;
if ($TSD_array_length == 0) {
  print "No TSDs found\n";
  print $log_out "No TSDs found\n";
  my $bad_out_path = File::Spec->catpath($volume, $out_path, $filename . ".no_tsd");
  error_out($bad_out_path, "$filename\tNo TSDs found for any copy", 1);
  my $gff_path = File::Spec->catpath($volume, $out_path, $filename . ".gff");
  if (defined $protein) {
      generate_gff($final_aln_obj, $gff_path, "TSD", \%element_info, $lflank_len, $rflank_len, $p_type);
  }
  else {
      generate_gff($final_aln_obj, $gff_path, "TSD", \%element_info, $flank, $flank);
  }
  exit 0;
}

#count the occurances of a TSD sequence and sort with the highest value first in an array
foreach my $row_ref (@put_TSD_names) {
  my @pos = @{$row_ref};
  $TSD_counts{ $pos[1] }++;
  $element_info{ $pos[0] }{"TSD"} = $pos[1];
}

#undef @put_TSD_names;
my @sorted_TSD_keys = sort { $TSD_counts{$b} <=> $TSD_counts{$a} } keys(%TSD_counts);
my $final_TSD_length;
my $final_TSD_seq;
my $final_TSD_fraction;
my $need_consensus = 0;
my @good_TSD_length;
my @final_tsd_info;
$final_align_len = $final_aln_obj->num_sequences();
$element_info{"copy_num"} = $final_align_len;
my $tsd_counts = @sorted_TSD_keys;
my $d;

foreach my $row (@sorted_TSD_keys) {
    $d++;
    my $TSD_fraction = $TSD_counts{$row} / $final_align_len;
    if ($TSD_fraction > 0.8 ) {
        $final_TSD_seq = $row;
        $final_TSD_length = length($row);
        $final_TSD_fraction = $TSD_fraction;
        push @good_TSD_length, length($row);
        push @final_tsd_info, [length($row), $row];
    }
    elsif ($TSD_fraction >= 0.01) {
        if (defined $final_TSD_seq) {
            $final_TSD_seq = $final_TSD_seq . ", " . $row;
            my $row_len = length($row);
            $final_TSD_length = $final_TSD_length . ", " . $row_len;
            $final_TSD_fraction = $final_TSD_fraction . ", " . $TSD_fraction;
            push @good_TSD_length, $row_len;
            push @final_tsd_info, [$row_len, $row];
        }
        else {
            if (!defined $final_TSD_seq) {
                $need_consensus = 1;
                last;
            }
        }
    }
    else {
        if ($d == 1) {
            print "TSDs found in fewer than 1% of copies, treating as not having any.\n";
            print $log_out "TSDs found in fewer than 1% of copies, treating as not having any.\n";
            my $bad_out_path = File::Spec->catpath($volume, $out_path, $filename . ".no_tsd");
            error_out($bad_out_path, "$filename\tTSDs in fewer than 1% of copies, treating as not having any.", 1);
            my $gff_path = File::Spec->catpath($volume, $out_path, $filename . ".gff");
            if (defined $protein) {
                generate_gff($final_aln_obj, $gff_path, "TSD", \%element_info, $lflank_len, $rflank_len, $p_type);
            }
            else {
                generate_gff($final_aln_obj, $gff_path, "TSD", \%element_info, $flank, $flank);
            }
            exit 0;
        }
        else {
            last;
        }
    }
}
if ($need_consensus == 0) {
    $element_info{"TSD_fraction"} = $final_TSD_fraction;
    $element_info{"TSD_seq"}      = $final_TSD_seq;
    $element_info{"TSD_len"}      = $final_TSD_length;
}
else {
    my %TSD_length_counts;
    foreach my $row (@putative_TSD) {
        $TSD_length_counts{length($row)}++;
    }
    my @sorted_TSD_length_keys = sort { $TSD_length_counts{$b} <=> $TSD_length_counts{$a} } keys(%TSD_length_counts);
    foreach my $tsd_len (@sorted_TSD_length_keys) {
        my $TSD_fraction = $TSD_length_counts{$tsd_len} / $final_align_len;
        if ($TSD_fraction > 0.01 ) {
            if (!defined $final_TSD_fraction) {
                $final_TSD_length = $tsd_len;
                $final_TSD_fraction = $TSD_fraction;
                push @good_TSD_length, $tsd_len;
            }
            else {
                if ($TSD_fraction > 0.05 ) {
                    $final_TSD_length = $final_TSD_length . ", " . $tsd_len;
                    $final_TSD_fraction = $final_TSD_fraction . ", " . $TSD_fraction;
                    push @good_TSD_length, $tsd_len;
                }
                else {
                    last;
                }
            }
        }
        else {
            if (!defined $final_TSD_fraction) {
                $final_TSD_length = 'NA';
                $final_TSD_fraction = 'NA';
                last;
            }
        }
    }
    $element_info{"TSD_len"}      = $final_TSD_length;
    $element_info{"TSD_fraction"} = $final_TSD_fraction;
}
#undef @putative_TSD;

foreach my $row (@good_TSD_length) {
    my $insertion_site_file = $filename . ".insertion-site" . $row . ".fa";
    my $insertion_site_out_path = File::Spec->catpath($volume, $out_path, $insertion_site_file);
    my $tsd_out_path = File::Spec->catpath($volume, $out_path, $filename . ".tsd" . $row . ".fa");

    open(my $tsd_info_out, ">", $tsd_out_path) or die "Error creating $tsd_out_path. $!\n";
    open(my $insertion_site_out, ">", $insertion_site_out_path) or die "Error creating $insertion_site_out_path. $!\n";
    foreach my $row_ref (@TSD_info) {
        #print "Printing insertion site info\n";
        #print $log_out "Printing insertion site info\n";
        my @pos = @{$row_ref};
        if (length($pos[2]) == $row) {
            print $tsd_info_out ">$pos[0]\n$pos[2]\n";
            if ($pos[1] !~ m/n/i) {
                print $insertion_site_out ">$pos[0]\n$pos[1]\n";
            }
        }
    }
    close($insertion_site_out);
    close($tsd_info_out);
    close($no_TSD_found_out);
    my $out_fix           = $out_path . "/";
    my $Blogo_config_path = $FindBin::Bin . "/../lib/blogo/Blogo.conf";
    my $Blogo_path        = $FindBin::Bin . "/../lib/blogo";
    system("$Blogo_path/Blogo_batch.pl file_path=$out_fix file_names=$insertion_site_file img_abs_dir=$out_fix conf_file=$Blogo_config_path");

    if ($need_consensus == 1) {
        #read in tsd file as MSA to generate IUPAC consensus
        my $tsd_in_obj =  Bio::AlignIO->new(-file => $tsd_out_path, -format => 'fasta');
        my $tsd_aln_obj   = $tsd_in_obj->next_aln();
        my $tsd_consensus = $tsd_aln_obj->consensus_iupac();
        if (!exists $element_info{"TSD_seq"} ) {
            $element_info{"TSD_seq"} = $tsd_consensus;
            push @final_tsd_info, [$row, $tsd_consensus];
        }
        else {
            $element_info{"TSD_seq"} = $element_info{"TSD_seq"} . ", " . $tsd_consensus;
            push @final_tsd_info, [$row, $tsd_consensus];
        }
    }
}

Classification:

#import list of transposase superfamilies found in genome MSA is from if specified
my %known_tpase;
if ($protein_match) {
    open(my $known_in, "<", $protein_match) or die "Error reading $protein_match . $!";
    while (my $line = <$known_in>) {
        chomp $line;
        my @split = split("\t", $line);
        my $super_name = $split[0];
        if ($super_name ne "retro" and $super_name ne "Unknown" and $super_name ne "unpredicted" and $super_name ne "Maverick" and $super_name ne "Helitron") {
            $known_tpase{$super_name} = 1;
        }
    }
}


#import TE characteristics from table and generate regular expressions to classify element by TIR and TSD if possible
print "\nProcessing TE charateristics table\n";
print $log_out "\nProcessing TE charateristics table\n";
my $TE_char_path = $FindBin::Bin . "/DNA_TE_TIR-TSD.table";
open(my $TE_in, "<", $TE_char_path) or die "Error reading $TE_char_path . $!";
my %element_char_hash;

while (my $line = <$TE_in>) {
  chomp $line;
  my @split = split("\t", $line);
  my $ele_name = $split[0];
  if (!$protein_match) {
      $known_tpase{$ele_name} = 1;
  }
  if ($split[1] ne ".") {
    my $seq     = $split[1];
    my $seq_obj = Bio::Seq->new(-seq => $seq, -alphabet => 'dna');
    my $iupac   = Bio::Tools::IUPAC->new(-seq => $seq_obj);
    my $regexp  = $iupac->regexp();
    $element_char_hash{$ele_name}{'tir_con'} = $regexp;
  }
  else {
    $element_char_hash{$ele_name}{'tir_con'} = '';
  }
  if ($split[2] =~ m/,/) {
    my @tsd_lengths   = split(',', $split[2]);
    my $tsd_array_len = @tsd_lengths;
    my $count         = 1;
    my $regexp        = '';
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
      my @tsd_cons     = split(',', $split[3]);
      my $tsd_cons_len = @tsd_cons;
      my $count        = 1;
      my $regexp       = '';
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
print "Checking element for DNA TE charateristics\n";
print $log_out "Checking element for DNA TE charateristics\n";
#print "final_tsd_info array dump\n";
#print $log_out "final_tsd_info array dump\n";
#print Dumper(\@final_tsd_info);

foreach my $tsd_info_ref (@final_tsd_info) {
    my %element_hits;
    my @tsd_info = @{$tsd_info_ref};
    foreach my $ele_name (keys %element_char_hash) {
        if (exists $known_tpase{$ele_name}) {
            if ($tsd_info[0] =~ m/$element_char_hash{$ele_name}{"tsd_length"}/) {
                $element_hits{$ele_name}++;
                if ( $element_char_hash{$ele_name}{"tsd_con"} ne '') {
                    if ($tsd_info[1] =~ m/$element_char_hash{$ele_name}{"tsd_con"}/i) {
                        $element_hits{$ele_name}++;
                    }
                    else {
                        delete $element_hits{$ele_name};
                        next;
                    }
                }
                if ($element_char_hash{$ele_name}{"tir_con"} ne '') {
                    my $right_match = $element_info{'right_tir_seq'};
                    $right_match =~ tr/ATGCatgc/TACGtacg/;
                    $right_match = reverse($right_match);
                    if (($element_info{'left_tir_seq'} =~ m/^$element_char_hash{$ele_name}{"tir_con"}/i) or ($right_match =~ m/^$element_char_hash{$ele_name}{"tir_con"}/i)) {
                        $element_hits{$ele_name}++;
                    }
    
                    else {
                        delete $element_hits{$ele_name};
                        next;
                    }
                }
                else {
                    next;
                }
            }
            else {
                next;
            }
        }
    }

    my @sorted;
    foreach my $key (sort { $element_hits{$b} <=> $element_hits{$a} } keys(%element_hits)) {
        my $count = $element_hits{$key};
        push @sorted, [ $key, $count ];
    }

    my $classification = '';

    #store all classifications and the number of hits to each
    foreach my $row_ref (@sorted) {
        my @info = @{$row_ref};
        if ($classification eq '') {
            $classification = $classification . $info[0] . "_" . $info[1];
        }
        else {
            $classification = $classification . ", " . $info[0] . "_" . $info[1];
        }
    }
    if ($classification eq '' or !defined $classification) {
        $classification = "Unknown";
    }
    if (!exists $element_info{"classification"}) {
        $element_info{"classification"} = $classification;
    }
    else {
        if ($element_info{"classification"} !~ m/$classification/) {
            $element_info{"classification"} =  $element_info{"classification"} . "|" . $classification
        }
    }
}
if (defined $protein and defined $p_type) {
    $element_info{"classification"} = $p_type . "||" . $element_info{"classification"};
}
print "Element classification finished\n";
print $log_out "Element classification finished\n";
my $element_info_out_path = File::Spec->catpath($volume, $out_path, $filename . ".element_info");
open(my $element_info_out, ">", $element_info_out_path)
  or die "Error creating $element_info_out_path. $!\n";
print $element_info_out join("\t", $fname_fin, $element_info{'copy_num'}, $element_id, $element_info{'left_tir_seq'},
  $element_info{'left_tir_id'}, $element_info{'right_tir_seq'}, $element_info{'left_tir_id'}, $element_info{'TSD_len'},
  $element_info{'TSD_seq'},     $element_info{'TSD_fraction'}, $element_info{"classification"}), "\n";
close($element_info_out);

print "\nSetting up gff path\n";
print $log_out "\nSetting up gff path\n";
my $gff_path = File::Spec->catpath($volume, $out_path, $filename . ".gff");
print "Calling gff printer\n";
print $log_out "Calling gff printer\n";
if (defined $protein) {
    generate_gff($final_aln_obj, $gff_path, "final", \%element_info, $lflank_len, $rflank_len, $p_type);
}
else {
    generate_gff($final_aln_obj, $gff_path, "final", \%element_info, $flank, $flank);
}
print "Exited gff printer\n";
print $log_out "Exited gff printer\n";

if (!defined $all) {
  print "Cleaning up files\n";
  print $log_out "Cleaning up files\n";
  clean_files($out_path);
}
exit 0;

#--------------------------Subroutines---------------------------------#

sub match_tirs {
  my $self       = shift;    ## seq_obj
  my $input_path = shift;
  my $round      = shift;
  $self->throw("Need Bio::LocatableSeq argument")
    unless ref $self && $self->isa('Bio::LocatableSeq');
  if (!-f $input_path) {
    die "Supplied filepath is not valid";
  }
  my $seq = $self->seq();
  $seq =~ s/-//g;
  my $seq_len = length($seq);
  my $check_seq = $seq;
  $check_seq =~ tr/ATGCatgc/TACGtacg/;
  $check_seq = reverse($check_seq);
  my $half_seq = substr($seq, 0, $seq_len/2);
  my $seq_name = $self->id();
  my $fa_aln_obj = Bio::SearchIO->new(-format => 'fasta', -file => $input_path);
  my @result;
  #print "ggsearch results path: $input_path\n";
  #print $log_out "ggsearch results path: $input_path\n";

  #go through the FASTA input object . . down to the HSP
  while (my $result = $fa_aln_obj->next_result) {

    #print "numHits: ",$result->num_hits,"\n";
    #print $log_out "numHits: ",$result->num_hits,"\n";
    if ($result->num_hits == 0) {
      push @result, (0, [0]);
      last;
    }
    while (my $hit = $result->next_hit) {
      while (my $hsp = $hit->next_hsp) {

        #grab the query, hit, and homology strings
        my $homo_string = $hsp->homology_string;
        my $query_str   = $hsp->query_string;
        my $hit_str     = $hsp->hit_string;
        my $len_homo    = length($homo_string);
        my $len_query   = length($query_str);
        my $len_hit     = length($hit_str);

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
        my $match_cutoff = 9;
        if ($round == 2) {
            $match_cutoff = 10;
        }
        my $first_match = 0;

        #parse homology string, keeping track of match length and mismatches or gaps
        for (my $count = 0 ; $count < length($homo_string) ; $count++) {
          my $homo_char  = substr($homo_string, $count, 1);
          my $query_char = substr($query_str,   $count, 1);
          my $hit_char   = substr($hit_str,     $count, 1);
          if ($round == 1 and $count == 8 and $total_mis_aln >= 4) {
            if ($match_len < 3) {
              $match_len = 0;
              $start_pos = '';
              $match_query = '';
              $match_hit = '';
              $end_pos = '';
              #print "No TIRs found near start of sequences, resetting counts and ending\n";
              last;
            }
          }
          if ($round == 2) {
            if (($count == 20 and $total_mis_aln >= 12) or ($count == 30 and $match_len <= $match_cutoff) ) {
                $match_len = 0;
                $start_pos = '';
                $match_query = '';
                $match_hit = '';
                $end_pos = '';
                last;
            }
            elsif ($match_len == 4 and $match_mis_aln > 1) {
                $match_len = 0;
                $start_pos = '';
                $match_query = '';
                $match_hit = '';
                $end_pos = '';
                last;
            }
          }
          ## skip any seqs that have 1 or more mismatches in the first 3 bases of the TIR
          #if ($round == 3 or $round == 1) {
          if ($round == 1) {
            if (($count == 3 and $total_mis_aln >= 1) or ($count == 6 and $total_mis_aln >= 2)) {
              $match_len = 0;
              $start_pos = '';
              $match_query = '';
              $match_hit = '';
              $end_pos = '';
              last;
            }
          }
          if ($round == 3) {
              if ($count == 9 and $total_mis_aln >= 4) {
                  $match_len = 0;
                  $start_pos = '';
                  $match_query = '';
                  $match_hit = '';
                  $end_pos = '';
                  last;
              }
          }

          if ($match_len == 0) {
            #if match length equals 0 and position is not a match, continue to next position
            if ($homo_char eq " ") {
                if ($round == 1 or $round == 2) {
                    if ($hit_char ne '-' and $query_char ne '-') {
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
                else {
                    $total_mis_aln++;
                    next;
                }
            }

            #if position is a match, store info, continue to next position
            elsif ($homo_char eq ":") {
              if ($hit_char eq 'N' or $query_char eq 'N') {
                  next;
              }
              $start_pos = $count;
              $last_good = $count;
              $match_len++;
              $match_query .= $query_char;
              $match_hit .= $hit_char;
              if ($first_match == 0) {
                  $first_match = 1;
              }

              #print "Initial match at $start_pos\n";
              next;
            }
          }
          elsif ($match_len >= 1 and $match_len <= ($match_cutoff - 1)) {

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

                # print "First Mismatch at $count\n";
                next;
              }

              #more than one mismatch, reset counters and other info, continue
              elsif ($match_mis_aln > 1 and $match_len < 5) {
                $match_len = 0;
                $start_pos = '';
                $match_query = '';
                $match_hit = '';
                $match_mis_aln = 0;

                #print "Another Mismatch at $count, resetting counts\n";
                next;
              }
              elsif ($match_mis_aln < 3 and $match_len >= 5) {
                $match_len++;
                $last_good = $count;
                $match_query .= $query_char;
                $match_hit .= $hit_char;

                # print "Another Mismatch at $count\n";
                next;
              }
              else {
                if ($round ==1) {
                    $match_len = 0;
                    $start_pos = '';
                    $match_query = '';
                    $match_hit = '';
                    $end_pos = '';
                    #print "Another Mismatch at $count, resetting counts and ending\n";
                    last;
                }
                else {
                    $match_len = 0;
                    $start_pos = '';
                    $match_query = '';
                    $match_hit = '';
                    $match_mis_aln = 0;
                    #print "Another Mismatch at $count, resetting counts\n";
                    next;
                }
              }
            }

            #position is a match, store info and continue
            elsif ($homo_char eq ":") {
              $last_good = $count;
              $match_len++;
              $match_query .= $query_char;
              $match_hit .= $hit_char;

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
                if ($count == (length($homo_string)-1)) {
                    $end_pos = $last_good;
                    $match_query =~ s/-//g;
                    $match_hit   =~ s/-//g;
                    
                    
                    #$match_query =~ tr/ATGCatgc/TACGtacg/;
                    #$match_query = reverse($match_query);
                    my $match_query_len = length($match_query);
                    my $match_hit_len   = length($match_hit);
                    
                    #find the position in the full sequence of the hit and query match sequences
                    $hit_pos = index(uc($half_seq), uc($match_hit)) + 1;
                    #$hit_pos = index(uc($seq), uc($match_hit), 40) + 1;
                    my $initial_query_pos = index(uc($check_seq), uc($match_query));
                    $query_pos = $seq_len - $initial_query_pos;
                    
                    #print "\nSeq: $seq\n";
                    #reverse complement the match query sequence
                    $match_query =~ tr/ATGCatgc/TACGtacg/;
                    $match_query = reverse($match_query);
                    
                    print "1st catch:\n$seq_name hit_pos:$hit_pos query_pos:$query_pos  $match_hit  $match_query\n";
                    print $log_out "1st catch:\n$seq_name hit_pos:$hit_pos query_pos:$query_pos  $match_hit  $match_query\n";
                    #store sequence name and the hit and query info
                    my @match = (
                    $seq_name,
                    {
                        "hit"   => [ $hit_pos,   $match_hit_len ],
                        "query" => [ $query_pos, $match_query_len ],
                        "seq" => [$match_hit, $match_query]
                    }
                    );
                    ## add a catch for $hit_pos or $query_pos == 0 [index returns -1 when $match_hit was not found in $seq]
                    if ($hit_pos == 0 or $query_pos == 0) {
                        $match_len   = 0;
                        $start_pos   = '';
                        $match_query = '';
                        $match_hit   = '';
                        $end_pos     = '';
                        last;
                    }
                    else {
                        push @result, (1, [@match]);
                    }
                    last;
                }

                # print "Another Mismatch at $count, proceeding\n";
                # print $log_out "Another Mismatch at $count, proceeding\n";
                next;
              }

              #mismatches 3 or more, store final info for alignment match and end parsing
              elsif ($match_mis_aln >= 3) {
                $end_pos = $last_good;
                $match_query =~ s/-//g;
                $match_hit   =~ s/-//g;

                #$match_query =~ tr/ATGCatgc/TACGtacg/;
                #$match_query = reverse($match_query);
                my $match_query_len = length($match_query);
                my $match_hit_len   = length($match_hit);
                
                #find the position in the full sequence of the hit and query match sequences
                $hit_pos = index(uc($half_seq), uc($match_hit)) + 1;
                #$hit_pos = index(uc($seq), uc($match_hit), 40) + 1;
                my $initial_query_pos = index(uc($check_seq), uc($match_query));
                $query_pos = $seq_len - $initial_query_pos;
                    
                #print "\nSeq: $seq\n";
                #reverse complement the match query sequence
                $match_query =~ tr/ATGCatgc/TACGtacg/;
                $match_query = reverse($match_query);
                
                print "2nd catch:\n$seq_name hit_pos:$hit_pos query_pos:$query_pos  $match_hit  $match_query\n";
                print $log_out "2nd catch:\n$seq_name hit_pos:$hit_pos query_pos:$query_pos  $match_hit  $match_query\n";
                #store sequence name and the hit and query info
                my @match = (
                  $seq_name,
                  {
                    "hit"   => [ $hit_pos,   $match_hit_len ],
                    "query" => [ $query_pos, $match_query_len ],
                    "seq" => [$match_hit, $match_query]
                  }
                );
                ## add a catch for $hit_pos or $query_pos == 0 [index returns -1 when $match_hit was not found in $seq]
                if ($hit_pos == 0 or $query_pos == 0) {
                    $match_len   = 0;
                    $start_pos   = '';
                    $match_query = '';
                    $match_hit   = '';
                    $end_pos     = '';
                    last;
                }
                else {
                    push @result, (1, [@match]);
                }
                last;
              }
            }

            #position is a match, store info and continue
            elsif ($homo_char eq ":") {
              $last_good = $count;
              $match_len++;
              $match_query .= $query_char;
              $match_hit   .= $hit_char;
              if ($count == (length($homo_string)-1)) {
                    $end_pos = $last_good;
                    $match_query =~ s/-//g;
                    $match_hit   =~ s/-//g;
                    
                    #$match_query =~ tr/ATGCatgc/TACGtacg/;
                    #$match_query = reverse($match_query);
                    my $match_query_len = length($match_query);
                    my $match_hit_len   = length($match_hit);
                    
                    #find the position in the full sequence of the hit and query match sequences
                    $hit_pos = index(uc($half_seq), uc($match_hit)) + 1;
                    #$hit_pos = index(uc($seq), uc($match_hit), 40) + 1;
                    my $initial_query_pos = index(uc($check_seq), uc($match_query));
                    $query_pos = $seq_len - $initial_query_pos;
                    
                    #print "\nSeq: $seq\n";
                    #reverse complement the match query sequence
                    $match_query =~ tr/ATGCatgc/TACGtacg/;
                    $match_query = reverse($match_query);
                    
                    print "3rd catch:\n$seq_name hit_pos:$hit_pos query_pos:$query_pos  $match_hit  $match_query\n";
                    print $log_out "3rd catch:\n$seq_name hit_pos:$hit_pos query_pos:$query_pos  $match_hit  $match_query\n";
                    #store sequence name and the hit and query info
                    my @match = (
                    $seq_name,
                    {
                        "hit"   => [ $hit_pos,   $match_hit_len ],
                        "query" => [ $query_pos, $match_query_len ],
                        "seq" => [$match_hit, $match_query]
                    }
                    );
                    ## add a catch for $hit_pos or $query_pos == 0 [index returns -1 when $match_hit was not found in $seq]
                    if ($hit_pos == 0 or $query_pos == 0) {
                        $match_len   = 0;
                        $start_pos   = '';
                        $match_query = '';
                        $match_hit   = '';
                        $end_pos     = '';
                        last;
                    }
                    else {
                        push @result, (1, [@match]);
                    }
                    last;
              }
              #print "Another match at $count. Length is $match_len\n";
              #print $log_out "Another match at $count. Length is $match_len\n";
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
  ## generate_gff ($aln_obj,$path,$round(final/TSD))
  my $self         = shift;    ## $aln_obj
  my $path         = shift;
  my $round        = shift;
  my $ele_info_ref = shift;
  my $left_flank = shift;
  my $right_flank = shift;
  my $type = shift;
  
  print "Entered .gff printer\n";
  print $log_out "Entered .gff printer\n";
  $self->throw("Need Bio::Align::AlignI argument")
    unless ref $self && $self->isa('Bio::Align::AlignI');

  #not fully implemented: round changes the output based on when the call is made
  if ($round eq 'final' or $round eq 'TSD') {
    open(my $out, ">", $path) or die "Error creating $path. $!\n";
    print "Round = final\n";
    print $log_out "Round = final\n";
    my $count = 1;
    foreach my $seq_obj ($self->each_seq()) {
      my $seq_name = $seq_obj->id();
      my $seq      = $seq_obj->seq();
      $seq =~ s/-//g;
      my $seq_len    = length($seq);
      my $ori_end    = $seq_len - $right_flank;
      my $left_comp  = $ele_info_ref->{$seq_name}{"left_tir_start"} - ($left_flank+1);
      my $right_comp = $ele_info_ref->{$seq_name}{"right_tir_start"} - $ori_end;
      my $copy_num;
      my $eleid;
      my $seqid;
      my $type = 'terminal_inverted_repeat_element';
      my $start;
      my $end;
      my $strand;
      my $input;

      #grab copy information from TARGeT output
      if ($seq_name =~ /^([0-9]*).+_Query:(.*)_Sbjct:(.*)_Length.+Location:\(([0-9]*)_\-_([0-9]*)\)_Direction:(.+)/) {
        $copy_num = $1;
        $eleid    = $2;
        $seqid    = $3;
        $start = $4;
        $end = $5;
        $strand   = $6;
        if ($strand eq 'plus') {
            $strand = '+';
            $start = $start + $left_comp;
            $end = $end + $right_comp;
        }
        else {
            $strand = '-';
            $start = $start + $right_comp;
            $end = $end + $left_comp;
        }
      }
      elsif ($seq_name =~ /(.+hit[0-9]+)_(.+)_([0-9]+)_([0-9]+)_(plus|minus)/) {
        $input = $1;
        $seqid = $2;
        $start = $3;
        $end = $4;
        $strand = $5;
        $copy_num = $count;
        if ($filename =~ /(cluster[0-9]+)_/) {
            $eleid = $1;
        }
        if ($strand eq 'plus') {
            $strand = '+';
            $start = $start + $left_comp;
            $end = $end + $right_comp;
        }
        else {
            $strand = '-';
            $start = $start + $right_comp;
            $end = $end + $left_comp;
        }
          
          
      }

      #grab copy information from RSPB output
      elsif ($seq_name =~ /^([0-9]*)_(.+)_(.+)-(.+)_(.+)/) {
        $copy_num = $1;
        $seqid    = $2;
        $start    = $3 + $left_comp;
        $end      = $4 + $right_comp;
        $eleid    = $5;
        $strand   = "?";
      }
      else {
        print "Header doesn't match TARGeT or RSPB:  $seq_name\n";
        print $log_out "Header doesn't match TARGeT or RSPB:  $seq_name\n";
        next;
      }
      if ($eleid =~ m/(.+)_TSD/ or $eleid =~ m/(.+)_Unknow/) {
        $eleid = $1;
      }
      my $ltir_end = $start + length($ele_info_ref->{"left_tir_seq"}) - 1;
      my $rtir_start = $end - length($ele_info_ref->{"right_tir_seq"});
      my $ele_id    = $ele_info_ref->{"element_id"};
      my $tir_id    = $ele_info_ref->{"left_tir_id"};
      my $ele_class = $type . '||N/A';
      my $tsd_frac  = 'N/A';
      my $tsd_con   = 'N/A';

      if ($round eq 'final') {
        $ele_class = $ele_info_ref->{"classification"};
        $tsd_frac  = $ele_info_ref->{"TSD_fraction"};
        $tsd_con   = $ele_info_ref->{"TSD_seq"};
      }
      my $Name = $eleid . "Copy" . $copy_num;
      if (defined $input) {
          $Name = $input;
      }
      print $out join("\t", $seqid, $PROGRAM_NAME, $type, $start, $end, '.', $strand, '.', join(";",  "ID=$eleid-$copy_num", "Name=$Name", "element_id=$ele_id",  "element_classification=$ele_class","tir_id=$tir_id", "tsd_fraction=$tsd_frac", "tsd_consensus=$tsd_con")), "\n";
      print $out join("\t", $seqid, $PROGRAM_NAME, "five_prime_terminal_inverted_repeat", $start, $ltir_end, ".", ".", ".", "Parent=$eleid-$copy_num"), "\n";
      print $out join("\t", $seqid, $PROGRAM_NAME, "three_prime_terminal_inverted_repeat", $rtir_start, $end, ".", ".", ".", "Parent=$eleid-$copy_num"), "\n";

      if ($round eq 'final' and defined $ele_info_ref->{$seq_name}{"TSD"}) {
        my $ltsd_start = $start - length($ele_info_ref->{$seq_name}{"TSD"});
        my $ltsd_end   = $start - 1;
        my $rtsd_start = $end + 1;
        my $rtsd_end   = $end + length($ele_info_ref->{$seq_name}{"TSD"});
        print $out join("\t", $seqid, $PROGRAM_NAME, "target_site_duplication",, $ltsd_start,
          $ltsd_end, ".", ".", ".", "Derives_from=$eleid-$copy_num"), "\n";
        print $out join("\t", $seqid, $PROGRAM_NAME, "target_site_duplication", $rtsd_start,
          $rtsd_end, ".", ".", ".", "Derives_from=$eleid-$copy_num"), "\n";
      }
      $count++;
    }
    close($out);
  }
}

sub clean_files {
  my $out_path = shift;
  opendir(my $in_DIR, $out_path) or die "Cannot open directory: $!";
  while (my $file = readdir($in_DIR)) {
    next if ($file =~ m/^\./);
    if ($file =~ m/\.(final|info|log|fa|element_info|removed_sequences|fasta|no_tsd|abort|gff|tif|jpg|full_id|consensus|conserved_left_flank|conserved_right_flank|conserved_flanks)$/) {
      next;
    }
    else {
      my $remove_path = File::Spec->catpath($volume, $out_path, $file);
      system("rm $remove_path");
    }
  }
  closedir($in_DIR);
}

sub get_columns {
  #This will get the column numbers of the start and stop position of both TIRs using their nucleotide positions from a hash.
  #Initially, the hash contains information for the sequences present after TIR matching. Sequences may be added subsequently.
  #This is needed after sequences and any resulting gap only columns have been removed from the alignment
  my $self          = shift;    # aln_obj
  my $tir_positions = shift;    # \%tir_positions
  my $round         = shift;    # 1=starts no ends, 2=starts and ends

  $self->throw("Need Bio::Align::AlignI argument")
    unless ref $self && $self->isa('Bio::Align::AlignI');
  my ($left_tir_start, $left_tir_end, $right_tir_start, $right_tir_end) = (0, 0, 0, 0);
  

  foreach my $seq_obj ($self->each_seq()) {
    my $seq_name = $seq_obj->id();
    #skip sequences not in the hash
    if (exists $tir_positions->{$seq_name}) {
      $left_tir_start = $self->column_from_residue_number($seq_name, $tir_positions->{$seq_name}{'left_tir_start'});
      $right_tir_start = $self->column_from_residue_number($seq_name, $tir_positions->{$seq_name}{'right_tir_start'});
      if ($round == 2 and !exists $tir_positions->{$seq_name}{'left_tir_end'} and !exists $tir_positions->{$seq_name}{'right_tir_end'}) {
        delete $tir_positions->{$seq_name};
      }
      elsif ($round == 2) {
        ##print "Round = $round\n";
        ##print $log_out "Round = $round\n";
        $left_tir_end = $self->column_from_residue_number($seq_name, $tir_positions->{$seq_name}{'left_tir_end'});
        $right_tir_end = $self->column_from_residue_number($seq_name, $tir_positions->{$seq_name}{'right_tir_end'});
        last;
      }
    }
  }

  #print "in get_col subrout: left tir start: $left_tir_start, right_tir_start: $right_tir_start, left_tir_end: $left_tir_end, right_tir_end: $right_tir_end\n";
  #print $log_out "in get_col subrout: left tir start: $left_tir_start, right_tir_start: $right_tir_start, left_tir_end: $left_tir_end, right_tir_end: $right_tir_end\n";
  return ($left_tir_start, $right_tir_start, $left_tir_end, $right_tir_end);
}

sub adjust_tir_starts {
  ### adjust tir starts.
  print "Attempting to adjust tir starts\n";
  print $log_out "Attempting to adjust tir starts\n";
  my $aln_obj         = shift;
  my $left_tir_start  = shift;
  my $right_tir_start = shift;
  my $round = shift;
  my $catch = 0.8;
  if (!defined $round) {
      $round = 0;
  }
  if ($round == 2) {
      $catch = 0.9;
  }
  print "In adjust_tir left_tir_start = $left_tir_start  right_tir_start = $right_tir_start\n";
  print $log_out "In adjust_tir left_tir_start = $left_tir_start  right_tir_start = $right_tir_start\n";

  my $num_seqs = $aln_obj->num_sequences;
  my %count;
  my %modified;
  my $orig_right_tir_start = $right_tir_start;
  my $orig_left_tir_start  = $left_tir_start;
  foreach my $seq_obj ($aln_obj->each_seq) {
    my $seq    = $seq_obj->seq;
    my $seq_id = $seq_obj->id;
    my ($gaps_left, $gaps_right) = ('', '');
    ## for right TIR
    my $left2right = $seq_obj->subseq($left_tir_start, $right_tir_start);
    my $next2end;
    if ($right_tir_start < length $seq) {
        $next2end = $seq_obj->subseq($right_tir_start + 1, length $seq);
    }
    if ($next2end =~ /^(-+)/) {
      $gaps_right = $1;
      $next2end =~ s/^(-+)//;
    }
    my @next2end = split '', $next2end;
    for (my $i = 0 ; $i < @next2end ; $i++) {
      my $nt = $next2end[$i];
      $count{right}{$i}{$nt}++;
    }

    ## for left_TIR
    my $start2left;
    if ($left_tir_start > 1) {
        $start2left = $seq_obj->subseq(1, $left_tir_start - 1);
    }
    
    if ($start2left =~ /(-+)$/) {
      $gaps_left = $1;
      $start2left =~ s/(-+)$//;
      #print "after:  $start2left\n";
      #print $log_out "after:  $start2left\n";
    }
    my @gaps = split '', $gaps_left;
    my @start2left = split '', $start2left;
    unshift @start2left , @gaps;
    for (my $i = ((scalar @start2left) - 1) ; $i >= 0 ; $i--) {
      my $nt = $start2left[$i];
      $count{left}{$i}{$nt}++;
      #print substr($seq_id,0,5) ," $i $nt " , $count{left}{$i}{$nt} , "\n";
    }
    my $new_seq = $gaps_left . $start2left . $left2right . $next2end . $gaps_right;
    $modified{$seq_id} = $new_seq;
  }

  #warn Dumper \%count;
  my $done = 0;
  ## for right tir
  foreach my $pos (sort { $a <=> $b } keys %{ $count{right} }) {
    my $nt = (sort { $count{right}{$pos}{$b} <=> $count{right}{$pos}{$a} }
        keys %{ $count{right}{$pos} } )[0];
    if ($nt eq '-') {
        last;
    }
    my $nt_count   = $count{right}{$pos}{$nt};
    my $percent_nt = $nt_count / $num_seqs;

    #print "adjust_tirs RT $nt($pos):nt_count (nextposIN:$count{left}{$pos-1}) %=$percent_nt\n";
    #print $log_out "adjust_tirs RT $nt($pos):nt_count (nextposIN:$count{left}{$pos-1}) %=$percent_nt\n";
    if ($percent_nt >= $catch) {
      $right_tir_start++;
    }
    else {
      $done = 1;
    }
    last if $done;
  }
  $done = 0;
  ## for left tir
  foreach my $pos (sort { $b <=> $a } keys %{ $count{left} }) {
    my $nt = (sort{$count{left}{$pos}{$b} <=> $count{left}{$pos}{$a}} keys %{$count{left}{$pos}})[0];
    if ($nt eq '-') {
        last;
    }
    my $nt_count   = $count{left}{$pos}{$nt};
    my $percent_nt = $nt_count / $num_seqs;

    #print "adjust_tirs LT $nt($pos):nt_count  ($nt_count/$num_seqs)=$percent_nt%\n";
    #print $log_out "adjust_tirs LT $nt($pos):nt_count  ($nt_count/$num_seqs)=$percent_nt%\n";
    if ($percent_nt >= $catch) {
      $left_tir_start--;
    }
    else {
      $done = 1;
    }
    last if $done;
  }

  if ($orig_right_tir_start != $right_tir_start or $orig_left_tir_start != $left_tir_start) {
    print "origR: $orig_right_tir_start\n";
    print $log_out "origR: $orig_right_tir_start\n";
    print "after adjR: $right_tir_start\n";
    print $log_out "after adjR: $right_tir_start\n";
    print "origL: $orig_left_tir_start\n";
    print $log_out "origL: $orig_left_tir_start\n";
    print "after adjL: $left_tir_start\n";
    print $log_out "after adjL: $left_tir_start\n";
    foreach my $seq_obj ($aln_obj->each_seq) {
      ## replace org seq with mod seq
      my $seq_id = $seq_obj->id;
      if (exists $modified{$seq_id}) {
        my $new_seq = $modified{$seq_id};
        $seq_obj->seq($new_seq);
      }
    }
  }
  return ($left_tir_start, $right_tir_start);
}

sub consensus_filter {
#### start
  my $gap_seq_pos_remove = shift;              ## \@gap_seq_pos_remove
  my $aln_obj            = shift;
  my $left_tir_start     = shift;
  my $right_tir_start    = shift;
  my $tir_positions      = shift;              # \%tir_positions
  my $round              = shift;
  my $try = shift;
  if (!defined $try) {
      $try = 1;
  }

  my $aln_len            = $aln_obj->length;
  my $num_seqs = $aln_obj->num_sequences;
  if ($num_seqs <= 1) {
    print "consensus_filter: <= 1 sequences\n";
    print $log_out "consensus_filter: <= 1 sequences\n";
    return ($left_tir_start, $right_tir_start, $gap_seq_pos_remove, $aln_obj, $tir_positions);
  }
  $round = defined $round ? $round : 'other';
  print "in consensu_filt sub: round=$round\n";
  print $log_out "in consensu_filt sub: round=$round\n";
  my %trim_gap_seq_remove;
  ## this will generate a sequence with '?' at every position in which there is less
  ## than 80% consensus
  my $consensus = $aln_obj->consensus_string(80);
  
  print "consensus:\n$consensus\n";
  print $log_out "consensus:\n$consensus\n";
  
  my $limit;
  if (defined $protein) {
      if ($num_seqs <= 6) {
          $limit = 16;
      }
      else {
        $limit = 13;
    }
  }
  else {
    if ($num_seqs <= 6) {
        $limit = 13;
    }
    else {
        $limit = 10;
    }
  }
  print "Limit: $limit\n";
  print "left_tir_start: $left_tir_start  right_tir_start: $right_tir_start\n";
  ## first round of tir finding
  if ($left_tir_start == 0 && $right_tir_start == 0 or defined $protein) {
    my @consensus    = split '', $consensus;
    my $nt_count     = 0;
    my $bad_count = 0;
    my $last_nomatch = 0;

    ## left tir
    for (my $i = 0 ; $i < ((scalar @consensus) * .8) ; $i++) {
      my $nt = $consensus[$i];
      if ($nt eq '?' and $nt_count < $limit) {
        if ($bad_count > 0) {
            if ($bad_count == 1) {
                if  ($try == 2) {
                    if ($nt_count >= 3) {
                        $bad_count++;
                        $nt_count++;
                    }
                    else {
                        $last_nomatch = $i;
                        $nt_count     = 0;
                        $bad_count = 0;
                    }
                }
                else {
                    $last_nomatch = $i;
                    $nt_count     = 0;
                    $bad_count = 0;
                }
            }
            else {
                $last_nomatch = $i;
                $nt_count     = 0;
                $bad_count = 0;
            }
        }
        else {
            if ($try == 2) {
                if ($nt_count >= 3) {
                    $bad_count++;
                    $nt_count++;
                }
                else {
                    $last_nomatch = $i;
                    $nt_count     = 0;
                }
            }
            else {
                if ($nt_count > 3) {
                    $bad_count++;
                    $nt_count++;
                }
                else {
                    $last_nomatch = $i;
                    $nt_count = 0;
                }
            }
        }
      }
      elsif ($nt ne '?' and $nt_count < $limit) {
        if ($nt eq 'n') {
            $last_nomatch = $i;
            $nt_count     = 0;
            $bad_count = 0;
            next;
        }
        $nt_count++;
      }
      else {    #($nt ne '?' and $nt_count => 3){
        $left_tir_start = $last_nomatch + 2;
        last;
      }
    }
    ## right tir
    $last_nomatch = 0;
    $nt_count     = 0;
    $bad_count = 0;
    for (my $i = ((scalar @consensus)-1); $i > ((scalar @consensus) * .2); $i--) {
      my $nt = $consensus[$i];
      if ($nt eq '?' and $nt_count < $limit) {
        if ($bad_count > 0) {
            if ($bad_count == 1) {
                if  ($try == 2) {
                    if ($nt_count >= 3) {
                        $bad_count++;
                        $nt_count++;
                    }
                    else {
                        $last_nomatch = $i;
                        $nt_count     = 0;
                        $bad_count = 0;
                    }
                }
                else {
                    $last_nomatch = $i;
                    $nt_count     = 0;
                    $bad_count = 0;
                }
            }
            else {
                $last_nomatch = $i;
                $nt_count     = 0;
                $bad_count = 0;
            }
        }
        else {
            if ($try == 2) {
                if ($nt_count >= 3) {
                    $bad_count++;
                    $nt_count++;
                }
                else {
                    $last_nomatch = $i;
                    $nt_count     = 0;
                }
            }
            else {
                if ($nt_count > 3) {
                    $bad_count++;
                    $nt_count++;
                }
                else {
                    $last_nomatch = $i;
                    $nt_count = 0;
                }
            }
        }
      }
      elsif ($nt ne '?' and $nt_count < $limit) {
        if ($nt eq 'n') {
            $last_nomatch = $i;
            $nt_count     = 0;
            $bad_count = 0;
            next;
        }
        $nt_count++;
      }
      else {    #($nt ne '?' and $nt_count => 3){
        $right_tir_start = $last_nomatch;
        last;
      }
    }
  }
  print "in con_fil sub: leftTIR: $left_tir_start\n";
  print $log_out "in con_fil sub: leftTIR: $left_tir_start\n";
  print "in con_fil sub: rightTIR: $right_tir_start\n";
  print $log_out "in con_fil sub: rightTIR: $right_tir_start\n";
  
  my $aln_depth = $aln_obj->num_sequences;
  my %con_mm;
  my %to_remove;
  for (my $i = 0 ; $i < $aln_len ; $i++) {
    if (($i >= $left_tir_start and $i <= $left_tir_start + 25) or ($i >= $right_tir_start - 25 and $i <= $right_tir_start)) {
      ## check each col for seqs that have a run of too many mismatches
      my $col_aln_obj = $aln_obj->slice($i + 1, $i + 1, 1);
      foreach my $seq_obj ($col_aln_obj->each_seq()) {
        my $seq_name = $seq_obj->id;
        my $nt       = $seq_obj->seq;
        my $con_nt   = substr($consensus, $i, 1);
        next if $con_nt =~ /N/i;
        if ($nt ne $con_nt) {
          $con_mm{$seq_name}++;
        }
      }
      foreach my $key (keys %con_mm) {
        next if $con_mm{$key} < 5;
        my $seq_obj_to_remove = $aln_obj->get_seq_by_id($key);
        $trim_gap_seq_remove{$key} = $seq_obj_to_remove;
        my @info = ($key, $i + 1, $seq_obj_to_remove);
        push @$gap_seq_pos_remove, [@info];
        $to_remove{$key} = 1;
      }
      ## end check each col for too many mismatches
    }
  }
  my $to_remove_count = keys %to_remove;
  print "in con_fil sub before get_tir_nt_starts: leftTIR: $left_tir_start\n";
  print $log_out "in con_fil sub before get_tir_nt_starts: leftTIR: $left_tir_start\n";
  print "in con_fil sub before get_tir_nt_starts: rightTIR: $right_tir_start\n";
  print $log_out "in con_fil sub before get_tir_nt_starts: rightTIR: $right_tir_start\n";
  ## skip seq removal if the number to remove is about the depth of the aln
  print "there are $aln_depth seqs in aln before cons_filter\n";
  print $log_out "there are $aln_depth seqs in aln before cons_filter\n";
  print "and there will be  $aln_depth - $to_remove_count = ", $aln_depth - $to_remove_count, " after filtering\n";
  print $log_out "and there will be  $aln_depth - $to_remove_count = ",
    $aln_depth - $to_remove_count, " after filtering\n";
  if ($aln_depth > ($to_remove_count + 5)) {
    print "removing seqs with cons_filter\n";
    print $log_out "removing seqs with cons_filter\n";
    foreach my $key (keys %trim_gap_seq_remove) {
      my $seq_obj = $trim_gap_seq_remove{$key};
      if (defined $seq_obj) {
        $aln_obj->remove_seq($seq_obj);
      }
    }
  }
  ##how close are the first round of TIR starts to the ends of the aln
  my $consensus_len = $right_tir_start - $left_tir_start + 1;
  
  my $message = '';
  my $catch = 0;
  if ($consensus_len >= ($aln_len - 50)){ ## ($consensus_len > ($aln_len - ($flank*2) + 50))
    $message = "$filename: too much of the alignment is conserved. Flanks are too similar\n";
    $catch = 1;
  }
  elsif (!$left_tir_start and !$right_tir_start) {
      $message = "$filename: no TIR start was found on either end\n";
  }
  elsif (!$left_tir_start) {
    $message = "$filename: no TIR start was found on left end\n";
  }
  elsif (!$right_tir_start){
    $message = "$filename: no TIR start was found on right end\n";
  }
  if ($message){
    if ($try == 1) {
        goto Cleaning_MSA;
    }
    if ($catch == 0) {
        my $abort_out_path = File::Spec->catpath($volume, $out_path, $filename . ".abort");
        error_out($abort_out_path, $message, 0);
    }
    else {
        my $abort_out_path = File::Spec->catpath($volume, $out_path, $filename . ".conserved_flanks");
        error_out($abort_out_path, $message, 0);
    }
  }
  my $ref2remove_these;
  if ($round ne 'final') {
    ($tir_positions,$ref2remove_these) = get_tir_nt_starts($aln_obj, $tir_positions, $left_tir_start, $right_tir_start);
  }
  foreach my $key (keys %$ref2remove_these) {
    my $seq_obj = $$ref2remove_these{$key};
    if (exists $tir_positions{$key}) {
        delete $tir_positions{$key};
    }
    if (defined $seq_obj) {
        $aln_obj->remove_seq($seq_obj);
        my @info = ($key, 0, $seq_obj);
        push @$gap_seq_pos_remove, [@info];
    }
  }
  $aln_obj = $aln_obj->remove_gaps('-', 1);
  ($left_tir_start, $right_tir_start) = get_columns($aln_obj, $tir_positions, 1);
  
  print "in con_fil sub after getCol: leftTIR: $left_tir_start\n";
  print $log_out "in con_fil sub after getCol: leftTIR: $left_tir_start\n";
  print "in con_fil sub after getCol: rightTIR: $right_tir_start\n";
  print $log_out "in con_fil sub after getCol: rightTIR: $right_tir_start\n";
  return ($left_tir_start, $right_tir_start, $gap_seq_pos_remove, $aln_obj, $tir_positions);
}

sub gap_filter {
  my $gap_seq_pos_remove = shift;              ## \@gap_seq_pos_remove
  my $aln_obj            = shift;
  my $left_tir_start     = shift;
  my $right_tir_start    = shift;
  my $aln_len            = $aln_obj->length;
  my $num_seqs = $aln_obj->num_sequences;

  my %trim_gap_seq_remove;
  for (my $i = 0 ; $i < $aln_len ; $i++) {
    if (($i >= $left_tir_start - 1 and $i <= $left_tir_start + 23) or ($i >= $right_tir_start - 25 and $i <= $right_tir_start - 1)) {

      my $gap_cols        = $aln_obj->gap_col_matrix();
      my @gap_col_array   = @{$gap_cols};
      my $gap_col_hashref = $gap_col_array[$i];
      my %gap_col_hash    = %{$gap_col_hashref};
      my $base_count      = 0;
      my $total_count     = 0;
      foreach my $key (keys %gap_col_hash) {
        if ($gap_col_hash{$key} != 1) {
          $base_count++;
        }
        $total_count++;
      }
      my $present_fraction = $base_count / $total_count;
      ## if half or more of the seqs in this col are gaps, mark seqs with a gaps
      ## to be be removed
      if ($present_fraction < 0.5) {
        foreach my $key (keys %gap_col_hash) {
          if ($gap_col_hash{$key} != 1) {
            my $seq_obj = $aln_obj->get_seq_by_id($key);
            $trim_gap_seq_remove{$key} = $seq_obj;
            my @info = ($key, $i + 1, $seq_obj);
            push @$gap_seq_pos_remove, [@info];
          }
        }
      }
      ## if more than half of the seqs are not gaps, but one seq has a gap,
      ## mark it for removal
      else {
        foreach my $key (keys %gap_col_hash) {
          if ($gap_col_hash{$key} == 1) {
            my $seq_obj = $aln_obj->get_seq_by_id($key);
            my $seq     = $seq_obj->seq();
            my $seq_pos = $seq_obj->location_from_column($i + 1);
            $trim_gap_seq_remove{$key} = $seq_obj;
            my @info = ($key, $i + 1, $seq_obj);
            push @$gap_seq_pos_remove, [@info];
          }
        }
      }
    }
  }
  
  my $remove_count = 0;
  foreach my $key (keys %trim_gap_seq_remove) {
      $remove_count++;
  }
  
  if (($num_seqs - $remove_count) > 5) {
      foreach my $key (keys %trim_gap_seq_remove) {
        my $seq_obj = $trim_gap_seq_remove{$key};
        if (defined $seq_obj) {
            $aln_obj->remove_seq($seq_obj);
        }
      }
  }
  else {
      undef $gap_seq_pos_remove;
  }
  return ($gap_seq_pos_remove, $aln_obj);
}

sub get_tir_nt_starts {
  my $aln_obj         = shift;
  my $tir_positions   = shift;    # \%tir_positions
  my $left_tir_start  = shift;
  my $right_tir_start = shift;
  my $remove_these;
  print "get_tir_nt_starts(top): lts:$left_tir_start rts:$right_tir_start\n";
  print $log_out "get_tir_nt_starts(top): lts:$left_tir_start rts:$right_tir_start\n";
  foreach my $seq_obj ($aln_obj->each_seq()) {
    my $seq_name            = $seq_obj->id();
    #my $seq                 = $seq_obj->seq();
    #if ($seq =~ /^-+$/g){
    #  #print "$seq_name: no seq found\n";
    #  #print $log_out "$seq_name: no seq found\n";
    #}
    if (defined $seq_obj->location_from_column($left_tir_start)){
      my $left_tir_start_obj  = $seq_obj->location_from_column($left_tir_start);
      my $left_tir_start_pos  = $left_tir_start_obj->start();
      if ($left_tir_start_pos > 0){
        $tir_positions->{$seq_name}{'left_tir_start'}  = $left_tir_start_pos;
      }
      else{
        delete  $tir_positions->{$seq_name};
        $$remove_these{$seq_name}=$seq_obj;
        next;
      }
    }
    else{
        delete  $tir_positions->{$seq_name};
        $$remove_these{$seq_name}=$seq_obj;
        next;
    }

    if (defined $seq_obj->location_from_column($right_tir_start)){
      my $right_tir_start_obj = $seq_obj->location_from_column($right_tir_start);
      my $right_tir_start_pos = $right_tir_start_obj->start();
      if ($right_tir_start_pos > 0){
        $tir_positions->{$seq_name}{'right_tir_start'}  = $right_tir_start_pos;
      }
      else{
        delete  $tir_positions->{$seq_name};
        $$remove_these{$seq_name}=$seq_obj;
        next;
      }
    }
    else{
        delete  $tir_positions->{$seq_name};
        $$remove_these{$seq_name}=$seq_obj;
        next;
    }
  }
  return ($tir_positions,$remove_these);
}

sub get_tir_nt_positions {
  ##($ref2_tir_positions)= get_tir_nt_positions ($aln_obj,\%tir_positions, $sorted_hitcolumn_keys[0],$sorted_hit_len_keys[0],$sorted_querycolumn_keys[0],$sorted_query_len_keys[0]);
  my $aln_obj         = shift;
  my $tir_positions   = shift;    # \%tir_positions
  my $left_tir_start  = shift;    # $sorted_hitcolumn_keys[0]
  my $left_tir_len    = shift;    # $sorted_hit_len_keys[0]
  my $right_tir_start = shift;    # $sorted_querycolumn_keys[0]
  my $right_tir_len   = shift;    # $sorted_query_len_keys[0]
  my $remove_these;
  print "in get_tir_nt_positions: lts=$left_tir_start,rts=$right_tir_start\n";
  print $log_out "in get_tir_nt_positions: lts=$left_tir_start,rts=$right_tir_start\n";
  foreach my $seq_obj ($trimmed_aln_obj->each_seq()) {
    my $seq_name           = $seq_obj->id();
    if (defined $seq_obj->location_from_column($left_tir_start)){
     my $left_tir_start_obj = $seq_obj->location_from_column($left_tir_start);
     my $left_tir_start_pos = $left_tir_start_obj->start();
     my $left_tir_end_pos   = $left_tir_start_pos + $left_tir_len - 1;
     if ($left_tir_start_pos > 0 and $left_tir_end_pos){
       $tir_positions->{$seq_name}{'left_tir_start'}  = $left_tir_start_pos;
       $tir_positions->{$seq_name}{'left_tir_end'}    = $left_tir_end_pos;
     }
     else{
        delete  $tir_positions->{$seq_name};
        $$remove_these{$seq_name}=$seq_obj;
        next;
      }
    }
    else{
      delete  $tir_positions->{$seq_name};
      $$remove_these{$seq_name}=$seq_obj;
      next;
    }
    if (defined $seq_obj->location_from_column($right_tir_start)){
      my $right_tir_start_obj = $seq_obj->location_from_column($right_tir_start);
      my $right_tir_start_pos = $right_tir_start_obj->start();
      my $right_tir_end_pos = $right_tir_start_pos - ($right_tir_len - 1);
      if ($right_tir_start_pos > 0 and $right_tir_end_pos > 0){
        $tir_positions->{$seq_name}{'right_tir_start'} = $right_tir_start_pos;
        $tir_positions->{$seq_name}{'right_tir_end'}   = $right_tir_end_pos;
      }
      else{
        delete  $tir_positions->{$seq_name};
        $$remove_these{$seq_name}=$seq_obj;
        next;
      }
    }
    else{
      delete  $tir_positions->{$seq_name};
      $$remove_these{$seq_name}=$seq_obj;
      next;
    }
  }
  return ($tir_positions, $remove_these);
}

sub error_out {
  my $bad_out_path = shift;    ## error file path
  my $message      = shift;
  my $exit         = shift;    ## 1=no exit; 0 or undef=exit
  open(my $bad_out, ">", $bad_out_path);
  print $bad_out "$message\n";
  print "ERROR:$message\n";
  print $log_out "ERROR:$message\n";
  #warn "ERROR:$message\n";
  close($bad_out);
  if (!defined $all) {
    print "Cleaning up files\n";
    print $log_out "Cleaning up files\n";
    clean_files($out_path);
  }
  if (!$exit or $exit == 0) {
    exit 0;
  }
}

sub print_fasta {
  my $filename = shift;
  my $aln_obj  = shift;
  my $seqIO_out_obj = Bio::SeqIO->new(-format => 'fasta', -file => ">$filename");
  foreach my $seq_obj ($aln_obj->each_seq) {
    my $seq = $seq_obj->seq;
    $seq =~ s/-//g;
    $seq_obj->seq($seq);
    $seqIO_out_obj->write_seq($seq_obj);
  }
}


sub get_percentID_perCol {
## @percent_id = get_percentID_perCol(file.aln,outfile);
  my $msa = shift;
  my $out = shift;
  open MSA,    $msa        or die "Can't open $msa\n";
  <MSA>;    #throw out first header
  my $seq;
  my $seq_count = 0;
  my @seq_order;
  my %seqs;
  while (my $line = <MSA>) {
    chomp $line;
    if ($line =~ /^>(\S+)/) {
      push @seq_order , $1;
      $seqs{$seq_count}= $seq;
      $seq_count++;
      $seq = '';
    }
    else {
      $seq .= $line;
    }
  }
## for last seq
  $seqs{$seq_count}= $seq;
  $seq_count++;
##
  close MSA;


  my @percent_id;
  #my $first_seq = `head -1 $msa.mod`;
  my $first_seq_len = length $seqs{0};
  if ($first_seq_len < 5){
   die "Error retrieving first line of $msa.mod\n";
  }
  open OUT, ">$out" or die "Can't opne $out $!\n";
  for (my $i = 0 ; $i < $first_seq_len ; $i++) {
    my $total_nt;
    my %nt_count;
    my @col;
    foreach my $seq_num (sort {$a <=> $b} keys %seqs){
      push @col , substr $seqs{$seq_num} , $i , 1;
    }
    foreach my $nt (@col) {
      $nt_count{total}++;
      if ($nt =~ /-/) {
        $nt_count{dash}++;
        next;
      }
      ## only increment if the character is not a '-'
      ## Ns will be included in the total
      #$nt_count{total}++;
      #next if $nt =~ /N/i;
      $nt_count{each}{$nt}++;
    }
    if ((scalar keys %{ $nt_count{each} }) == 0) {
      push @percent_id, [ $i+1, 0, 0 ];
      print OUT join("\t", $i+1, 0, 0), "\n";
    }
    else {
      my $most_freq_nt =
        (sort { $nt_count{each}{$b} <=> $nt_count{each}{$a} }
          keys %{ $nt_count{each} })[0];
      if ($most_freq_nt =~ /N/i) {
        push @percent_id, [ $i+1, 0, 0 ];
        print OUT join("\t", $i+1, 0, 0), "\n";
      }
      my $count      = $nt_count{each}{$most_freq_nt};
      my $col_total  = $nt_count{total};
      my $dash_total = exists $nt_count{dash} ? $nt_count{dash} : 0;
      ## rigth now, col_total is eq to seq_count
      my $pid = ($count / $col_total) * 100;

      #push @percent_id, ($count/$col_total);
      push @percent_id,
        [ $i+1, $pid, (($col_total - $dash_total) / $seq_count) ];

      # print out this info to a file
      print OUT
        join("\t", $i+1, $pid, (($col_total - $dash_total) / $seq_count)),
        "\n";
    }
  }
  return @percent_id;
}

sub remove_most {

#my ($left_tir_start1,$right_tir_start1,$tmp_aln_obj,$ref2tp,$ref2gsr, $ref2gspr) = remove_most ($full_aln_obj,\%tir_positions, \@full_id_array);
  my $org_aln_obj   = shift;
  my $tir_positions = shift;
  my $id_array      = shift;
  my $aln_obj       = $org_aln_obj;
  my $aln_len       = $aln_obj->length;
  my $ref2gsr;    ## %gap_seq_remove
  my $try = shift;
  
  for (my $i = 1 ; $i < $aln_len ; $i++) {
    my $id_row_ref      = $$id_array[$i];
    my $gap_col_hashref = $gap_cols->[$i];
    my %gap_col_hash    = %{$gap_col_hashref};
    ## remove entire seq if %gaps in column is >= 80
    ## and (%id <= 50 or %nt_aligned <= 10%)
    if ($gap_id_array[$i] >= 80.0 and ($id_row_ref->[1] <= 50 or $id_row_ref->[2] <= .1)) {
      foreach my $key (keys %gap_col_hash) {
        ## value of 1 means gap
        ## if there is no gap at this position where there are many gaps in the column
        ## we will make a note to remove this seq
        if ($gap_col_hash{$key} != 1) {
          my $seq_obj = $aln_obj->get_seq_by_id($key);
          my @info = ($key, $i + 1, $seq_obj);
          push @$ref2gspr, [@info];
          $$ref2gsr{$key} = $seq_obj;
        }
      }
    }
    elsif ($i > ($flank/2) and $i < ($aln_len - ($flank/2)) and $gap_id_array[$i] < 30.0) {
        foreach my $key (keys %gap_col_hash) {
            if ($gap_col_hash{$key} == 1) {
                my $seq_obj = $aln_obj->get_seq_by_id($key);
                my @info = ($key, $i + 1, $seq_obj);
                push @$ref2gspr, [@info];
                $$ref2gsr{$key} = $seq_obj;
            }
        }
    }
  }
  ## removes any sequences found in last step
  foreach my $key (keys %$ref2gsr) {
    my $seq_obj = $$ref2gsr{$key};
    if (defined $seq_obj) {
        $aln_obj->remove_seq($seq_obj);
    }
  }
  ## removes any columns that are now all gaps
  $aln_obj = $aln_obj->remove_gaps('-', 1);

  my $test_aln_out = File::Spec->catpath($volume, $out_path, $filename . ".removeMost_trim_0");
  $out = Bio::AlignIO->new(
    -file             => ">$test_aln_out",
    -format           => 'fasta',
    -displayname_flat => 0
 );
  $out->write_aln($aln_obj);

  my ($left_tir_start, $right_tir_start);
  ($left_tir_start, $right_tir_start, $ref2gspr, $aln_obj, $tir_positions) = consensus_filter($ref2gspr, $aln_obj, 0, 0, $tir_positions, "other", $try);

  #@gap_seq_pos_remove = @{$ref2gspr};

  print "remove most after cons_filter filename.trim: $left_tir_start, $right_tir_start\n";
  print $log_out "remove most after cons_filter filename.trim: $left_tir_start, $right_tir_start\n";

#my ($left_tir_start1,$right_tir_start1,$tmp_aln_obj,$ref2tp,$ref2gsr, $ref2gspr) = remove_most ($full_aln_obj,\%tir_positions, \@full_id_array);
  return ($left_tir_start, $right_tir_start, $aln_obj, $tir_positions, $ref2gsr, $ref2gspr);

} ## end remove_most

sub remove_least {

#my ($left_tir_start1,$right_tir_start1,$tmp_aln_obj,$ref2tp,$ref2gsr, $ref2gspr) = remove_most ($full_aln_obj,\%tir_positions, \@full_id_array);
  my $org_aln_obj   = shift;
  my $tir_positions = shift;
  my $id_array      = shift;
  my $aln_obj       = $org_aln_obj;
  my $aln_len       = $aln_obj->length;
  my $ref2gsr;    ## %gap_seq_remove
  
  if (defined $protein) {
      my ($left_tir_start, $right_tir_start);
      ($left_tir_start, $right_tir_start, $ref2gspr, $aln_obj, $ref2tp) = consensus_filter($ref2gspr, $aln_obj, 0, 0, $tir_positions);
      print "after protein cons_filter filename.trim: $left_tir_start, $right_tir_start\n";
      print $log_out "after protein cons_filter filename.trim: $left_tir_start, $right_tir_start\n";
      
      return ($left_tir_start, $right_tir_start, $aln_obj, $tir_positions, $ref2gsr, $ref2gspr);
      print "Shouldn't print this!\n";
  }

  foreach my $arrayref (@$id_array) {
    my $col_id = $$arrayref[0];
    my $pos_id = $$arrayref[1];
    if ($first_col_80 == 1 and $pos_id >= 80) {
      $first_col_80 = $col_id;
    }
    elsif ($pos_id >= 80) {
      $last_col_80 = $col_id;
    }
  }

  print "first_col_80:$first_col_80 last_col_80:$last_col_80\n";
  print $log_out "first_col_80:$first_col_80 last_col_80:$last_col_80\n";
  for (my $i = $first_col_80 ; $i < $last_col_80 ; $i++) {
    next if $i == $first_col_80 + ($aln_len * .4) or $i < $last_col_80 - ($aln_len * .4);
    my $id_row_ref      = $$id_array[$i];
    my $gap_col_hashref = $gap_cols->[$i];
    my %gap_col_hash    = %{$gap_col_hashref};
    ## remove entire seq if %gaps in column is >= 80
    ## and (%id <= 50 or %nt_aligned <= 10%)
    if ($gap_id_array[$i] >= 80.0
      and ($id_row_ref->[1] <= 50 or $id_row_ref->[2] <= .1))
    {
      foreach my $key (keys %gap_col_hash) {
        ## value of 1 means gap
        ## if there is no gap at this position where there are may gaps in the column
        ## we will make a note to remove this seq
        if ($gap_col_hash{$key} != 1) {
          my $seq_obj = $aln_obj->get_seq_by_id($key);
          my @info = ($key, $i + 1, $seq_obj);
          push @$ref2gspr, [@info];
          $$ref2gsr{$key} = $seq_obj;
        }
      }
    }
  }
## removes any sequences found in last step
  foreach my $key (keys %$ref2gsr) {
    my $seq_obj = $$ref2gsr{$key};
    if (defined $seq_obj) {
        $aln_obj->remove_seq($seq_obj);
    }
  }
## removes any columns that are now all gaps
  $aln_obj = $aln_obj->remove_gaps('-', 1);

  my $test_aln_out = File::Spec->catpath($volume, $out_path,
    $filename . ".removeLeast_trim_0");
  my $out = Bio::AlignIO->new(
    -file             => ">$test_aln_out",
    -format           => 'fasta',
    -displayname_flat => 0
 );
  $out->write_aln($aln_obj);

  my ($left_tir_start, $right_tir_start);
  ($left_tir_start, $right_tir_start, $ref2gspr, $aln_obj, $ref2tp) = consensus_filter($ref2gspr, $aln_obj, 0, 0, $tir_positions, "other", $try);

  #@gap_seq_pos_remove = @{$ref2gspr};

  print "remove least after cons_filter filename.trim: $left_tir_start, $right_tir_start\n";
  print $log_out "remove least after cons_filter filename.trim: $left_tir_start, $right_tir_start\n";

  return ($left_tir_start, $right_tir_start, $aln_obj, $tir_positions,
    $ref2gsr, $ref2gspr);

} ## end remove_least

sub get_org_aln {
  my $infile = shift;
  #create input object and an alignment object from it
  warn "Reading in original alignment file\n";
  my $in_obj = Bio::AlignIO->new(-file => $infile, -format => 'fasta' , -alphabet => 'dna');
  my $aln_obj;
  unless ($aln_obj = $in_obj->next_aln()) {
    die "Cannot find any alignment in file $infile\n";
  }
  
  return ($aln_obj);
}

sub tir_mismatch_filter {
    my $aln_obj = shift;
    my $left_tir_start = shift;
    my $left_tir_end = shift;
    my $right_tir_start = shift;
    my $right_tir_end = shift;
    my $round = shift;
    my $recheck_hash = shift;

    my @to_remove;
    my $max_mismatch = 0;
    my $max_half = 0;
    if (!defined $round or $round == 0) {
        $round = 0;
        $left_tir_end = $left_tir_start+20;
        $right_tir_end = $right_tir_start-20;
        $max_mismatch = 5;
        $max_half = $max_mismatch/2;
    }
    if ($round == 1) {
        $max_mismatch = int(((($left_tir_end-$left_tir_start) + ($right_tir_start-$right_tir_end))*0.15));
        $max_half = int($max_mismatch/2)+1;
    }
    if ($round == 2) {
        $max_mismatch = int(((($left_tir_end-$left_tir_start) + ($right_tir_start-$right_tir_end))*0.1))+1;
        $max_half = int($max_mismatch/2)+1;
    }
    if ($max_half <= 1) {
        $max_half = 2;
    }
    print "Max mismatch = $max_mismatch  Max half: $max_half\n";
    print $log_out "Max mismatch = $max_mismatch  Max half: $max_half\n";
    #check for too many mismatches in putative TIRSs against the consensus sequence
    my $aln_consensus = $aln_obj->consensus_string(80);
    #print "Consensus for mismatches:\n$aln_consensus\n\n";
    #print $log_out "Consensus for mismatches:\n$aln_consensus\n\n";

    if ($recheck_hash) {
        foreach my $key (keys %$recheck_hash) {
            #print "Adding back: $key\n";
            #print $log_out "Adding back: $key\n";
            my $seq_obj = $$recheck_hash{$key};
            $aln_obj->add_seq($seq_obj);
        }
    }
    #print "left start: $left_tir_start  left end: $left_tir_end  right end: $right_tir_end  right start: $right_tir_start\n";
    #print $log_out "left start: $left_tir_start  left end: $left_tir_end  right end: $right_tir_end  right start: $right_tir_start\n";
    foreach my $seq_obj ($aln_obj->each_seq()) {
        my $mismatches = 0;
        my $left_mis = 0;
        my $right_mis = 0;
        my $seq_name = $seq_obj->id();
        my $seq = $seq_obj->seq();
        #print "Seq for $seq_name:\n$seq\n";
        #print $log_out "Seq for $seq_name:\n$seq\n";
        for (my $i = 1; $i < length($aln_consensus); $i++) {
            my $pos_seq = substr($seq, $i-1, 1);
            my $consensus_pos = substr($aln_consensus, ($i-1), 1);
            #print "compare: $pos_seq to $consensus_pos\n";
            #print $log_out "compare: $pos_seq to $consensus_pos\n";
            if ($i >= $left_tir_start and $i <= $left_tir_end) {
                print "In left TIR, compare at $i: $pos_seq to $consensus_pos\n";
                print $log_out "In left TIR, compare at $i: $pos_seq to $consensus_pos\n";
                next if ($round == 0 or $round == 1) and $pos_seq eq "-";
                if ($round == 2 and $pos_seq eq "-") {
                    $mismatches++;
                    $left_mis++;
                    next;
                }
                next if $consensus_pos eq "?";
                if ($pos_seq ne $consensus_pos) {
                    $mismatches++;
                    $right_mis++;
                    if ($round == 0) {
                        if ($i <= ($left_tir_start + 5)) {
                            push @to_remove, $seq_name;
                            last;
                        }
                    }
                    else {
                        if ($i <= ($left_tir_start)) {
                            push @to_remove, $seq_name;
                            last;
                        }
                    }
                }
                else {
                    #print "Equal\n";
                    #print $log_out "Equal\n";
                }
            }

            elsif ($i >= $right_tir_end and $i <= $right_tir_start) {
                print "In right TIR, compare at $i: $pos_seq to $consensus_pos\n";
                print $log_out "In right TIR, compare at $i: $pos_seq to $consensus_pos\n";
                next if ($round == 0 or $round == 1) and $pos_seq eq "-";
                if ($round == 2 and $pos_seq eq "-") {
                    $mismatches++;
                    $right_mis++;
                    next;
                }
                next if $consensus_pos eq "?";
                if ($pos_seq ne $consensus_pos) {
                    $mismatches++;
                    $right_mis++;
                    if ($round == 0) {
                        if ($i >= ($right_tir_start - 5)) {
                            push @to_remove, $seq_name;
                            last;
                        }
                    }
                    else {
                        if ($i >= ($right_tir_start)) {
                            push @to_remove, $seq_name;
                            last;
                        }
                    }
                }
                else {
                    #print "Equal\n";
                    #print $log_out "Equal\n";
                }
            }
            else {
                #print "$i not in TIR sequences\n";
                #print $log_out "$i not in TIR sequences\n";
                next;
            }
            if ($mismatches >= $max_mismatch or $left_mis >= $max_half or $right_mis >= $max_half) {
                push @to_remove, $seq_name;
                last;
            }
        }
        print "Mismatches = $mismatches  $seq_name\n\n";
        print $log_out "Mismatches = $mismatches  $seq_name\n\n";
    }
    return(@to_remove);
}

sub match_tirs2 {
    my $self       = shift;    ## seq_obj
    my $input_path = shift;
    
    $self->throw("Need Bio::LocatableSeq argument")
        unless ref $self && $self->isa('Bio::LocatableSeq');
    if (!-f $input_path) {
        die "Supplied filepath is not valid";
    }
    my $seq = $self->seq();
    $seq =~ s/-//g;
    my $seq_len = length($seq);
    my $check_seq = $seq;
    $check_seq =~ tr/ATGCatgc/TACGtacg/;
    $check_seq = reverse($check_seq);
    my $half_seq = substr($seq, 0, $seq_len/2);
    my $seq_name = $self->id();
    my $fa_aln_obj = Bio::SearchIO->new(-format => 'fasta', -file => $input_path);
    my @result;
    #print "ggsearch results path: $input_path\n";
    #print $log_out "ggsearch results path: $input_path\n";
    
    #go through the FASTA input object . . down to the HSP
    while (my $result = $fa_aln_obj->next_result) {
        if ($result->num_hits == 0) {
            push @result, (0, [0]);
            last;
        }
        while (my $hit = $result->next_hit) {
            #initialize variables
            my $matches     = 0;
            my $mismatches = 0;
            my $first_three_count = 0;
            while (my $hsp = $hit->next_hsp) {
                
                #grab the query, hit, and homology strings
                my $homo_string = $hsp->homology_string;
                my $query_str   = $hsp->query_string;
                my $hit_str     = $hsp->hit_string;
                my $len_homo    = length($homo_string);
                my $len_query   = length($query_str);
                my $len_hit     = length($hit_str);
                
                
                my $match_query   = '';
                my $match_hit     = '';
                my $last_good     = 0;
                my $hit_pos;
                my $query_pos;
                
                for (my $count = 0 ; $count < 16 ; $count++) {
                    my $homo_char  = substr($homo_string, $count, 1);
                    my $query_char = substr($query_str,   $count, 1);
                    my $hit_char   = substr($hit_str,     $count, 1);
                    
                    if ($homo_char eq ":") {
                        $matches++;
                        $last_good = $count;
                        if ($count < 3 and $query_char ne "N" and $query_char ne "n" and $hit_char ne "N" and $hit_char ne "n") {
                            $first_three_count++;
                        }
                    }
                    else {
                        $mismatches++;
                    }
                }
                if ($first_three_count == 3) {
                    if ($mismatches/$matches <= 0.25) {
                        $match_hit = substr($hit_str, 0, $last_good + 1);
                        $match_query = substr($query_str, 0, $last_good + 1);
                        my $match_hit_len = $last_good + 1;
                        my $match_query_len = $last_good + 1;
                        $match_query =~ s/-//g;
                        $match_hit   =~ s/-//g;
                        
                        $hit_pos = index(uc($half_seq), uc($match_hit)) + 1;
                        my $initial_query_pos = index(uc($check_seq), uc($match_query));
                        $query_pos = $seq_len - $initial_query_pos;
                        
                        $match_query =~ tr/ATGCatgc/TACGtacg/;
                        $match_query = reverse($match_query);
                        
                        my @match = (
                        $seq_name,
                        {
                            "hit"   => [ $hit_pos,   $match_hit_len ],
                            "query" => [ $query_pos, $match_query_len ],
                            "seq" => [$match_hit, $match_query]
                        }
                        );
                        ## add a catch for $hit_pos or $query_pos == 0 [index returns -1 when $match_hit was not found in $seq]
                        if ($hit_pos == 0 or $query_pos == 0) {
                            push @result, [0, [0]];
                        }
                        else {
                            push @result, (1, [@match]);
                        }
                    }
                    else {
                        push @result, [0, [0]];
                    }
                }
                else {
                    push @result, [0, [0]];
                }
            }
        }
    }
    return (@result);
}
