#!/usr/bin/env perl

=head1 NAME

activeTE_msa - evaluate MSAs for characteristics of active TEs

=head1 SYNOPSIS

activeTE_msa family.aln

=head1 DESCRIPTION

=head1 Command line arguments

 -flank   [integer] Length of flanking sequence (default 100)
 -all     Process all

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

#set default flank length and reset if option is given

my $flank = 100;
my $trimal =
  'trimal';   # specify the full path to this if it can't be found automatically
my $all;
my $PROGRAM_NAME = "activeTE";

GetOptions(
  'flank:i'  => \$flank,
  'all'      => \$all,
  'trimal:s' => \$trimal,
);

if ( defined $all ) {
  print "Intermediate files will be kept\n";
}
else {
  print "Intermediate files will be removed\n";
}

# make sure only one valid input file is specified
if ( @ARGV != 1 ) {
  print "Please specify only one input file\n";
  help();
}
my $infile = shift @ARGV;

if ( !-f $infile ) {
  print "Supplied filepath is not valid: $!\n";
  help();
}
if ( $flank <= 25 ) {
  print "Please use alignments with flanks  >25 nucleotides long\n";
  help();
}

sub help {
  print "

usage:

activeTE_msa.pl -a -f <int> <multiple alignment file>

-a keep all intermediate files  DEFAULT = remove
-f length of sequence flanks  DEFAULT = 100

The test alignments all have f = 100.

";
  exit 1;
}

#breakdown directory and filename, create output directory
my ( $volume, $in_dir, $filename ) = File::Spec->splitpath($infile);

#my ($fname_fin) = split( '\\.', $filename );    # first bit after first '.' is the name
my @fname_fin =
  split( '\\.', $filename );    # first bit after first '.' is the name
pop @fname_fin;
my $fname_fin    = join( '.', @fname_fin );
my $out_dir_name = "activeTE_out_" . $fname_fin;
my $out_path     = File::Spec->catdir( $in_dir, $out_dir_name );
if ( !-d $out_path ) {
  mkdir($out_path)
    or die "Can't create dir:$out_path $!\n";    # use built-in Perl mkdir
}

my %element_info;    # this hash is for storing element information
my %tir_positions;

#generate trimal alignment gap summary file, using sed to remove the header infmal -in $infile -sgco

open( my $trimal_run => "$trimal -in $infile -sgc | sed -n '4,\$p' |" )
  || die "Cannot run trimal!";

# parse the trimal output file with the gap %  of each position in the full alignment
my @gap_id_array;
while ( my $line = <$trimal_run> ) {
  chomp $line;       # probably not necessary
  if ( $line =~ /\S+\s+(\S+)\s+\S+/ ) {
    push @gap_id_array, $1;
  }
}
close($trimal_run);

#initialize a global array to store results later

my @good_aln = ();
#my ( $initial_left_TSD, $initial_right_TSD );
## replacing white space with underscores in MSA
my @in_array;
open( my $in => "<$infile" );
while ( my $line = <$in> ) {
  $line =~ s/ /_/g;

  # did you want to skip lines that don't have any whitespace?
  push @in_array, $line;
}
close($in);
open( my $fix_out => ">$infile" );
print $fix_out @in_array;
close($fix_out);
@in_array = ();

# the above could just be re-written
# `perl -i -p -e 's/ /_/g;' $infile`;

my $full_aln_obj = get_org_aln ($infile);
my $full_aln_len = $full_aln_obj->length();
my $gap_cols = $full_aln_obj->gap_col_matrix();
my $full_aln_num_seqs = $full_aln_obj->num_sequences;
if ($full_aln_num_seqs < 2){
  my $bad_out_path =
    File::Spec->catpath( $volume, $out_path, $filename . ".bad" );
  error_out( $bad_out_path,
    "$filename\thas only $full_aln_num_seqs sequence(s)\n" );
}

#calculate the % nucleotide identity and the fraction of copies with sequence at each position the full alignment and print to file

my $full_id_out_path =
  File::Spec->catpath( $volume, $out_path, $filename . ".full_id" );

my @full_id_array;

my $first_col_80 = 1;
my $last_col_80  = $full_aln_len;
if ( !-e $full_id_out_path ) {
  print "Starting Full ID calculation\n";
  @full_id_array = get_percentID_perCol( $infile, $full_id_out_path );
}
else {
  open( my $id_out, '<', "$full_id_out_path" );
  ##head 169.msa.full_id
  ##1       0       0.000893655049151028
  ##2       0       0.000893655049151028
  while ( my $line = <$id_out> ) {
    chomp $line;
    my ( $i, $pos_id, $pos_present ) = split /\t/, $line;
    push @full_id_array, [ $i, $pos_id, $pos_present ];
  }
}
## remove as many seqs as possible for the cleanest aln
my ( $left_tir_start1, $right_tir_start1, $tmp_aln_obj, $ref2tp, $ref2gsr,
  $ref2gspr )
  = remove_most( $full_aln_obj, \%tir_positions, \@full_id_array );
%tir_positions = %$ref2tp;
my %gap_seq_remove     = %$ref2gsr;
my @gap_seq_pos_remove = @$ref2gspr;

## if this removes too many, remove a few as possible
if ( $left_tir_start1 == 0 or $right_tir_start1 == 0 ) {
  my $aln_obj = get_org_aln ($infile);
  print "run less stringent intial filtering\n";
  (
    $left_tir_start1, $right_tir_start1, $tmp_aln_obj, $ref2tp, $ref2gsr,
    $ref2gspr
  ) = remove_least( $aln_obj, \%tir_positions, \@full_id_array );
  %tir_positions      = %$ref2tp;
  %gap_seq_remove     = %$ref2gsr;
  @gap_seq_pos_remove = @$ref2gspr;
}
##overwrite full_aln with tmp_aln
$full_aln_obj = $tmp_aln_obj;

## printing out the new cleaned up alignments
## if more than 5 sequences, make a new alignment with muscle
my $trim_aln_out =
  File::Spec->catpath( $volume, $out_path, $filename . ".trim" );
my $out = Bio::AlignIO->new(
  -file             => ">$trim_aln_out",
  -format           => 'fasta',
  -displayname_flat => 0
);
$out->write_aln($full_aln_obj);

## removes any columns that are now all gaps
$full_aln_obj = $full_aln_obj->remove_gaps( '-', 1 );
( $left_tir_start1, $right_tir_start1 ) =
  get_columns( $full_aln_obj, \%tir_positions, 1 );

#print "after removing gap-onlycolumsn from .trim get_columns (filename.trim_3): $left_tir_start1, $right_tir_start1\n";

#$trim_aln_out =
#  File::Spec->catpath( $volume, $out_path, $filename . ".trim_3" );
#$out = Bio::AlignIO->new( -file => ">$trim_aln_out", -format => 'fasta' , -displayname_flat => 0);
#$out->write_aln($full_aln_obj);

my $trimmed_aln_obj = $full_aln_obj;

if ( !defined $left_tir_start1
  or !defined $right_tir_start1
  or $left_tir_start1 == 0
  or $right_tir_start1 == 0 )
{

  #open a file to store info on why analysis of an element was aborted
  my $bad_out_path =
    File::Spec->catpath( $volume, $out_path, $filename . ".bad" );
  error_out( $bad_out_path,
    "$filename\tTIRs not found in MSA by initial search" );
}

my ( $org_right_tir_start1, $org_left_tir_start1 ) =
  ( $right_tir_start1, $left_tir_start1 );
( $right_tir_start1, $left_tir_start1 ) =
  adjust_tir_starts( $trimmed_aln_obj, $right_tir_start1, $left_tir_start1 );

my $left_tir_adjusted = $org_left_tir_start1 - $left_tir_start1;
my $right_tir_adjusted = $right_tir_start1 - $org_right_tir_start1;
#check whether any sequences have gaps within 50bp of start of putative TIRSs and remove them
my ( $left_tir_start, $right_tir_start, $ref2array, $ref2hash );
( $left_tir_start, $right_tir_start, $ref2array, $trimmed_aln_obj, $ref2hash ) =
  consensus_filter(
  \@gap_seq_pos_remove, $trimmed_aln_obj, $left_tir_start1,
  $right_tir_start1,    \%tir_positions
  );
@gap_seq_pos_remove = @$ref2array;
%tir_positions      = %$ref2hash;

print
  "new tir starts after consensus filter: $left_tir_start, $right_tir_start\n";
( $ref2array, $trimmed_aln_obj ) =
  gap_filter( \@gap_seq_pos_remove, $trimmed_aln_obj, $left_tir_start,
  $right_tir_start );
@gap_seq_pos_remove = @$ref2array;
$trimmed_aln_obj = $trimmed_aln_obj->remove_gaps( '-', 1 );
##Sofia added 06192013
print "before get_Col: Trim2 Left TIR start column: $left_tir_start1\n";
( $left_tir_start1, $right_tir_start1 ) =
  get_columns( $trimmed_aln_obj, \%tir_positions, 1 );

my $test_len = $trimmed_aln_obj->length();
if ( $test_len == 0 ) {

  #open a file to store info on why analysis of an element was aborted
  my $bad_out_path =
    File::Spec->catpath( $volume, $out_path, $filename . ".bad" );
  error_out( $bad_out_path,
    "$filename\tTIRs not found in MSA. All sequences filtered out." );
}

my $trim_aln_out2 =
  File::Spec->catpath( $volume, $out_path, $filename . ".trim2" );
my $out2 = Bio::AlignIO->new(
  -file             => ">$trim_aln_out2",
  -format           => 'fasta',
  -displayname_flat => 0
);
$out2->write_aln($trimmed_aln_obj);
print "after get_Col: Trim2 Left TIR start column: $left_tir_start1\n";
print "Trim2 Right TIR start column: $right_tir_start1 (filename.trim2)\n";

#Store the column positions of the potential TIRs in the trimmed alignment and get the columns of them in the original alignminent to remove sequences that cause gaps in the tirs
my %trim_left_pos_hash;
my %trim_right_pos_hash;
my %left_tir_start_check_counts;
my %right_tir_start_check_counts;

foreach my $seq_obj ( $trimmed_aln_obj->each_seq() ) {
  my $left_res_pos_obj  = $seq_obj->location_from_column($left_tir_start1);
  my $right_res_pos_obj = $seq_obj->location_from_column($right_tir_start1);
  my $seq               = $seq_obj->seq();
  my $seq_len           = length($seq);
  my $seq_name          = $seq_obj->id();
  my $left_res_pos      = $left_res_pos_obj->start();
  my $left_res          = substr( $seq, $left_res_pos, 1 );
  my $right_res_pos     = $right_res_pos_obj->start();
  my $right_res         = substr( $seq, $right_res_pos, 1 );

  #my $right_flank_pos   = $seq_len - $right_res_pos;
  $trim_left_pos_hash{$seq_name}  = $left_res_pos;
  $trim_right_pos_hash{$seq_name} = $right_res_pos;
  $left_tir_start_check_counts{$left_res_pos}++;
  $right_tir_start_check_counts{$right_res_pos}++;

#print "TIR starts for $seq_name:  Left - $left_res_pos  Right - $right_res_pos\n";
}

#sort the left and right TIR starts by largest count to smallest
my @sorted_left_tir_start_keys =
  sort { $left_tir_start_check_counts{$b} <=> $left_tir_start_check_counts{$a} }
  keys(%left_tir_start_check_counts);
my @sorted_right_tir_start_keys = sort {
  $right_tir_start_check_counts{$b} <=> $right_tir_start_check_counts{$a}
} keys(%right_tir_start_check_counts);

print
"After removing copies with problems in TIRs these are the most common nucleotide positions of the TIRs:\n\tLeft TIR  start: $sorted_left_tir_start_keys[0]\tRight TIR start: $sorted_right_tir_start_keys[0]\n";
print "This is after the .trim to .trim2 transition\n";

#check if TIRs start in flanks, indicating the flanks are highly similar
my $left_flank_catch  = 0;
my $right_flank_catch = 0;

if ( !defined $sorted_left_tir_start_keys[0] ) {
  if ( !defined $sorted_right_tir_start_keys[0] ) {

    #open a file to store info on why analysis of an element was aborted
    my $bad_out_path =
      File::Spec->catpath( $volume, $out_path, $filename . ".bad" );
    error_out( $bad_out_path,
      "$filename\tTIRs not found or both flanks similar" );
  }
  $left_flank_catch++;

  #open a file to store info on why analysis of an element was aborted
  my $bad_out_path =
    File::Spec->catpath( $volume, $out_path, $filename . ".bad" );
  error_out( $bad_out_path, "$filename\tLeft flank similar" );
}
elsif ( !defined $sorted_right_tir_start_keys[0] ) {

  #open a file to store info on why analysis of an element was aborted
  my $bad_out_path =
    File::Spec->catpath( $volume, $out_path, $filename . ".bad" );
  error_out( $bad_out_path, "$filename\tRight flank similar" );
}

if ( $sorted_left_tir_start_keys[0] <= $flank - 25 ) {
  $left_flank_catch++;
}
if ( $sorted_right_tir_start_keys[0] <= $flank - 25 ) {
  $right_flank_catch++;
}

if ( $left_flank_catch != 0 ) {
  if ( $right_flank_catch != 0 ) {

    #open a file to store info on why analysis of an element was aborted
    my $bad_out_path =
      File::Spec->catpath( $volume, $out_path, $filename . ".bad" );
    error_out( $bad_out_path, "$filename\tBoth flanks similar" );
  }
  else {

    #open a file to store info on why analysis of an element was aborted
    my $bad_out_path =
      File::Spec->catpath( $volume, $out_path, $filename . ".bad" );
    error_out( $bad_out_path, "$filename\tLeft flank similar" );
  }
}
elsif ( $right_flank_catch != 0 ) {

  #open a file to store info on why analysis of an element was aborted
  my $bad_out_path =
    File::Spec->catpath( $volume, $out_path, $filename . ".bad" );
  error_out( $bad_out_path, "$filename\tRight flank similar" );
}

my @bad_remove;
my @bad_aln;

my $element_len;
my $ele_half_len;

#grab all sequences from trimmed alignment
foreach my $seq_obj ( $trimmed_aln_obj->each_seq() ) {

  #grab sequence name and shorten it for use in some output filenames
  my $seq_name    = $seq_obj->id();
  my @seqn_part   = split( ":", $seq_name );
  my $fname_start = $seqn_part[0];

  #grab sequence & strip hyphens
  my $seq = $seq_obj->seq();
  $seq =~ s/-//g;

  #initialize variables
  my $sub_seq;
  my $first;
  my $last;
  my $first_path = File::Spec->catpath( $volume, $out_path, "first_half.txt" );
  my $last_path  = File::Spec->catpath( $volume, $out_path, "last_half.txt" );

  #get end  sequences
  $element_len =
    $trim_right_pos_hash{$seq_name} - $trim_left_pos_hash{$seq_name} + 1;
  $ele_half_len = int( $element_len / 2 );

  $first = substr( $seq, $trim_left_pos_hash{$seq_name} - 1, $ele_half_len );
  $last = substr( $seq, $trim_right_pos_hash{$seq_name} - $ele_half_len,
    $ele_half_len );

  #save the two ends as files to use as inputs for a ggsearch search
  open( my $first_out, ">", $first_path )
    or die "Error creating $first_path. $!\n";
  open( my $last_out, ">", $last_path )
    or die "Error creating $last_path. $!\n";
  print $first_out ">first\n$first\n";
  print $last_out ">last\n$last\n";
  close($first_out);
  close($last_out);

  #create fasta output filename then call ggsearch
  my $out_opt =
    File::Spec->catpath( $volume, $out_path, $fname_start . ".ggsearch.out" );
  system("ggsearch36 -n -i -T 1 -d 1 $last_path $first_path > $out_opt");

  my @tir_match_result = match_tirs( $seq_obj, $out_opt, 1 );

  if ( $tir_match_result[0] == 1 ) {
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

if ( $good_aln_len == 0 ) {
  print "No TIRs found from MSA: First round\n";
}
else {
  print "Found TIRs in MSA: 1st round\n";

  #go through each entry of @good_aln array
  foreach my $row_ref (@good_aln) {
    my $query_aln_pos;
    my $hit_aln_pos;
    @entry = @{$row_ref};

#get column position in alignment of left (hit) TIR, using sequence name and index position of TIR in full sequence. Increment the count of that column position in %hit_column_counts. Also, increment the count of TIR length in %hit_match_len
    $hit_aln_pos =
      $trimmed_aln_obj->column_from_residue_number( $entry[0],
      ${ $entry[1]{"hit"} }[0] );
    $hit_column_counts{$hit_aln_pos}++;
    $hit_match_len{ ${ $entry[1]{"hit"} }[1] }++;

    #do the same for the query TIR
    $query_aln_pos =
      $trimmed_aln_obj->column_from_residue_number( $entry[0],
      ${ entry [1]{"query"} }[0] );
    $query_column_counts{$query_aln_pos}++;
    $query_match_len{ ${ $entry[1]{"query"} }[1] }++;

#Storing column info for $entry[0] in first ggsearch round\nLeft TIR start: \${entry[1]{'hit'}}[0]\tRight TIR start: \${entry[1]{'query'}}[0]\n";
  }

#sort the hit and query column and match length hashes by largest count to smallest
  @sorted_hitcolumn_keys =
    sort { $hit_column_counts{$b} <=> $hit_column_counts{$a} }
    keys(%hit_column_counts);
  @sorted_querycolumn_keys =
    sort { $query_column_counts{$b} <=> $query_column_counts{$a} }
    keys(%query_column_counts);
  @sorted_hit_len_keys =
    sort { $hit_match_len{$b} <=> $hit_match_len{$a} } keys(%hit_match_len);
  @sorted_query_len_keys =
    sort { $query_match_len{$b} <=> $query_match_len{$a} }
    keys(%query_match_len);

  print
"\@sorted_hitcolumn_keys: $sorted_hitcolumn_keys[0]  \@sorted_hit_len_keys: $sorted_hit_len_keys[0]\n\@sorted_querycolumn_keys: $sorted_querycolumn_keys[0]  \@sorted_hit_len_keys: $sorted_hit_len_keys[0]\n";
}

#repeat above but shifted out 3bp on each end to catch cases where found TIR is slightly off
my $adjustment = 50;
my %trim_left_pos_hash2;
my %trim_right_pos_hash2;
my @good_aln2;

foreach my $seq_obj ( $trimmed_aln_obj->each_seq() ) {
  my $left_res_pos_obj  = $seq_obj->location_from_column($left_tir_start1);
  my $right_res_pos_obj = $seq_obj->location_from_column($right_tir_start1);
  my $seq               = $seq_obj->seq();
  my $seq_len           = length($seq);
  my $left_res_pos      = ( $left_res_pos_obj->start() ) - $adjustment;
  $left_res_pos = $left_res_pos >= 0 ? $left_res_pos : 0;
  my $left_res = substr( $seq, $left_res_pos, 1 );
  my $right_res_pos = ( $right_res_pos_obj->start() ) + $adjustment;
  $right_res_pos = $right_res_pos < $seq_len ? $right_res_pos : $seq_len;
  my $right_res       = substr( $seq, $right_res_pos, 1 );
  my $seq_name        = $seq_obj->id();
  my $right_flank_pos = $seq_len - $right_res_pos;
  $trim_left_pos_hash2{$seq_name}  = $left_res_pos;
  $trim_right_pos_hash2{$seq_name} = $right_res_pos;
}

my @bad_remove2;
my @bad_aln2;

#grab all sequences from trimmed alignment
foreach my $seq_obj ( $trimmed_aln_obj->each_seq() ) {

  #grab sequence name and shorten it for use in some output filenames
  my $seq_name    = $seq_obj->id();
  my @seqn_part   = split( ":", $seq_name );
  my $fname_start = $seqn_part[0];

  #grab sequence, make a copy, strip leading and trailing hyphens
  my $seq = $seq_obj->seq();
  $seq =~ s/-//g;

  #initialize variables
  my $sub_seq;
  my $first;
  my $last;
  my $first_path = File::Spec->catpath( $volume, $out_path, "first_half.txt" );
  my $last_path  = File::Spec->catpath( $volume, $out_path, "last_half.txt" );

  #get end  sequences
  $first =
    substr( $seq, $trim_left_pos_hash2{$seq_name} - 50, $ele_half_len + 50 );
  $last = substr(
    $seq,
    $trim_right_pos_hash2{$seq_name} - $ele_half_len,
    $ele_half_len + 50
  );

  #save the two ends as files to use as inputs for a ggsearch search
  open( my $first_out, ">", $first_path )
    or die "Error creating $first_path. $!\n";
  open( my $last_out, ">", $last_path )
    or die "Error creating $last_path. $!\n";
  print $first_out ">first\n$first\n";
  print $last_out ">last\n$last\n";
  close($first_out);
  close($last_out);

  #create fasta output filename then call ggsearch
  my $out_opt =
    File::Spec->catpath( $volume, $out_path, $fname_start . ".ggsearch.out2" );
  system("ggsearch36 -n -i -T 1 -d 1 $last_path $first_path > $out_opt");

  my @tir_match_result = match_tirs( $seq_obj, $out_opt, 2 );

  if ( $tir_match_result[0] == 1 ) {
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

if ( $good_aln2_len == 0 ) {
  print "No TIRs found from MSA: 2nd round\n";
}
else {
  print "Found TIRs in MSA: 2nd round\n";

  #go through each entry of @good_aln2 array
  foreach my $row_ref (@good_aln2) {
    @entry2 = @{$row_ref};

#get column positions in of TIRs, using sequence name and index position of TIR in full sequence.
    $hit_aln_pos2 =
      $trimmed_aln_obj->column_from_residue_number( $entry2[0],
      ${ $entry2[1]{"hit"} }[0] );
    $hit_column_counts2{$hit_aln_pos2}++;
    $hit_match_len2{ ${ $entry2[1]{"hit"} }[1] }++;

    #do the same for the query (left) TIR
    $query_aln_pos2 =
      $trimmed_aln_obj->column_from_residue_number( $entry2[0],
      ${ $entry2[1]{"query"} }[0] );
    $query_column_counts2{$query_aln_pos2}++;
    $query_match_len2{ ${ $entry2[1]{"query"} }[1] }++;
  }

#sort the hit and query column and match length hashes by largest count to smallest
  @sorted_hitcolumn_keys2 =
    sort { $hit_column_counts2{$b} <=> $hit_column_counts2{$a} }
    keys(%hit_column_counts2);
  @sorted_querycolumn_keys2 =
    sort { $query_column_counts2{$b} <=> $query_column_counts2{$a} }
    keys(%query_column_counts2);
  @sorted_hit_len_keys2 =
    sort { $hit_match_len2{$b} <=> $hit_match_len2{$a} } keys(%hit_match_len2);
  @sorted_query_len_keys2 =
    sort { $query_match_len2{$b} <=> $query_match_len2{$a} }
    keys(%query_match_len2);
}

my $good_aln_len2 = @good_aln2;
my $bad_aln_len   = @bad_aln;
my $bad_aln_len2  = @bad_aln2;

print "Good aln len:  $good_aln_len  Good aln len2: $good_aln_len2\n";
print "Bad aln len:  $bad_aln_len  Bad aln len2: $bad_aln_len2\n";

#check which run generated the longer TIRs
#print "test len: $good_aln_len2 ?? $good_aln_len\n";
if ( $good_aln_len2 == 0 and $good_aln_len == 0 ) {
  my $bad_out_path =
    File::Spec->catpath( $volume, $out_path, $filename . ".bad" );
  error_out( $bad_out_path,
    "$filename\tThe two ggsearch runs failed to find TIRs" );
}
elsif ( $good_aln_len2 == 0 and $good_aln_len != 0 ) {
  print "Alignment set 1 better\n\n";
}
elsif ( $good_aln_len2 != 0 and $good_aln_len == 0 ) {
  print "Alignment set 2 better\n\n";
  @sorted_hitcolumn_keys   = @sorted_hitcolumn_keys2;
  @sorted_querycolumn_keys = @sorted_querycolumn_keys2;
  @sorted_hit_len_keys     = @sorted_hit_len_keys2;
  @sorted_query_len_keys   = @sorted_query_len_keys2;
  @bad_remove              = @bad_remove2;
}
elsif ( $good_aln_len2 > $good_aln_len ) {
  print "Alignment set 2 better\n";
  @sorted_hitcolumn_keys   = @sorted_hitcolumn_keys2;
  @sorted_querycolumn_keys = @sorted_querycolumn_keys2;
  @sorted_hit_len_keys     = @sorted_hit_len_keys2;
  @sorted_query_len_keys   = @sorted_query_len_keys2;
  @bad_remove              = @bad_remove2;
}
elsif ( $good_aln_len2 == $good_aln_len ) {
  if (
    (
          $sorted_hit_len_keys2[0] > $sorted_hit_len_keys[0]
      and $sorted_query_len_keys2[0] > $sorted_query_len_keys[0]
    )
    and ( $sorted_hitcolumn_keys2[0] >= $sorted_hitcolumn_keys[0]
      and $sorted_querycolumn_keys2[0] >= $sorted_querycolumn_keys[0] )
    )
  {
    print "Alignment set 2 better\n";
    @sorted_hitcolumn_keys   = @sorted_hitcolumn_keys2;
    @sorted_querycolumn_keys = @sorted_querycolumn_keys2;
    @sorted_hit_len_keys     = @sorted_hit_len_keys2;
    @sorted_query_len_keys   = @sorted_query_len_keys2;
    @bad_remove              = @bad_remove2;
  }
  else {
    print "Alignment set 1 better\n";
  }
}
else {
  print "Alignment set 1 better\n\n";
}

#open a file to store info on why analysis of a copy was aborted
my $removed_out_path =
  File::Spec->catpath( $volume, $out_path, $filename . ".removed_sequences" );
open( my $removed_out, ">", $removed_out_path );
foreach my $seq_obj (@bad_remove) {
  my $seq_name = $seq_obj->id();
  print $removed_out "$seq_name\tNo TIRs found by first two ggsearch runs\n";
  $trimmed_aln_obj->remove_seq($seq_obj);
}

print
"before get_tir_nt_positions: ( $sorted_hitcolumn_keys[0],$sorted_hit_len_keys[0],$sorted_querycolumn_keys[0],$sorted_query_len_keys[0])\n";
my ($ref2_tir_positions) = get_tir_nt_positions(
  $trimmed_aln_obj,            \%tir_positions,
  $sorted_hitcolumn_keys[0],   $sorted_hit_len_keys[0],
  $sorted_querycolumn_keys[0], $sorted_query_len_keys[0]
);
%tir_positions = %$ref2_tir_positions;

#print "check errfile for tir_positions just before rereading org aln:\n";
#warn "tir_positions just before rereading org aln:\n";
#warn Dumper \%tir_positions;

#reread the original alignment file
my $in_obj2     = Bio::AlignIO->new( -file => $infile, -format => 'fasta' );
my $ori_aln_obj = $in_obj2->next_aln();
my $ori_aln_len = $ori_aln_obj->length();

## because we adjusted the left tir the tir_pos is off by the number of gaps removed. but gaps were only removed from a part of the whole.
## need to go back to org, best tir, before adjust, substract the number of nt offset by adjust and by gaps, then use that col
if ($left_tir_adjusted or $right_tir_adjusted){
  my %col;
  my $tir;
  foreach my $seq_obj($trimmed_aln_obj->each_seq){
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
    print substr($seq_id,0,4),":$left_gaps...$before_left_tir...$element...$after_right_tir...$right_gaps\n";
    my $new_seq = $left_gaps.$before_left_tir.$element.$after_right_tir.$right_gaps;
    $seq_obj->seq($new_seq);
  }
  
}
#get the column positions of tirs in the original alignment
print "before get_col: lts:$left_tir_start, rts:$right_tir_start,\n";
my ( $left_tir_end, $right_tir_end );
( $left_tir_start, $right_tir_start, $left_tir_end, $right_tir_end ) = get_columns( $ori_aln_obj, \%tir_positions, 2 );
print
"after get_col: $left_tir_start, $right_tir_start,$left_tir_end, $right_tir_end \n";

## cleaning up aln using tir_starts found with get_columns
## and overwriting the tir_starts with the adj values
( $right_tir_start, $left_tir_start ) =
  adjust_tir_starts( $ori_aln_obj, $right_tir_start, $left_tir_start );
print "rt:$right_tir_start , lt:$left_tir_start\n";

my $tir_length = $left_tir_end - $left_tir_start;

#print
#"First column grab just after reimporting the original alignment\nLeft Tir Start: $left_tir_start  Left Tir End: $left_tir_end\nRight Tir End: $right_tir_end  Right TIR Start: $right_tir_start\n";
#print "Starting MSA length = $ori_aln_len TIR length = $tir_length (filename.all_adj)\n";

#go back through the full alignment and remove only the sequences that create gaps in the TIRs or that failed the first two TIR searches
my %gap_seq_remove2;
my %search_tirs;

for ( my $i = 0 ; $i < $ori_aln_len ; $i++ ) {
  if ( ( $i >= $left_tir_start - 1 and $i <= $left_tir_end - 1 )
    or ( $i >= $right_tir_end - 1 and $i <= $right_tir_start - 1 ) )
  {
    my $id_row_ref      = $full_id_array[$i];
    my @id_info         = @{$id_row_ref};
    my @gap_col_array   = @{$gap_cols};
    my $gap_col_hashref = $gap_col_array[$i];
    my %gap_col_hash    = %{$gap_col_hashref};
    my $base_count      = 0;
    my $total_count     = 0;
    foreach my $key ( keys %gap_col_hash ) {

      if ( $gap_col_hash{$key} != 1 ) {
        $base_count++;
      }
      $total_count++;
    }
    my $pos_present = $base_count / $total_count;
    if ( $gap_id_array[$i] >= 50.0 and $id_info[1] <= 50 or $pos_present < .5 )
    {
      foreach my $key ( keys %gap_col_hash ) {
        if ( $gap_col_hash{$key} != 1 ) {
          $gap_seq_remove2{$key}++;
        }
      }
    }
    elsif ( $gap_id_array[$i] < 50.0 ) {
      foreach my $key ( keys %gap_col_hash ) {
        if ( $gap_col_hash{$key} == 1 ) {
          $gap_seq_remove2{$key}++;
        }
      }
    }
  }
}

foreach my $row_ref (@gap_seq_pos_remove) {
  @entry = @{$row_ref};
  if ( exists $gap_seq_remove2{ $entry[0] } ) {
    next;
  }
  else {
    $search_tirs{ $entry[0] }++;
  }
}

foreach my $key ( keys %gap_seq_remove ) {
  if ( exists $gap_seq_remove2{$key} ) {
    next;
  }
  else {
    $search_tirs{$key}++;
  }
}

#keep track of removed sequences to print to file
my %removed_seq_hash;
foreach my $key ( keys %gap_seq_remove2 ) {
  my $remove = $ori_aln_obj->get_seq_by_id($key);
  my $seq_id = $remove->id();
  print $removed_out
"$seq_id\tStep2:Sequence caused or contained gaps in at least one of the TIRs\n";
  $removed_seq_hash{$seq_id}++;
  $ori_aln_obj->remove_seq($remove);
}

my $final_aln_obj = $ori_aln_obj->remove_gaps( '-', 1 );
my $int_aln_out =
  File::Spec->catpath( $volume, $out_path, $filename . ".intermediate" );
$out = Bio::AlignIO->new(
  -file             => ">$int_aln_out",
  -format           => 'fasta',
  -displayname_flat => 0
);
$out->write_aln($final_aln_obj);
my $final_len = $final_aln_obj->length();
print "MSA length is $final_len now\n";

#get the column positions of tirs in the intermediate alignment
( $left_tir_start, $right_tir_start, $left_tir_end, $right_tir_end ) =
  get_columns( $final_aln_obj, \%tir_positions, 2 );
print
"2nd column grab after removing some TIR disrupting copies and removing gap only columns\nLeft Tir Start: $left_tir_start  Left Tir End: $left_tir_end\nRight Tir End: $right_tir_end  Right TIR Start: $right_tir_start (filename.intermediate)\n";

my $new_gap_cols = $final_aln_obj->gap_col_matrix();
my %new_gap_seq_remove;

for ( my $i = 0 ; $i < $final_len ; $i++ ) {
  if ( ( $i >= $left_tir_start - 1 and $i <= $left_tir_end - 1 )
    or ( $i >= $right_tir_end - 1 and $i <= $right_tir_start - 1 ) )
  {
    my $trim_gap_col_hashref = $new_gap_cols->[$i];
    my $base_count           = 0;
    my $total_count          = 0;
    foreach my $key ( keys %{$trim_gap_col_hashref} ) {
      if ( $trim_gap_col_hashref->{$key} != 1 ) {
        $base_count++;
      }
      $total_count++;
    }
    my $pos_present = $base_count / $total_count;
    if ( $pos_present < 0.5 ) {
      foreach my $key ( keys %{$trim_gap_col_hashref} ) {
        if ( $trim_gap_col_hashref->{$key} != 1 ) {
          my $seq_obj = $final_aln_obj->get_seq_by_id($key);
          $new_gap_seq_remove{$key} = $seq_obj;
        }
      }
    }
    else {
      foreach my $key ( keys %{$trim_gap_col_hashref} ) {
        if ( $trim_gap_col_hashref->{$key} == 1 ) {
          my $seq_obj = $final_aln_obj->get_seq_by_id($key);
          $new_gap_seq_remove{$key} = $seq_obj;
        }
      }
    }
  }
}

foreach my $key ( keys %new_gap_seq_remove ) {
  my $seq_obj = $new_gap_seq_remove{$key};
  my $seq_id  = $seq_obj->id();
  print $removed_out
    "$seq_id\tSequence caused or contained gaps in at least one of the TIRs\n";
  $final_aln_obj->remove_seq($seq_obj);
}
$final_aln_obj = $final_aln_obj->remove_gaps( '-', 1 );
my $check_len = $final_aln_obj->length();
print "MSA length now: $check_len\n";

#get the column positions of tirs in the intermediate alignment
( $left_tir_start, $right_tir_start, $left_tir_end, $right_tir_end ) =
  get_columns( $final_aln_obj, \%tir_positions, 2 );
print
"3rd column grab after removing more copies with TIR isssues\nLeft Tir Start: $left_tir_start  Left Tir End: $left_tir_end\nRight Tir End: $right_tir_end  Right TIR Start: $right_tir_start ($filename.intermediate)\n";
my %to_remove_from_final;
foreach my $seq_name ( keys %search_tirs ) {
  my $seq_obj = $final_aln_obj->get_seq_by_id($seq_name);
  next if !defined $seq_obj;
  if ( exists $gap_seq_remove2{$seq_name} ) {
    $to_remove_from_final{$seq_name} = $seq_obj;
    next;
  }
  if ( exists $removed_seq_hash{$seq_name} ) {
    $to_remove_from_final{$seq_name} = $seq_obj;
    next;
  }
  if ( exists $new_gap_seq_remove{$seq_name} ) {
    $to_remove_from_final{$seq_name} = $seq_obj;
    next;
  }
  my @seqn_part   = split( ":", $seq_name );
  my $fname_start = $seqn_part[0];
  my $seq         = $seq_obj->seq();
  my $seq_len     = length($seq);
  print "Error $seq_name:seq_len is 0\n" if $seq_len == 0;

  my $left_test_pos  = substr( $seq, $left_tir_start - 1,  1 );
  my $right_test_pos = substr( $seq, $right_tir_start - 1, 1 );
  if ( $left_test_pos eq "-" or $right_test_pos eq "-" ) {
    print $removed_out
      "$seq_name\tContains a gap at the start of one or both TIRs\n";
    $final_aln_obj->remove_seq($seq_obj);
    delete $to_remove_from_final{$seq_name}
      if exists $to_remove_from_final{$seq_name};
    next;
  }

  my $left_tir_start_obj  = $seq_obj->location_from_column($left_tir_start);
  my $left_tir_start_pos  = $left_tir_start_obj->start();
  my $right_tir_start_obj = $seq_obj->location_from_column($right_tir_start);
  my $right_tir_start_pos = $right_tir_start_obj->start();

  $seq =~ s/-//g;
  my $first;
  my $last;
  my $first_path = File::Spec->catpath( $volume, $out_path, "first_half.txt" );
  my $last_path  = File::Spec->catpath( $volume, $out_path, "last_half.txt" );

  $first = substr( $seq, $left_tir_start_pos - 1, $ele_half_len );
  $last = substr( $seq, $right_tir_start_pos - $ele_half_len, $ele_half_len );

  #continue TIR search
  open( my $first_out, ">", $first_path )
    or die "Error creating $first_path. $!\n";
  open( my $last_out, ">", $last_path )
    or die "Error creating $last_path. $!\n";
  print $first_out ">first\n$first\n";
  print $last_out ">last\n$last\n";
  close($first_out);
  close($last_out);
  my $out_opt =
    File::Spec->catpath( $volume, $out_path, $fname_start . ".ggsearch3.out" );
  system("ggsearch36 -n -i -T 8 -d 1 $last_path $first_path > $out_opt");
  my @tir_match_result = match_tirs( $seq_obj, $out_opt, 3 );

  if ( $tir_match_result[0] == 0 or !@tir_match_result ) {
    my $seq_id = $seq_obj->id();
    $final_aln_obj->remove_seq($seq_obj);
    delete $to_remove_from_final{$seq_name}
      if exists $to_remove_from_final{$seq_name};
    print $removed_out "$seq_id\tNo TIR matches found by last ggsearch run\n";
  }
}
foreach my $seq_id ( keys %to_remove_from_final ) {
  my $seq_obj = $final_aln_obj->get_seq_by_id($seq_id);
  if ( defined $seq_obj ) {
    $final_aln_obj->remove_seq($seq_obj) if defined $seq_obj;
  }
  else {
    print "$seq_id in final_aln_obj is not defined\n";
  }
}

$final_aln_obj = $final_aln_obj->remove_gaps( '-', 1 );
my $int_aln_out2 =
  File::Spec->catpath( $volume, $out_path, $filename . ".intermediate2" );
$out = Bio::AlignIO->new(
  -file             => ">$int_aln_out2",
  -format           => 'fasta',
  -displayname_flat => 0
);
$out->write_aln($final_aln_obj);
my $last_len = $final_aln_obj->length();

#get the column positions of tirs in the final alignment
( $left_tir_start, $right_tir_start, $left_tir_end, $right_tir_end ) =
  get_columns( $final_aln_obj, \%tir_positions, 2 );

#print
#"Column grab after third ggsearch run\nLeft Tir Start: $left_tir_start  Right: $right_tir_start ($filename.intermediate2)\n";
#my @gap_seq_pos_remove3;
#($left_tir_start, $right_tir_start , $ref2array, $trimmed_aln_obj, $ref2hash )
#   = consensus_filter(\@gap_seq_pos_remove3, $final_aln_obj,$left_tir_start,$right_tir_start, \%tir_positions, 'final');
#@gap_seq_pos_remove3 = @$ref2array;
#%tir_positions = %$ref2hash;
#print "after last consenus filter on final_aln. LTStart: $left_tir_start  Right: $right_tir_start ($filename.final)\n";
#
#foreach my $info ( @gap_seq_pos_remove3 ) {
#  my $seq_id = ${$info}[0];
#  print $removed_out
#    "$seq_id\tStep3:Sequence contained too many mismataches in at least one of the TIRs\n";
#}

#Extract the left and right TIRs as new alignment objects
my $left_TIR_aln_obj =
  $final_aln_obj->slice( $left_tir_start, $left_tir_end, 1 );
my $right_TIR_aln_obj =
  $final_aln_obj->slice( $right_tir_end, $right_tir_start, 1 );
my $element_aln_obj =
  $final_aln_obj->slice( $left_tir_start, $right_tir_start, 1 );

#find TSDs based off column # of TIRs. Must check for 2-bp and 4-bp TSD included in TIRS and then a longer TSD. Longer TSDs that start and finish with same 2- or 4-bp as the 2- or 4-bp TSDs supersede the shorter TSDs. Also, grab flanks and save them to file.
my @putative_TSD;
my @TSD_info;
my @no_TSD_found;
my @put_TSD_names;
my $flanks_out_path =
  File::Spec->catpath( $volume, $out_path, $filename . "_flanks.fa" );
open( my $flanks_out, ">", $flanks_out_path )
  or die "Error creating $flanks_out_path. $!\n";
my $tsd1_count = 0;
my $tsd2_count = 0;
my $tsd3_count = 0;

print "Starting code to find TSDs\n";


my $no_TSD_found_out_path =
  File::Spec->catpath( $volume, $out_path, $filename . ".TSD_issues.info" );
open( my $no_TSD_found_out, ">", $no_TSD_found_out_path );

foreach my $seq_obj ( $final_aln_obj->each_seq() ) {
  my $seq_name    = $seq_obj->id();
  my @seqn_part   = split( ":", $seq_name );
  my $fname_start = $seqn_part[0];
  my $seq_ori     = $seq_obj->seq();
  my $seq         = $seq_ori;
  $seq =~ s/-//g;
  my $left_tsd;
  my $left_tsd_substr;
  my $starting_left_flank;
  my $right_tsd;
  my $right_tsd_substr;
  my $starting_right_flank;
  my $last = 0;

#get seq_pos of the start of TIRs from column # in alignment. Adjust to include the 4-bp at outside ends of TIRs and then split into the 2-, 4-, and 10-bp sequences to search for TSDs
  my $left_tsd_loc_obj = $seq_obj->location_from_column($left_tir_start);
  my $left_tsd_end_pos = $left_tsd_loc_obj->start();
  $starting_left_flank = substr( $seq, 0, $left_tsd_end_pos + 4 );
  $left_tsd = substr( $seq, $left_tsd_end_pos - 21, 24 );
  if (length $left_tsd < 20){
    my $message = "Left flank too short to look at TSDs.\n";
    print $no_TSD_found_out $message;
    next;
  }
  $left_tsd_substr = substr( $left_tsd, 10, 10 );

  my $right_tsd_loc_obj   = $seq_obj->location_from_column($right_tir_start);
  my $right_tsd_start_pos = $right_tsd_loc_obj->start();
  $starting_right_flank = substr( $seq,       $right_tsd_start_pos - 4 );
  $right_tsd            = substr( $seq,       $right_tsd_start_pos - 4, 24 );
  if (length $right_tsd < 14){
    my $message = "Right flank too short to look at TSDs.";
    print $no_TSD_found_out $message;
    next;
  }
  $right_tsd_substr     = substr( $right_tsd, 4, 10 );

  my $left_2bp  = substr( $left_tsd,  -4, 2 );
  my $right_2bp = substr( $right_tsd, 2,  2 );
  my $left_4bp  = substr( $left_tsd,  -4, 4 );
  my $right_4bp = substr( $right_tsd, 0,  4 );
  my $tsd1;
  my $tsd2;
  my $tsd3;

  if ( $left_2bp eq $right_2bp ) {
    print "Round 1 TSD search - Yes! $fname_start\n";
    $tsd1 = $left_2bp;
  }
  if ( $left_4bp eq $right_4bp ) {
    print "Round 2 TSD search - Yes! $fname_start\n";
    $tsd2 = $left_4bp;
  }
  for ( my $i = 0 ; $i < 9 ; $i++ ) {
    if (
      substr( $left_tsd_substr, $i ) eq substr( $right_tsd_substr, 0, -($i) ) )
    {
      print "Round 3 TSD search - Yes! $i $fname_start\n";
      $tsd3 = substr( $left_tsd_substr, $i );
      next;
    }
  }

  #Save found TSD to an array or report that a TSD wasn't found
  if ($tsd1) {
    if ( !$tsd2 and !$tsd3 ) {
      my $insertion_site =
        substr( $left_tsd, -14, 12 ) . substr( $right_tsd, 4, 10 );
      push @TSD_info, [ $seq_name, $insertion_site, $tsd1 ];
      push @putative_TSD, $tsd1;
      push @put_TSD_names, [ $seq_name, $tsd1 ];
      my $left_flank  = substr( $starting_left_flank,  0, -2 );
      my $right_flank = substr( $starting_right_flank, 2 );
      my $flanks      = $left_flank . $right_flank;
      print $flanks_out ">$fname_start\n$flanks\n";
      print "tsd1: $fname_start\n";
      $tsd1_count++;
    }
    elsif ( !$tsd2 and $tsd3 ) {
      if (  ( length($tsd3) > length($tsd1) )
        and ( substr( $tsd3, 0, 2 ) eq $tsd1 )
        and ( substr( $tsd3, -2 ) eq $tsd1 ) )
      {
        my $insertion_site = substr(
          $left_tsd,
          ( -4 - length($tsd3) - 10 ),
          ( length($tsd3) + 10 )
        ) . substr( $right_tsd, 4 + length($tsd3), 10 );
        push @TSD_info, [ $seq_name, $insertion_site, $tsd3 ];
        push @putative_TSD, $tsd3;
        push @put_TSD_names, [ $seq_name, $tsd3 ];
        my $left_flank =
          substr( $starting_left_flank, 0, ( -2 - length($tsd3) ) );
        my $right_flank =
          substr( $starting_right_flank, ( 2 + length($tsd3) ) );
        my $flanks = $left_flank . $right_flank;
        print "tsd3: $fname_start\n";
        print $flanks_out ">$fname_start\n$flanks\n";
        $tsd3_count++;
      }
      else {
        my $insertion_site =
          substr( $left_tsd, -12 ) . substr( $right_tsd, 2, 10 );
        push @TSD_info, [ $seq_name, $insertion_site, $tsd1 ];
        push @putative_TSD, $tsd1;
        push @put_TSD_names, [ $seq_name, $tsd1 ];
        my $left_flank  = substr( $starting_left_flank,  0, -2 );
        my $right_flank = substr( $starting_right_flank, 2 );
        my $flanks      = $left_flank . $right_flank;
        print "tsd1: $fname_start\n";
        print $flanks_out ">$fname_start\n$flanks\n";
        $tsd1_count++;
      }
    }
    elsif ( $tsd2 and !$tsd3 ) {
      if ( substr( $tsd2, 0, 2 ) eq $tsd1 and substr( $tsd2, -2 ) eq $tsd1 ) {
        my $insertion_site =
          substr( $left_tsd, (-14) ) . substr( $right_tsd, 4, 10 );
        push @TSD_info, [ $seq_name, $insertion_site, $tsd2 ];
        push @putative_TSD, $tsd2;
        push @put_TSD_names, [ $seq_name, $tsd2 ];
        my $left_flank = substr( $starting_left_flank, 0, -4 - length($tsd2) );
        my $right_flank = substr( $starting_right_flank, 4 + length($tsd2) );
        my $flanks = $left_flank . $right_flank;
        print "tsd2: $fname_start\n";
        print $flanks_out ">$fname_start\n$flanks\n";
        $tsd2_count++;
      }
      else {
        my $insertion_site =
          substr( $left_tsd, -14, 12 ) . substr( $right_tsd, 2, 10 );
        push @TSD_info, [ $seq_name, $insertion_site, $tsd1 ];
        push @putative_TSD, $tsd1;
        push @put_TSD_names, [ $seq_name, $tsd1 ];
        my $left_flank  = substr( $starting_left_flank,  0, -4 );
        my $right_flank = substr( $starting_right_flank, 4 );
        my $flanks      = $left_flank . $right_flank;
        print "tsd1: $fname_start\n";
        print $flanks_out ">$fname_start\n$flanks\n";
        $tsd1_count++;
      }
    }
    elsif ( $tsd2 and $tsd3 ) {
      if (  ( substr( $tsd3, 0, 2 ) eq $tsd1 )
        and ( substr( $tsd3, -2 ) eq $tsd1 ) )
      {
        my $insertion_site = substr(
          $left_tsd,
          ( -4 - length($tsd3) - 10 ),
          ( length($tsd3) + 10 )
        ) . substr( $right_tsd, 4 + length($tsd3), 10 );
        push @TSD_info, [ $seq_name, $insertion_site, $tsd3 ];
        push @putative_TSD, $tsd3;
        push @put_TSD_names, [ $seq_name, $tsd3 ];
        my $left_flank =
          substr( $starting_left_flank, 0, ( -4 - length($tsd3) ) );
        my $right_flank =
          substr( $starting_right_flank, ( 4 + length($tsd3) ) );
        my $flanks = $left_flank . $right_flank;
        print "tsd3: $fname_start\n";
        print $flanks_out ">$fname_start\n$flanks\n";
        $tsd3_count++;
      }
      elsif ( ( substr( $tsd2, 0, 2 ) eq $tsd1 )
        and ( substr( $tsd2, -2 ) eq $tsd1 ) )
      {
        my $insertion_site = substr(
          $left_tsd,
          ( -2 - length($tsd2) - 10 ),
          ( length($tsd2) + 10 )
        ) . substr( $right_tsd, 2 + length($tsd2), 10 );
        push @TSD_info, [ $seq_name, $insertion_site, $tsd2 ];
        push @putative_TSD, $tsd2;
        push @put_TSD_names, [ $seq_name, $tsd2 ];
        my $left_flank = substr( $starting_left_flank, 0, -2 - length($tsd2) );
        my $right_flank = substr( $starting_right_flank, 2 + length($tsd2) );
        my $flanks = $left_flank . $right_flank;
        print "tsd2: $fname_start\n";
        print $flanks_out ">$fname_start\n$flanks\n";
        $tsd2_count++;
      }
      else {
        my $insertion_site =
          substr( $left_tsd, -14, 12 ) . substr( $right_tsd, 4, 10 );
        push @TSD_info, [ $seq_name, $insertion_site, $tsd1 ];
        push @putative_TSD, $tsd1;
        push @put_TSD_names, [ $seq_name, $tsd1 ];
        my $left_flank  = substr( $starting_left_flank,  0, -4 );
        my $right_flank = substr( $starting_right_flank, 4 );
        my $flanks      = $left_flank . $right_flank;
        print "tsd1: $fname_start\n";
        print $flanks_out ">$fname_start\n$flanks\n";
        $tsd1_count++;
      }
    }
  }
  elsif ($tsd3) {
    if ( !$tsd2 ) {
      my $insertion_site =
        substr( $left_tsd, ( -4 - length($tsd3) - 10 ), ( length($tsd3) + 10 ) )
        . substr( $right_tsd, 4 + length($tsd3), 10 );
      push @TSD_info, [ $seq_name, $insertion_site, $tsd3 ];
      push @putative_TSD, $tsd3;
      push @put_TSD_names, [ $seq_name, $tsd3 ];
      my $left_flank = substr( $starting_left_flank, 0, -4 - length($tsd3) );
      my $right_flank = substr( $starting_right_flank, ( 4 + length($tsd3) ) );
      my $flanks = $left_flank . $right_flank;
      print "tsd3: $fname_start\n";
      print $flanks_out ">$fname_start\n$flanks\n";
      $tsd3_count++;
    }
    else {
      if (  ( length($tsd3) > length($tsd2) )
        and ( substr( $tsd3, 0, 4 ) eq $tsd2 )
        and ( substr( $tsd3, -4 ) eq $tsd2 ) )
      {
        my $insertion_site = substr(
          $left_tsd,
          ( -4 - length($tsd3) - 10 ),
          ( length($tsd3) + 10 )
        ) . substr( $right_tsd, 4 + length($tsd3), 10 );
        push @TSD_info, [ $seq_name, $insertion_site, $tsd3 ];
        push @putative_TSD, $tsd3;
        push @put_TSD_names, [ $seq_name, $tsd3 ];
        my $left_flank = substr( $starting_left_flank, 0, -4 - length($tsd3) );
        my $right_flank =
          substr( $starting_right_flank, ( 4 + length($tsd3) ) );
        my $flanks = $left_flank . $right_flank;
        print "tsd3: $fname_start\n";
        print $flanks_out ">$fname_start\n$flanks\n";
        $tsd3_count++;
      }
      else {
        my $insertion_site =
          substr( $left_tsd, -14 ) . substr( $right_tsd, 4, 10 );
        push @TSD_info, [ $seq_name, $insertion_site, $tsd2 ];
        push @putative_TSD, $tsd2;
        push @put_TSD_names, [ $seq_name, $tsd2 ];
        my $left_flank  = substr( $starting_left_flank,  0, -4 );
        my $right_flank = substr( $starting_right_flank, 4 );
        my $flanks      = $left_flank . $right_flank;
        print "tsd2: $fname_start\n";
        print $flanks_out ">$fname_start\n$flanks\n";
        $tsd2_count++;
      }
    }
  }
  elsif ($tsd2) {
    my $insertion_site = substr( $left_tsd, -14 ) . substr( $right_tsd, 4, 10 );
    push @TSD_info, [ $seq_name, $insertion_site, $tsd2 ];
    push @putative_TSD, $tsd2;
    push @put_TSD_names, [ $seq_name, $tsd2 ];
    my $left_flank  = substr( $starting_left_flank,  0, -4 );
    my $right_flank = substr( $starting_right_flank, 4 );
    my $flanks      = $left_flank . $right_flank;
    print "tsd2: $fname_start\n";
    print $flanks_out ">$fname_start\n$flanks\n";
    $tsd2_count++;
  }
  else {

#store or print to report that TSD was not found for this copy use $seq_name, etc
#use the end position of tirs to grab the flanking sequences
    my $left_flank  = substr( $starting_left_flank,  0, -4 );
    my $right_flank = substr( $starting_right_flank, 4 );
    my $flanks      = $left_flank . $right_flank;
    push @no_TSD_found, $seq_name;
    print "no TSD: $fname_start\n";
    print $flanks_out ">$fname_start\n$flanks\n";
  }
}

close($flanks_out);
print "Finished code to find TSDs\n\n";

my $final_align_len = $final_aln_obj->num_sequences();
$element_info{"copy_num"} = $final_align_len;
my $final_aln_out =
  File::Spec->catpath( $volume, $out_path, $filename . ".final" );
$out = Bio::AlignIO->new(
  -file             => ">$final_aln_out",
  -format           => 'fasta',
  -displayname_flat => 0
);
$out->write_aln($final_aln_obj);

my $no_TSD_found_len = @no_TSD_found;
if ( $no_TSD_found_len >= 1 ) {
  print "Copies without TSDs removed and printed to file\n";
  foreach my $item (@no_TSD_found) {
    print $no_TSD_found_out "$item\tNo TSDs found\n";
  }
}

#if necessary, change TIRS to remove the 2bp TSD that are included in them & calculate the overall percent identity of each TIR
if ( $tsd1_count > $tsd2_count and $tsd1_count > $tsd3_count ) {
  print "Adjusting TIRs by 2bp\n\n";
  $left_tir_start += 2;
  $right_tir_start -= 2;
  $left_TIR_aln_obj =
    $final_aln_obj->slice( $left_tir_start, $left_tir_end, 1 );
  $right_TIR_aln_obj =
    $final_aln_obj->slice( $right_tir_end, $right_tir_start, 1 );
  $element_aln_obj =
    $final_aln_obj->slice( $left_tir_start, $right_tir_start, 1 );
  foreach my $seq_obj ( $final_aln_obj->each_seq() ) {
    my $seq_name = $seq_obj->id();
    if ( !defined $tir_positions{$seq_name}{'left_tir_start'} ) {
      my $left_pos_obj  = $seq_obj->location_from_column($left_tir_start);
      my $right_pos_obj = $seq_obj->location_from_column($right_tir_start);
      my $left_pos      = $left_pos_obj->start();
      my $right_pos     = $right_pos_obj->start();
      $tir_positions{$seq_name}{'left_tir_start'}  = $left_pos;
      $tir_positions{$seq_name}{'right_tir_start'} = $right_pos;
    }
    $element_info{$seq_name}{"left_tir_start"} =
      $tir_positions{$seq_name}{'left_tir_start'} + 2;
    $element_info{$seq_name}{"right_tir_start"} =
      $tir_positions{$seq_name}{'right_tir_start'} - 2;
  }
}
elsif ( $tsd2_count > $tsd1_count and $tsd2_count > $tsd3_count ) {
  print "Adjusting TIRs by 4bp\n\n";
  $left_tir_start += 4;
  $right_tir_start -= 4;
  $left_TIR_aln_obj =
    $final_aln_obj->slice( $left_tir_start, $left_tir_end, 1 );
  $right_TIR_aln_obj =
    $final_aln_obj->slice( $right_tir_end, $right_tir_start, 1 );
  $element_aln_obj =
    $final_aln_obj->slice( $left_tir_start, $right_tir_start, 1 );
  foreach my $seq_obj ( $final_aln_obj->each_seq() ) {
    my $seq_name = $seq_obj->id();
    if ( !defined $tir_positions{$seq_name}{'left_tir_start'} ) {
      my $left_pos_obj  = $seq_obj->location_from_column($left_tir_start);
      my $right_pos_obj = $seq_obj->location_from_column($right_tir_start);
      my $left_pos      = $left_pos_obj->start();
      my $right_pos     = $right_pos_obj->start();
      $tir_positions{$seq_name}{'left_tir_start'}  = $left_pos;
      $tir_positions{$seq_name}{'right_tir_start'} = $right_pos;
    }
    $element_info{$seq_name}{"left_tir_start"} =
      $tir_positions{$seq_name}{'left_tir_start'} + 4;
    $element_info{$seq_name}{"right_tir_start"} =
      $tir_positions{$seq_name}{'right_tir_start'} - 4;
  }
}
else {
  foreach my $seq_obj ( $final_aln_obj->each_seq() ) {
    my $seq_name = $seq_obj->id();
    if ( !defined $tir_positions{$seq_name}{'left_tir_start'} ) {
      my $left_pos_obj  = $seq_obj->location_from_column($left_tir_start);
      my $right_pos_obj = $seq_obj->location_from_column($right_tir_start);
      my $left_pos      = $left_pos_obj->start();
      my $right_pos     = $right_pos_obj->start();
      $tir_positions{$seq_name}{'left_tir_start'}  = $left_pos;
      $tir_positions{$seq_name}{'right_tir_start'} = $right_pos;
    }
    $element_info{$seq_name}{"left_tir_start"} =
      $tir_positions{$seq_name}{'left_tir_start'};
    $element_info{$seq_name}{"right_tir_start"} =
      $tir_positions{$seq_name}{'right_tir_start'};
  }
}
my $left_tir_id   = $left_TIR_aln_obj->percentage_identity;
my $right_tir_id  = $right_TIR_aln_obj->percentage_identity;
my $element_id    = $element_aln_obj->percentage_identity;
my $left_tir_seq  = $left_TIR_aln_obj->consensus_string();
my $right_tir_seq = $right_TIR_aln_obj->consensus_string();

$element_info{"element_id"}    = $element_id;
$element_info{"left_tir_seq"}  = $left_tir_seq;
$element_info{"left_tir_id"}   = $left_tir_id;
$element_info{"right_tir_seq"} = $right_tir_seq;
$element_info{"right_tir_id"}  = $right_tir_id;

#Determine the most common TSD by sequence first if possible or by length. If >80% of the TSDs are the same, then that sequence is stored as the TSD for output. Otherwise, look at the lengths of the TSDs and store if >80% of the TSDs have the same length.
my %TSD_counts;
my $TSD_array_length = @put_TSD_names;
if ( $TSD_array_length == 0 ) {
  print "No TSDs found\n";
  my $bad_out_path =
    File::Spec->catpath( $volume, $out_path, $filename . ".bad" );
  error_out( $bad_out_path, "$filename\tNo TSDs found for any copy", 1 );
  my $gff_path = File::Spec->catpath( $volume, $out_path, $filename . ".gff" );
  generate_gff( $final_aln_obj, $gff_path, "TSD", \%element_info );
  exit 0;
}

#count the occurances of a TSD sequence and sort with the highest value first in an array
foreach my $row_ref (@put_TSD_names) {
  my @pos = @{$row_ref};
  $TSD_counts{ $pos[1] }++;
  $element_info{ $pos[0] }{"TSD"} = $pos[1];
}

#undef @put_TSD_names;
my @sorted_TSD_keys =
  sort { $TSD_counts{$b} <=> $TSD_counts{$a} } keys(%TSD_counts);

#check whether the same TSD sequence was found in >80% of the copies. If not, count the occurances of TSD length, sort by highest occurances, and check if the length of the TSD is the same in >80% of the copies
my $final_TSD_length;
my $final_TSD_seq;
my $TSD_fraction   = $TSD_counts{ $sorted_TSD_keys[0] } / $TSD_array_length;
my $need_consensus = 0;

if ( ($TSD_fraction) > 0.8 ) {
  $element_info{"TSD_fraction"} = $TSD_fraction;
  $final_TSD_length             = length( $sorted_TSD_keys[0] );
  $final_TSD_seq                = $sorted_TSD_keys[0];
  $element_info{"TSD_seq"}      = $final_TSD_seq;
  $element_info{"TSD_len"}      = $final_TSD_length;
  $element_info{"TSD_fraction"} = $TSD_fraction;
}
else {
  $need_consensus = 1;
  my %TSD_length_counts;
  foreach my $row (@putative_TSD) {
    $TSD_length_counts{ length($row) }++;
  }
  my @sorted_TSD_length_keys =
    sort { $TSD_length_counts{$b} <=> $TSD_length_counts{$a} }
    keys(%TSD_length_counts);
  $final_TSD_length = $sorted_TSD_length_keys[0];
  $TSD_fraction     = $TSD_array_length / $final_align_len;

  $element_info{"TSD_len"}      = $final_TSD_length;
  $element_info{"TSD_fraction"} = $TSD_fraction;
}

#undef @putative_TSD;

my $insertion_site_file = $filename . ".insertion-site.fa";
my $insertion_site_out_path =
  File::Spec->catpath( $volume, $out_path, $insertion_site_file );
my $tsd_out_path =
  File::Spec->catpath( $volume, $out_path, $filename . "_tsd.fa" );
open( my $tsd_info_out, ">", $tsd_out_path )
  or die "Error creating $tsd_out_path. $!\n";
open( my $insertion_site_out, ">", $insertion_site_out_path )
  or die "Error creating $insertion_site_out_path. $!\n";
my $insertion_num = @TSD_info;
foreach my $row_ref (@TSD_info) {

  #print "Printing insertion site info\n";
  my @pos = @{$row_ref};
  if ( length( $pos[2] ) == $element_info{"TSD_len"} ) {
    print $tsd_info_out ">$pos[0]\n$pos[2]\n";
    if ( $pos[1] !~ m/n/i ) {
      print $insertion_site_out ">$pos[0]\n$pos[1]\n";
    }
  }
  else {
    print $no_TSD_found_out
      "$pos[0]\tTSD length different than >=80% of other copies\n";
  }
}
close($insertion_site_out);
close($tsd_info_out);
close($no_TSD_found_out);

if ( $need_consensus == 1 ) {

  #read in tsd file as MSA to generate IUPAC consensus
  my $tsd_in_obj =
    Bio::AlignIO->new( -file => $tsd_out_path, -format => 'fasta' );
  my $tsd_aln_obj   = $tsd_in_obj->next_aln();
  my $tsd_consensus = $tsd_aln_obj->consensus_iupac();
  $element_info{"TSD_seq"} = $tsd_consensus;
}

#import TE characteristics from table and generate regular expressions to classify element by TIR and TSD if possible
print "\nProcessing TE charateristics table\n";
my $TE_char_path = $FindBin::Bin . "/DNA_TE_TIR-TSD.table";
open( my $TE_in, "<", $TE_char_path ) or die "Error reading $TE_char_path . $!";
my %element_char_hash;

while ( my $line = <$TE_in> ) {
  chomp $line;
  my @split = split( "\t", $line );
  my $ele_name = $split[0];
  if ( $split[1] ne "." ) {
    my $seq     = $split[1];
    my $seq_obj = Bio::Seq->new( -seq => $seq, -alphabet => 'dna' );
    my $iupac   = Bio::Tools::IUPAC->new( -seq => $seq_obj );
    my $regexp  = $iupac->regexp();
    $element_char_hash{$ele_name}{'tir_con'} = $regexp;
  }
  else {
    $element_char_hash{$ele_name}{'tir_con'} = '';
  }
  if ( $split[2] =~ m/,/ ) {
    my @tsd_lengths   = split( ',', $split[2] );
    my $tsd_array_len = @tsd_lengths;
    my $count         = 1;
    my $regexp        = '';
    foreach my $tsd_len (@tsd_lengths) {
      if ( $count < $tsd_array_len ) {
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
  if ( $split[3] ne "." ) {
    if ( $split[3] =~ m/,/ ) {
      my @tsd_cons     = split( ',', $split[3] );
      my $tsd_cons_len = @tsd_cons;
      my $count        = 1;
      my $regexp       = '';
      foreach my $tsd_consensus (@tsd_cons) {
        my $seq_obj =
          Bio::Seq->new( -seq => $tsd_consensus, -alphabet => 'dna' );
        my $iupac = Bio::Tools::IUPAC->new( -seq => $seq_obj );
        my $regexp_ori = $iupac->regexp();
        if ( $count < $tsd_cons_len ) {
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
      my $seq_obj = Bio::Seq->new( -seq => $split[3], -alphabet => 'dna' );
      my $iupac = Bio::Tools::IUPAC->new( -seq => $seq_obj );
      my $regexp = $iupac->regexp();
      $element_char_hash{$ele_name}{'tsd_con'} = $regexp;
    }
  }
  else {
    $element_char_hash{$ele_name}{'tsd_con'} = '';
  }
}
print "Checking element for DNA TE charateristics\n";
my %element_hits;
foreach my $ele_name ( keys %element_char_hash ) {
  if (
    $element_info{'TSD_len'} =~ m/$element_char_hash{$ele_name}{"tsd_length"}/ )
  {
    $element_hits{$ele_name}++;
    if (  $element_char_hash{$ele_name}{"tsd_con"} ne ''
      and $element_info{'TSD_seq'} =~
      m/$element_char_hash{$ele_name}{"tsd_con"}/i )
    {
      $element_hits{$ele_name}++;
    }
    if ( $element_char_hash{$ele_name}{"tir_con"} ne '' ) {
      if ( $element_info{'left_tir_seq'} =~
        m/^$element_char_hash{$ele_name}{"tir_con"}/i )
      {
        $element_hits{$ele_name}++;
      }

      else {
        delete $element_hits{$ele_name};
      }
    }
  }
}

my @sorted;
foreach my $key (
  sort { $element_hits{$b} <=> $element_hits{$a} }
  keys(%element_hits)
  )
{
  my $count = $element_hits{$key};
  push @sorted, [ $key, $count ];
}

my $classification = '';

#store all classifications and the number of hits to each
foreach my $row_ref (@sorted) {
  my @info = @{$row_ref};
  if ( $info[1] == ${ $sorted[0] }[1] ) {
    if ( $classification eq '' ) {
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
if ( $classification eq '' ) {
  $classification = "Unknown";
}
$element_info{"classification"} = $classification;
print "Element classification finished\n";
my $element_info_out_path =
  File::Spec->catpath( $volume, $out_path, $filename . ".element_info" );
open( my $element_info_out, ">", $element_info_out_path )
  or die "Error creating $element_info_out_path. $!\n";
print $element_info_out join( "\t",
  $fname_fin,                   $element_info{'copy_num'},
  $element_id,                  $element_info{'left_tir_seq'},
  $element_info{'left_tir_id'}, $element_info{'right_tir_seq'},
  $element_info{'left_tir_id'}, $element_info{'TSD_len'},
  $element_info{'TSD_seq'},     $element_info{'TSD_fraction'},
  $classification ),
  "\n";
close($element_info_out);

print "Printing fasta\n";
my $element_fasta_out_path =
  File::Spec->catpath( $volume, $out_path, $filename . ".fasta" );
print_fasta( $element_fasta_out_path, $element_aln_obj );

my $out_fix           = $out_path . "/";
my $Blogo_config_path = $FindBin::Bin . "/../lib/blogo/Blogo.conf";
my $Blogo_path        = $FindBin::Bin . "/../lib/blogo";
system(
"$Blogo_path/Blogo_batch.pl file_path=$out_fix file_names=$insertion_site_file img_abs_dir=$out_fix conf_file=$Blogo_config_path"
);

print "\nSetting up gff path\n";
my $gff_path = File::Spec->catpath( $volume, $out_path, $filename . ".gff" );
print "Calling gff printer\n";
generate_gff( $final_aln_obj, $gff_path, "final", \%element_info );
print "Exited gff printer\n";

if ( !defined $all ) {
  print "Cleaning up files\n";
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
  if ( !-f $input_path ) {
    die "Supplied filepath is not valid";
  }
  my $seq = $self->seq();
  $seq =~ s/-//g;
  my $seq_name = $self->id();
  my $fa_aln_obj =
    Bio::SearchIO->new( -format => 'fasta', -file => $input_path );
  my @result;
  print "ggsearch results path: $input_path\n";

  #go through the FASTA input object . . down to the HSP
  while ( my $result = $fa_aln_obj->next_result ) {

    #print "numHits: ",$result->num_hits,"\n";
    if ( $result->num_hits == 0 ) {
      push @result, ( 0, [0] );
      last;
    }
    while ( my $hit = $result->next_hit ) {
      while ( my $hsp = $hit->next_hsp ) {

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
        my $match_cutoff = 8;

    #parse homology string, keeping track of match length and mismatches or gaps
        for ( my $count = 0 ; $count < length($homo_string) ; $count++ ) {
          my $homo_char  = substr( $homo_string, $count, 1 );
          my $query_char = substr( $query_str,   $count, 1 );
          my $hit_char   = substr( $hit_str,     $count, 1 );
          if ( $count == 8 and $total_mis_aln >= 5 ) {
            if ( $match_len < 3 ) {
              $match_len   = 0;
              $start_pos   = '';
              $match_query = '';
              $match_hit   = '';
              $end_pos     = '';
              print
"No TIRs found near start of sequences, resetting counts and ending\n";
              last;
            }
          }
          if ( $round == 2 ) {
            if ( $count == 6 and $total_mis_aln >= 4 ) {
              $match_len   = 0;
              $start_pos   = '';
              $match_query = '';
              $match_hit   = '';
              $end_pos     = '';
              last;
            }
          }
          ## skip any seqs that have more than 2 mismatches in the first 3 bases of the TIR
          if ( $round == 3 or $round == 1 ) {
            if ( $count == 3 and $total_mis_aln >= 2 ) {
              $match_len   = 0;
              $start_pos   = '';
              $match_query = '';
              $match_hit   = '';
              $end_pos     = '';
              last;
            }
          }

          if ( $match_len == 0 ) {

#if match length equals 0 and position is not a match, continue to next position
            if ( $homo_char eq " " ) {
              $total_mis_aln++;
              next;
            }

            #if position is a match, store info, continue to next position
            elsif ( $homo_char eq ":" ) {
              $start_pos = $count;
              $last_good = $count;
              $match_len++;
              $match_query .= $query_char;
              $match_hit   .= $hit_char;

              #print "Initial match at $start_pos\n";
              next;
            }
          }
          elsif ( $match_len >= 1 and $match_len < ( $match_cutoff - 1 ) ) {

#if match length is 1-3 and position is not a match, increment mismatch counter and check if more than one mismatch has occurred
            if ( $homo_char eq " " ) {
              $match_mis_aln++;
              $total_mis_aln++;

              #allow one mismatch, store info and continue
              if ( $match_mis_aln <= 1 ) {
                $match_len++;
                $last_good = $count;
                $match_query .= $query_char;
                $match_hit   .= $hit_char;

                # print "First Mismatch at $count\n";
                next;
              }

              #more than one mismatch, reset counters and other info, continue
              elsif ( $match_mis_aln > 1 and $match_len < 5 ) {
                $match_len     = 0;
                $start_pos     = '';
                $match_query   = '';
                $match_hit     = '';
                $match_mis_aln = 0;

                #print "Another Mismatch at $count, resetting counts\n";
                next;
              }
              elsif ( $match_mis_aln < 3 and $match_len >= 5 ) {
                $match_len++;
                $last_good = $count;
                $match_query .= $query_char;
                $match_hit   .= $hit_char;

                # print "Another Mismatch at $count\n";
                next;
              }
              elsif ( $total_mis_aln >= 3 ) {
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
            elsif ( $homo_char eq ":" ) {
              $last_good = $count;
              $match_len++;
              $match_query .= $query_char;
              $match_hit   .= $hit_char;

              #print "Another match at $count. Length is $match_len\n";
              next;
            }
          }
          elsif ( $match_len >= $match_cutoff - 1 ) {

#match length is 5 or higher. If position is not a match, increment mismatch counter and check if more than 2 mismatches have occurred. If a match, continue.
            if ( $homo_char eq " " ) {
              $match_mis_aln++;
              $total_mis_aln++;

              #mismatches under 3, store info and continue
              if ( $match_mis_aln <= 3 ) {
                $match_len++;
                $last_good = $count;
                $match_query .= $query_char;
                $match_hit   .= $hit_char;

                # print "Another Mismatch at $count, proceeding\n";
                next;
              }

     #mismatches 3 or more, store final info for alignment match and end parsing
              elsif ( $match_mis_aln >= 3 ) {
                $end_pos = $last_good;
                $match_query =~ s/-//g;
                $match_hit   =~ s/-//g;

                #reverse complement the match query sequence
                $match_query =~ tr/ATGC/TACG/;
                $match_query = reverse($match_query);
                my $match_query_len = length($match_query);
                my $match_hit_len   = length($match_hit);

    #find the position in the full sequence of the hit and query match sequences
                $hit_pos = index( uc($seq), uc($match_hit) ) + 1;
                #$hit_pos = index( uc($seq), uc($match_hit), 40 ) + 1;
                $query_pos = rindex( uc($seq), uc($match_query) ) + 1;

          #print "$seq_name hit_pos:$hit_pos query_pos:$query_pos $match_hit\n";
          #store sequence name and the hit and query info
                my @match = (
                  $seq_name,
                  {
                    "hit"   => [ $hit_pos,   $match_hit_len ],
                    "query" => [ $query_pos, $match_query_len ]
                  }
                );
                ## add a catch for $hit_pos or $query_pos == 0 [index returns -1 when $match_hit was not found in $seq]
                push @result, ( 1, [@match] );

#print "Another Mismatch at $count. Match is long, pushing match info for output\n";
                last;
              }
            }

            #position is a match, store info and continue
            elsif ( $homo_char eq ":" ) {
              $last_good = $count;
              $match_len++;
              $match_query .= $query_char;
              $match_hit   .= $hit_char;

              #print "Another match at $count. Length is $match_len\n";
              next;
            }
          }
        }

        #add in check to see if TIRs were found.
        if ( $end_pos eq '' ) {
          push @result, ( 0, [0] );
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
  print "Entered .gff printer\n";
  $self->throw("Need Bio::Align::AlignI argument")
    unless ref $self && $self->isa('Bio::Align::AlignI');

  #not implemented: round changes the output based on when the call is made
  if ( $round eq 'final' or $round eq 'TSD' ) {
    open( my $out, ">", $path ) or die "Error creating $path. $!\n";
    print "Round = final\n";
    foreach my $seq_obj ( $self->each_seq() ) {
      my $seq_name = $seq_obj->id();
      my $seq      = $seq_obj->seq();
      $seq =~ s/-//g;
      my $seq_len    = length($seq);
      my $ori_end    = $seq_len - $flank;
      my $left_comp  = $ele_info_ref->{$seq_name}{"left_tir_start"} - 101;
      my $right_comp = $ele_info_ref->{$seq_name}{"right_tir_start"} - $ori_end;
      my $copy_num;
      my $eleid;
      my $seqid;
      my $type = 'terminal_inverted_repeat_element';
      my $start;
      my $end;
      my $strand;

      #grab copy information from TARGeT output
      if ( $seq_name =~
/^([0-9]*).+_Query:(.*)_Sbjct:(.*)_Length.+Location:\(([0-9]*)_\-_([0-9]*)\)_Direction:(.+)/
        )
      {
        $copy_num = $1;
        $eleid    = $2;
        $seqid    = $3;
        $start    = $4 + $left_comp;
        $end      = $5 + $right_comp;
        $strand   = $6;
      }

      #grab copy information from RSPB output
      elsif ( $seq_name =~ /^([0-9]*)_(.+)_(.+)-(.+)_(.+)/ ) {
        $copy_num = $1;
        $seqid    = $2;
        $start    = $3 + $left_comp;
        $end      = $4 + $right_comp;
        $eleid    = $5;
        $strand   = "?";
      }
      else {
        print "Header doesn't match TARGeT or RSPB:  $seq_name\n";
        next;
      }
      if ( $eleid =~ m/(.+)_TSD/ or $eleid =~ m/(.+)_Unknow/ ) {
        $eleid = $1;
      }
      my $ltir_end = $start + length( $ele_info_ref->{"left_tir_seq"} ) - 1;
      my $rtir_start =
        $end - ( length( $ele_info_ref->{"right_tir_seq"} ) - 1 );
      my $ele_id    = $ele_info_ref->{"element_id"};
      my $tir_id    = $ele_info_ref->{"left_tir_id"};
      my $ele_class = 'N/A';
      my $tsd_frac  = 'N/A';
      my $tsd_con   = 'N/A';

      if ( $round eq 'final' ) {
        $ele_class = $ele_info_ref->{"classification"};
        $tsd_frac  = $ele_info_ref->{"TSD_fraction"};
        $tsd_con   = $ele_info_ref->{"TSD_seq"};
      }
      print $out join(
        "\t", $seqid,
        $PROGRAM_NAME,
        $type, $start, $end, '.', $strand, '.',
        join( ";",
          "ID=$eleid-$copy_num", "Name=$eleid Copy$copy_num",
          "element_id=$ele_id",  "element_classification=$ele_class",
          "tir_id=$tir_id",      "tsd_fraction=$tsd_frac",
          "tsd_consensus=$tsd_con" )
        ),
        "\n";
      print $out join( "\t",
        $seqid, $PROGRAM_NAME, "five_prime_terminal_inverted_repeat",
        $start, $ltir_end, ".", ".", ".", "Parent=$eleid-$copy_num" ),
        "\n";
      print $out join( "\t",
        $seqid, $PROGRAM_NAME, "three_prime_terminal_inverted_repeat",
        $rtir_start, $end, ".", ".", ".", "Parent=$eleid-$copy_num" ),
        "\n";

      if ( $round eq 'final' and defined $ele_info_ref->{$seq_name}{"TSD"} ) {
        my $ltsd_start = $start - length( $ele_info_ref->{$seq_name}{"TSD"} );
        my $ltsd_end   = $start - 1;
        my $rtsd_start = $end + 1;
        my $rtsd_end   = $end + length( $ele_info_ref->{$seq_name}{"TSD"} );
        print $out join( "\t",
          $seqid, $PROGRAM_NAME, "target_site_duplication",, $ltsd_start,
          $ltsd_end, ".", ".", ".", "Derives_from=$eleid-$copy_num" ),
          "\n";
        print $out join( "\t",
          $seqid, $PROGRAM_NAME, "target_site_duplication", $rtsd_start,
          $rtsd_end, ".", ".", ".", "Derives_from=$eleid-$copy_num" ),
          "\n";
      }
    }
    close($out);
  }
}

sub clean_files {
  my $out_path = shift;
  opendir( my $in_DIR, $out_path ) or die "Cannot open directory: $!";
  while ( my $file = readdir($in_DIR) ) {
    next if ( $file =~ m/^\./ );
    if ( $file =~ m/\.(final|info|fa|element_info|bad|gff|tif|jpg|full_id)$/ ) {
      next;
    }
    else {
      my $remove_path = File::Spec->catpath( $volume, $out_path, $file );
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
  my ( $left_tir_start, $left_tir_end, $right_tir_start, $right_tir_end ) =
    ( 0, 0, 0, 0 );

  foreach my $seq_obj ( $self->each_seq() ) {
    my $seq_name = $seq_obj->id();

    #print $seq_name,"\n";
    #skip sequences not in the hash
    if ( exists $tir_positions->{$seq_name} ) {

      #print "For get_columns: Entry exist in hash -", $seq_name, "\n";
      #print Dumper($tir_positions->{$seq_name});
      #print "Left start: $tir_positions->{$seq_name}{'left_tir_start'}\n";
      $left_tir_start =
        $self->column_from_residue_number( $seq_name,
        $tir_positions->{$seq_name}{'left_tir_start'} );
      $right_tir_start =
        $self->column_from_residue_number( $seq_name,
        $tir_positions->{$seq_name}{'right_tir_start'} );
      if (  $round == 2
        and !exists $tir_positions->{$seq_name}{'left_tir_end'}
        and !exists $tir_positions->{$seq_name}{'right_tir_end'} )
      {
        delete $tir_positions->{$seq_name};
      }
      elsif ( $round == 2 ) {
        ##print "Round = $round\n";
        $left_tir_end =
          $self->column_from_residue_number( $seq_name,
          $tir_positions->{$seq_name}{'left_tir_end'} );
        $right_tir_end =
          $self->column_from_residue_number( $seq_name,
          $tir_positions->{$seq_name}{'right_tir_end'} );
        last;
      }
    }
  }

#print "in get_col subrout: left tir start: $left_tir_start, right_tir_start: $right_tir_start, left_tir_end: $left_tir_end, right_tir_end: $right_tir_end\n";
  return ( $left_tir_start, $right_tir_start, $left_tir_end, $right_tir_end );
}

sub adjust_tir_starts {
### adjust tir starts.
  print "Attempting to adjust tir starts\n";
  my $aln_obj         = shift;
  my $right_tir_start = shift;
  my $left_tir_start  = shift;

  my $num_seqs = $aln_obj->num_sequences;
  my %count;
  my %modified;
  my $orig_right_tir_start = $right_tir_start;
  my $orig_left_tir_start  = $left_tir_start;
  foreach my $seq_obj ( $aln_obj->each_seq ) {
    my $seq    = $seq_obj->seq;
    my $seq_id = $seq_obj->id;
    my ( $gaps_left, $gaps_right ) = ( '', '' );
    ## for right TIR
    my $nt_col_pos = $right_tir_start;
    my $left2right = $seq_obj->subseq( $left_tir_start, $nt_col_pos );
    my $next2end   = $seq_obj->subseq( $nt_col_pos + 1, length $seq );
    if ( $next2end =~ /^(-+)/ ) {
      $gaps_right = $1;
      $next2end =~ s/^(-+)//;
    }
    my @next2end = split '', $next2end;
    for ( my $i = 0 ; $i < @next2end ; $i++ ) {
      my $nt = $next2end[$i];
      $count{right}{$i}{$nt}++;
    }

    ## for left_TIR
    $nt_col_pos = $left_tir_start;
    my $start2left = $seq_obj->subseq( 1, $nt_col_pos - 1 );

    #print "$seq_id\n";
    #print "before: $start2left\n";
    if ( $start2left =~ /(-+)$/ ) {
      $gaps_left = $1;
      $start2left =~ s/(-+)$//;
      #print "after:  $start2left\n";
    }
    my @gaps = split '', $gaps_left;
    my @start2left = split '', $start2left;
    unshift @start2left , @gaps;
    for ( my $i = ( (scalar @start2left ) - 1 ) ; $i >= 0 ; $i-- ) {
      my $nt = $start2left[$i];
      $count{left}{$i}{$nt}++;
      #print substr($seq_id,0,5) ," $i $nt " , $count{left}{$i}{$nt} , "\n";
    }
    my $new_seq =
      $gaps_left . $start2left . $left2right . $next2end . $gaps_right;
    $modified{$seq_id} = $new_seq;

#my $new_seq_p = "$gaps_left . $start2left . $left2right . $next2end . $gaps_right";
#print "new: $new_seq_p\n";
  }

  #warn Dumper \%count;
  my $done = 0;
  ## for right tir
  foreach my $pos ( sort { $a <=> $b } keys %{ $count{right} } ) {
    my $nt = (
      sort { $count{right}{$pos}{$b} <=> $count{right}{$pos}{$a} }
        keys %{ $count{right}{$pos} }
    )[0];
    my $nt_count   = $count{right}{$pos}{$nt};
    my $percent_nt = $nt_count / $num_seqs;

#print "adjust_tirs RT $nt($pos):nt_count (nextposIN:$count{left}{$pos-1}) %=$percent_nt\n";
    if ( $percent_nt >= .8 ) {
      $right_tir_start++;
    }
    else {
      $done = 1;
    }
    last if $done;
  }
  $done = 0;
  ## for left tir
  foreach my $pos ( sort { $b <=> $a } keys %{ $count{left} } ) {
    my $nt = (
      sort { $count{left}{$pos}{$b} <=> $count{left}{$pos}{$a} }
        keys %{ $count{left}{$pos} }
    )[0];
    my $nt_count   = $count{left}{$pos}{$nt};
    my $percent_nt = $nt_count / $num_seqs;

#print "adjust_tirs LT $nt($pos):nt_count  ($nt_count/$num_seqs)=$percent_nt%\n";
    if ( $percent_nt >= .8 ) {
      $left_tir_start--;
    }
    else {
      $done = 1;
    }
    last if $done;
  }

  if ( $orig_right_tir_start != $right_tir_start
    or $orig_left_tir_start != $left_tir_start )
  {
    print "origR: $orig_right_tir_start\n";
    print "after adjR: $right_tir_start\n";
    print "origL: $orig_left_tir_start\n";
    print "after adjL: $left_tir_start\n";
    foreach my $seq_obj ( $aln_obj->each_seq ) {
      ## replace org seq with mod seq
      my $seq_id = $seq_obj->id;
      if ( exists $modified{$seq_id} ) {
        my $new_seq = $modified{$seq_id};
        $seq_obj->seq($new_seq);
      }
    }
  }
  return ( $right_tir_start, $left_tir_start );
}

sub consensus_filter {
#### start
  my $gap_seq_pos_remove = shift;              ## \@gap_seq_pos_remove
  my $aln_obj            = shift;
  my $left_tir_start     = shift;
  my $right_tir_start    = shift;
  my $tir_positions      = shift;              # \%tir_positions
  my $round              = shift;
  my $aln_len            = $aln_obj->length;
  $round = defined $round ? $round : 'other';
  print "in consensu_filt sub: round=$round\n";
  my %trim_gap_seq_remove;
  ## this will generate a sequence with '?' at every position in which there is less
  ## than 80% consensus
  my $consensus = $aln_obj->consensus_string(80);
  print "$consensus\n";
  ## first round of tir finding
  if ( $left_tir_start == 0 and $right_tir_start == 0 ) {
    my @consensus    = split '', $consensus;
    my $nt_count     = 0;
    my $last_nomatch = 0;

    ## left tir
    for ( my $i = 0 ; $i < ( ( scalar @consensus ) * .5 ) ; $i++ ) {
      my $nt = $consensus[$i];
      if ( $nt eq '?' and $nt_count < 3 ) {
        $last_nomatch = $i;
        $nt_count     = 0;
      }
      elsif ( $nt ne '?' and $nt_count < 3 ) {
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
    for (
      my $i = ( ( scalar @consensus ) - 1 ) ;
      $i > ( ( scalar @consensus ) * .5 ) ;
      $i--
      )
    {
      my $nt = $consensus[$i];
      if ( $nt eq '?' and $nt_count < 3 ) {
        $last_nomatch = $i;
        $nt_count     = 0;
      }
      elsif ( $nt ne '?' and $nt_count < 3 ) {
        $nt_count++;
      }
      else {    #($nt ne '?' and $nt_count => 3){
        $right_tir_start = $last_nomatch;
        last;
      }
    }
  }
  print "in con_fil sub: leftTIR: $left_tir_start\n";
  print "in con_fil sub: rightTIR: $right_tir_start\n";
  my $aln_depth = $aln_obj->num_sequences;

  #print "consen: " ,$consensus , "\n";
  my %con_mm;
  my %to_remove;
  for ( my $i = 0 ; $i < $aln_len ; $i++ ) {
    if ( ( $i >= $left_tir_start and $i <= $left_tir_start + 20 )
      or ( $i >= $right_tir_start - 20 and $i <= $right_tir_start ) )
    {

      ## check each col for seqs that have a run of too many mismatches
      my $col_aln_obj = $aln_obj->slice( $i + 1, $i + 1, 1 );
      foreach my $seq_obj ( $col_aln_obj->each_seq() ) {
        my $seq_name = $seq_obj->id;
        my $nt       = $seq_obj->seq;
        my $con_nt   = substr( $consensus, $i, 1 );
        next if $con_nt =~ /N/i;
        if ( $nt ne $con_nt ) {
          $con_mm{$seq_name}++;
        }
      }
      foreach my $key ( keys %con_mm ) {
        next if $con_mm{$key} < 5;
        my $seq_obj_to_remove = $aln_obj->get_seq_by_id($key);
        $trim_gap_seq_remove{$key} = $seq_obj_to_remove;
        my @info = ( $key, $i + 1, $seq_obj_to_remove );
        push @$gap_seq_pos_remove, [@info];
        $to_remove{$key} = 1;
      }
      ## end check each col for too many mismatches
    }
  }
  my $to_remove_count = keys %to_remove;
  print "in con_fil sub before get_tir_nt_starts: leftTIR: $left_tir_start\n";
  print "in con_fil sub before get_tir_nt_starts: rightTIR: $right_tir_start\n";
  ## skip seq removal if the number to remove is about the depth of the aln
  print "there are $aln_depth seqs in aln before cons_filter\n";
  print "and there will be  $aln_depth - $to_remove_count = ",
    $aln_depth - $to_remove_count, " after filtering\n";
  if ( $aln_depth > ( $to_remove_count + 5 ) ) {
    print "removing seqs with cons_filter\n";
    foreach my $key ( keys %trim_gap_seq_remove ) {
      my $seq_obj = $trim_gap_seq_remove{$key};
      $aln_obj->remove_seq($seq_obj);
    }

  }## paren for remove least
  ##how close are the first round of TIR starts to the ends of the aln
  my $consensus_len = $right_tir_start - $left_tir_start + 1;
  #print "consensus string (len=$consensus_len) is ",$consensus_len/$aln_len," of the len of the aln (len=$aln_len)\n"; 
  my $message = '';
  if ($consensus_len > ($aln_len*.9)){ ## ($consensus_len > ($aln_len - ($flank*2) + 50 ))
    $message = "$filename: too much of the alignment is conserved. Flanks are too similar\n";
  }elsif (!$left_tir_start or !$right_tir_start){
    $message = "$filename: no TIR start was found on one or both ends\n";
  }
  if ($message){
    my $bad_out_path =
      File::Spec->catpath( $volume, $out_path, $filename . ".bad" );
    error_out( $bad_out_path, $message );
  }
  my $ref2remove_these;
  if ( $round ne 'final' ) {
    ($tir_positions,$ref2remove_these) =
      get_tir_nt_starts( $aln_obj, $tir_positions, $left_tir_start,
      $right_tir_start );
  }
  foreach my $key ( keys %$ref2remove_these ) {
    my $seq_obj = $$ref2remove_these{$key};
    $aln_obj->remove_seq($seq_obj);
    my @info = ( $key, 0, $seq_obj );
    push @$gap_seq_pos_remove, [@info];
  } 

    $aln_obj = $aln_obj->remove_gaps( '-', 1 );
    ( $left_tir_start, $right_tir_start ) =
      get_columns( $aln_obj, $tir_positions, 1 );
  #}    ##paren for remove most
  print "in con_fil sub after getCol: leftTIR: $left_tir_start\n";
  print "in con_fil sub after getCol: rightTIR: $right_tir_start\n";
  return ( $left_tir_start, $right_tir_start, $gap_seq_pos_remove, $aln_obj,
    $tir_positions );
}

sub gap_filter {
  my $gap_seq_pos_remove = shift;              ## \@gap_seq_pos_remove
  my $aln_obj            = shift;
  my $left_tir_start     = shift;
  my $right_tir_start    = shift;
  my $aln_len            = $aln_obj->length;

  my %trim_gap_seq_remove;
  for ( my $i = 0 ; $i < $aln_len ; $i++ ) {
    if ( ( $i >= $left_tir_start - 1 and $i <= $left_tir_start + 19 )
      or ( $i >= $right_tir_start - 21 and $i <= $right_tir_start - 1 ) )
    {

      my $gap_cols        = $aln_obj->gap_col_matrix();
      my @gap_col_array   = @{$gap_cols};
      my $gap_col_hashref = $gap_col_array[$i];
      my %gap_col_hash    = %{$gap_col_hashref};
      my $base_count      = 0;
      my $total_count     = 0;
      foreach my $key ( keys %gap_col_hash ) {
        if ( $gap_col_hash{$key} != 1 ) {
          $base_count++;
        }
        $total_count++;
      }
      my $present_fraction = $base_count / $total_count;
      ## if half or more of the seqs in this col are gaps, mark seqs with a gaps
      ## to be be removed
      if ( $present_fraction < 0.5 ) {
        foreach my $key ( keys %gap_col_hash ) {
          if ( $gap_col_hash{$key} != 1 ) {
            my $seq_obj = $aln_obj->get_seq_by_id($key);
            $trim_gap_seq_remove{$key} = $seq_obj;
            my @info = ( $key, $i + 1, $seq_obj );
            push @$gap_seq_pos_remove, [@info];
          }
        }
      }
      ## if more than half of the seqs are not gaps, but one seq has a gap,
      ## mark it for removal
      else {
        foreach my $key ( keys %gap_col_hash ) {
          if ( $gap_col_hash{$key} == 1 ) {
            my $seq_obj = $aln_obj->get_seq_by_id($key);
            my $seq     = $seq_obj->seq();
            my $seq_pos = $seq_obj->location_from_column( $i + 1 );
            $trim_gap_seq_remove{$key} = $seq_obj;
            my @info = ( $key, $i + 1, $seq_obj );
            push @$gap_seq_pos_remove, [@info];
          }
        }
      }
    }
  }

  foreach my $key ( keys %trim_gap_seq_remove ) {
    my $seq_obj = $trim_gap_seq_remove{$key};
    $aln_obj->remove_seq($seq_obj);
  }
  return ( $gap_seq_pos_remove, $aln_obj );
}

sub get_tir_nt_starts {
  my $aln_obj         = shift;
  my $tir_positions   = shift;    # \%tir_positions
  my $left_tir_start  = shift;
  my $right_tir_start = shift;
  my $remove_these;
  print "get_tir_nt_starts(top): lts:$left_tir_start rts:$right_tir_start\n";
  foreach my $seq_obj ( $aln_obj->each_seq() ) {
    my $seq_name            = $seq_obj->id();
    my $seq                 = $seq_obj->seq();
    if ($seq =~ /^-+$/g){
      #print "$seq_name: no seq found\n";
      $$remove_these{$seq_name}=$seq_obj;
    }
    my $left_tir_start_pos = 0;
    if (defined $seq_obj->location_from_column($left_tir_start) ){
      my $left_tir_start_obj  = $seq_obj->location_from_column($left_tir_start);
      $left_tir_start_pos  = $left_tir_start_obj->start();
    }
    $tir_positions->{$seq_name}{'left_tir_start'}  = $left_tir_start_pos;
    
    my $right_tir_start_pos = 0;
    if (defined $seq_obj->location_from_column($right_tir_start)){ 
      my $right_tir_start_obj = $seq_obj->location_from_column($right_tir_start);
      $right_tir_start_pos = $right_tir_start_obj->start();
    }
    $tir_positions->{$seq_name}{'right_tir_start'} = $right_tir_start_pos;
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

  foreach my $seq_obj ( $trimmed_aln_obj->each_seq() ) {
    my $seq_name           = $seq_obj->id();
    my $trim_seq           = $seq_obj->seq();
    my $left_tir_start_obj = $seq_obj->location_from_column($left_tir_start);
    my $left_tir_start_pos = $left_tir_start_obj->start();
    my $left_tir_end_pos   = $left_tir_start_pos + $left_tir_len - 1;
    my $right_tir_start_obj =
      $seq_obj->location_from_column(
      $right_tir_start + ( $right_tir_len - 1 ) );
    my $right_tir_start_pos = $right_tir_start_obj->start();
    my $right_tir_end_pos = $right_tir_start_pos - ( $right_tir_len - 1 );
    $tir_positions->{$seq_name}{'left_tir_start'}  = $left_tir_start_pos;
    $tir_positions->{$seq_name}{'left_tir_end'}    = $left_tir_end_pos;
    $tir_positions->{$seq_name}{'right_tir_start'} = $right_tir_start_pos;
    $tir_positions->{$seq_name}{'right_tir_end'}   = $right_tir_end_pos;
  }
  return ($tir_positions);
}

sub error_out {
  my $bad_out_path = shift;    ## error file path
  my $message      = shift;
  my $exit         = shift;    ## 1=no exit; 0 or undef=exit
  open( my $bad_out, ">", $bad_out_path );
  print $bad_out "$message\n";
  close($bad_out);
  if ( !defined $all ) {
    print "Cleaning up files\n";
    clean_files($out_path);
  }
  if ( !$exit ) {
    exit 0;
  }
}

sub print_fasta {
  my $filename = shift;
  my $aln_obj  = shift;
  my $seqIO_out_obj =
    Bio::SeqIO->new( -format => 'fasta', -file => ">$filename" );
  foreach my $seq_obj ( $aln_obj->each_seq ) {
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
  open MSAOUT, ">$msa.mod" or die "Can't open $msa.mod\n";
  <MSA>;    #throw out first header
  my $seq;
  my $seq_count = 0;
  while ( my $line = <MSA> ) {
    chomp $line;
    if ( $line =~ /^>/ ) {
      $seq_count++;
      print MSAOUT "$seq\n";
      $seq = '';
    }
    else {
      $seq .= $line;
    }
  }
## for last seq
  $seq_count++;
  print MSAOUT "$seq\n";
##

  close MSAOUT;
  close MSA;
  my @percent_id;
  my $first_seq = `head -1 $msa.mod`;
  if (length $first_seq < 5){
   die "Error retrieving first line of $msa.mod\n";
  } 
  chomp $first_seq;
  my $len = length $first_seq;
  open OUT, ">$out" or die "Can't opne $out $!\n";
  for ( my $i = 1 ; $i < $len + 1 ; $i++ ) {
    my $total_nt;
    my %nt_count;
    chomp( my @col = `cut -c $i $msa.mod` );
    foreach my $nt (@col) {
      $nt_count{total}++;
      if ( $nt =~ /-/ ) {
        $nt_count{dash}++;
        next;
      }
      ## only increment if the character is not a '-'
      ## Ns will be included in the total
      #$nt_count{total}++;
      #next if $nt =~ /N/i;
      $nt_count{each}{$nt}++;
    }
    if ( ( scalar keys %{ $nt_count{each} } ) == 0 ) {
      push @percent_id, [ $i, 0, 0 ];
      print OUT join( "\t", $i, 0, 0 ), "\n";
    }
    else {
      my $most_freq_nt =
        ( sort { $nt_count{each}{$b} <=> $nt_count{each}{$a} }
          keys %{ $nt_count{each} } )[0];
      if ( $most_freq_nt =~ /N/i ) {
        push @percent_id, [ $i, 0, 0 ];
        print OUT join( "\t", $i, 0, 0 ), "\n";
      }
      my $count      = $nt_count{each}{$most_freq_nt};
      my $col_total  = $nt_count{total};
      my $dash_total = exists $nt_count{dash} ? $nt_count{dash} : 0;
      ## rigth now, col_total is eq to seq_count
      my $pid = ( $count / $col_total ) * 100;

      #push @percent_id, ($count/$col_total);
      push @percent_id,
        [ $i, $pid, ( ( $col_total - $dash_total ) / $seq_count ) ];

      # print out this info to a file
      print OUT
        join( "\t", $i, $pid, ( ( $col_total - $dash_total ) / $seq_count ) ),
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
  for ( my $i = 1 ; $i < $aln_len ; $i++ ) {
    my $id_row_ref      = $$id_array[$i];
    my $gap_col_hashref = $gap_cols->[$i];
    my %gap_col_hash    = %{$gap_col_hashref};
    ## remove entire seq if %gaps in column is >= 80
    ## and (%id <= 50 or %nt_aligned <= 10%)
    if ( $gap_id_array[$i] >= 80.0
      and ( $id_row_ref->[1] <= 50 or $id_row_ref->[2] <= .1 ) )
    {
      foreach my $key ( keys %gap_col_hash ) {
        ## value of 1 means gap
        ## if there is no gap at this position where there are may gaps in the column
        ## we will make a note to remove this seq
        if ( $gap_col_hash{$key} != 1 ) {
          my $seq_obj = $aln_obj->get_seq_by_id($key);
          my @info = ( $key, $i + 1, $seq_obj );
          push @$ref2gspr, [@info];
          $$ref2gsr{$key} = $seq_obj;
        }
      }
    }
  }
## removes any sequences found in last step
  foreach my $key ( keys %$ref2gsr ) {
    my $seq_obj = $$ref2gsr{$key};
    $aln_obj->remove_seq($seq_obj);
  }
## removes any columns that are now all gaps
  $aln_obj = $aln_obj->remove_gaps( '-', 1 );

  my $test_aln_out =
    File::Spec->catpath( $volume, $out_path, $filename . ".removeMost_trim_0" );
  $out = Bio::AlignIO->new(
    -file             => ">$test_aln_out",
    -format           => 'fasta',
    -displayname_flat => 0
  );
  $out->write_aln($aln_obj);

  my ( $left_tir_start, $right_tir_start );
  ( $left_tir_start, $right_tir_start, $ref2gspr, $aln_obj, $tir_positions ) =
    consensus_filter( $ref2gspr, $aln_obj, 0, 0, $tir_positions );

  #@gap_seq_pos_remove = @{$ref2gspr};

  print
    "after cons_filter (filename.trim: $left_tir_start, $right_tir_start\n";

#my ($left_tir_start1,$right_tir_start1,$tmp_aln_obj,$ref2tp,$ref2gsr, $ref2gspr) = remove_most ($full_aln_obj,\%tir_positions, \@full_id_array);
  return ( $left_tir_start, $right_tir_start, $aln_obj, $tir_positions,
    $ref2gsr, $ref2gspr );

} ## end remove_most

sub remove_least {

#my ($left_tir_start1,$right_tir_start1,$tmp_aln_obj,$ref2tp,$ref2gsr, $ref2gspr) = remove_most ($full_aln_obj,\%tir_positions, \@full_id_array);
  my $org_aln_obj   = shift;
  my $tir_positions = shift;
  my $id_array      = shift;
  my $aln_obj       = $org_aln_obj;
  my $aln_len       = $aln_obj->length;
  my $ref2gsr;    ## %gap_seq_remove

  foreach my $arrayref (@$id_array) {
    my $col_id = $$arrayref[0];
    my $pos_id = $$arrayref[1];
    if ( $first_col_80 == 1 and $pos_id >= 80 ) {
      $first_col_80 = $col_id;
    }
    elsif ( $pos_id >= 80 ) {
      $last_col_80 = $col_id;
    }
  }

  print "first_col_80:$first_col_80 last_col_80:$last_col_80\n";
  for ( my $i = $first_col_80 ; $i < $last_col_80 ; $i++ ) {
    next
      if $i == $first_col_80 + ( $aln_len * .4 )
        or $i < $last_col_80 - ( $aln_len * .4 );
    my $id_row_ref      = $$id_array[$i];
    my $gap_col_hashref = $gap_cols->[$i];
    my %gap_col_hash    = %{$gap_col_hashref};
    ## remove entire seq if %gaps in column is >= 80
    ## and (%id <= 50 or %nt_aligned <= 10%)
    if ( $gap_id_array[$i] >= 80.0
      and ( $id_row_ref->[1] <= 50 or $id_row_ref->[2] <= .1 ) )
    {
      foreach my $key ( keys %gap_col_hash ) {
        ## value of 1 means gap
        ## if there is no gap at this position where there are may gaps in the column
        ## we will make a note to remove this seq
        if ( $gap_col_hash{$key} != 1 ) {
          my $seq_obj = $aln_obj->get_seq_by_id($key);
          my @info = ( $key, $i + 1, $seq_obj );
          push @$ref2gspr, [@info];
          $$ref2gsr{$key} = $seq_obj;
        }
      }
    }
  }
## removes any sequences found in last step
  foreach my $key ( keys %$ref2gsr ) {
    my $seq_obj = $$ref2gsr{$key};
    $aln_obj->remove_seq($seq_obj);
  }
## removes any columns that are now all gaps
  $aln_obj = $aln_obj->remove_gaps( '-', 1 );

  my $test_aln_out = File::Spec->catpath( $volume, $out_path,
    $filename . ".removeLeast_trim_0" );
  my $out = Bio::AlignIO->new(
    -file             => ">$test_aln_out",
    -format           => 'fasta',
    -displayname_flat => 0
  );
  $out->write_aln($aln_obj);

  my ( $left_tir_start, $right_tir_start );
  ( $left_tir_start, $right_tir_start, $ref2gspr, $aln_obj, $ref2tp ) =
    consensus_filter( $ref2gspr, $aln_obj, 0, 0, $tir_positions );

  #@gap_seq_pos_remove = @{$ref2gspr};

  print
    "after cons_filter (filename.trim: $left_tir_start, $right_tir_start\n";

#my ($left_tir_start1,$right_tir_start1,$tmp_aln_obj,$ref2tp,$ref2gsr, $ref2gspr) = remove_most ($full_aln_obj,\%tir_positions, \@full_id_array);
  return ( $left_tir_start, $right_tir_start, $aln_obj, $tir_positions,
    $ref2gsr, $ref2gspr );

} ## end remove_least

sub get_org_aln {
my $infile = shift;
#create input object and an alignment object from it
my $in_obj = Bio::AlignIO->new( -file => $infile, -format => 'fasta' , -alphabet => 'dna');
my $aln_obj;
unless ( $aln_obj = $in_obj->next_aln() ) {
  die "Cannot find any alignment in file $infile\n";
}
#generate a gap column matrix from the full alignment
return ($aln_obj);
}
