#!/usr/bin/env perl

package Blogo;
require Exporter;
use Image::Magick;
use Statistics::ChisqIndep;
use Statistics::Fisher::twotailed;
use POSIX;
@ISA = qw(Exporter);
@EXPORT = qw( 
	run_Blogo
	); 

@EXPORT_OK = qw(bgfre); ##??

=global variables:
$seq_type,$reg_type,$reg_start,$reg_length,$bgfrqtype
%bgfre: background frequencies of elements
%posfre: elements frequencies in a position
$bg_fre_str
$middle_total_height, $y_max,$y_min
$image
$logo_type_method
$seq_number
$shortest_len
$img_name
$warning_text
$error_text
$newline_str
@fonts
=cut

#use strict;
@NTs=qw(A C G T); #nucleotides
@AAs=qw(G S T Y C Q N K R H D E A V L I P W F M); #amino acids
@codons=qw(A1 C1 G1 T1 A2 C2 G2 T2 A3 C3 G3 T3);
$font = 'Arial';
$margin_left=50;
$margin_top=20;
$margin_bottom=20;
$char_width=25;
%NTs_color=("A"=>"green","C"=>"blue","G"=>"orange2","T"=>"red");
%codons_color=("A"=>"green","C"=>"blue","G"=>"orange2","T"=>"red");
%AAs_color=(
	"G"=>"green","S"=>"green","T"=>"green","Y"=>"green","C"=>"green","Q"=>"green","N"=>"green",
	"K"=>"blue","R"=>"blue","H"=>"blue",
	"D"=>"red","E"=>"red",
	"A"=>"black","V"=>"black","L"=>"black","I"=>"black","P"=>"black","W"=>"black","F"=>"black","M"=>"black"
);
%newline_strs=("html"=>"<br>","text"=>"\n");

sub new {
  my $class = shift; # Get the request class name
  my $this = {};
  bless $this, $class; # Use class name to bless() reference
  return $this;
} 

sub treat_raw_seq_str{
	($FORM_pnt,$seq_str_pnt)=@_;
	$maybe_fasta=($$seq_str_pnt=~/>/);
	if($maybe_fasta && ($$FORM_pnt{'seq_format'} ne "fasta")){
		$warning_text.=("$newline_str Your sequences seems to be in fasta format (contain \">\"), 
			please check and choose the correct format next time. 
			The result was made by taking your input as in fasta format.$newline_str");
		$$FORM_pnt{'seq_format'}="fasta";
	}elsif(! $maybe_fasta && ($$FORM_pnt{'seq_format'} eq "fasta")){
		$warning_text.= "$newline_str Your sequences (do NOT contain \">\") seems to be NOT in fasta format, 
			please check and choose the correct format next time. 
			The result was made by taking your input as in plain text format.$newline_str ";
		$$FORM_pnt{'seq_format'}="plain_txt";
	}
	
	#treat the raw sequences:##########
	$$seq_str_pnt=~tr/uUa-z/TTA-Z/;
	if($$FORM_pnt{'seq_format'} eq "fasta"){
		$$seq_str_pnt=~s/\+//g;
		$$seq_str_pnt=~s/[^>]*>/>/;
		$$seq_str_pnt=~s/>.*\n/>/g;
		$$seq_str_pnt=~s/\s//g;
		#print"seq_str=$$seq_str_pnt\n";
		@sequences=split(/>/,$$seq_str_pnt);
		#@sequences=split(/>.*\n/,$$seq_str_pnt);
	}else{  # eq "plain_txt"
		$$seq_str_pnt=~s/\+//g;
		$$seq_str_pnt=~s/\+\t//g;
		$$seq_str_pnt=~s/[ \t]//g; 
		@sequences=split(/\s+/,$$seq_str_pnt);
	}
	#foreach (@sequences){print;print",$newline_str";}
	if(!$sequences[0]){shift(@sequences);}
}

sub pre_check_seq {
 #pre-check (1) whether the input are protein sequences or not; 
 # 	   (2)whether the sequences number is enouph for a type2 logo:
	($FORM_pnt,$sequences_pnt)=@_;
	$seq_number=@$sequences_pnt;
	my $is_right_input="";
	if($seq_type eq "AAs"){
		$is_right_input=guess_seq($$sequences_pnt[0]);
		if($is_right_input eq "DNA_RNA"){
			$error_text.="Your input seems to be DNA or RNA sequences
				 (NOT protein sequences), try again!";
			print"$newline_str$error_text$newline_str";
			die("Wrong input!");
		}
		if($$FORM_pnt{'logo_type'}eq"type2" && $seq_number<=100){
			$warning_text.="Sequences number $seq_number is NOT enouph for a meaningful Protein type2 logo, Input more sequences;";
		}
	}else{ # eq "NTs" "trNT" or "codons"
		$is_right_input=guess_seq($$sequences_pnt[0]);
		if($is_right_input eq "Protein"){
			$error_text.="Your input seems to be Protein sequences
				 (NOT DNA or RNAs), try again!";
			print"$newline_str$error_text$newline_str";
			die("Wrong input!");
		}
		if($$FORM_pnt{'logo_type'}eq"type2" && $seq_number<=20){
			$warning_text.="Sequences number $seq_number is NOT enouph for a meaningful DNA/RNA type2 logo, Input more sequences;";
		}
	}
}

sub translate_DNAs {
	my $sequences_pnt=@_[0];
	foreach $seq (@$sequences_pnt){
		$seq=translate_one_DNA($seq);
	}
}

sub translate_one_DNA {
	my $seq=@_[0];
	$seq=~s/\w{3}/$& /g;
	$seq=~s/ \w{1,2}$//g;
	$seq=~s/TTT|TTC/F/gi;
	$seq=~s/TAT|TAC/Y/gi;
	$seq=~s/CAT|CAC/H/gi;
	$seq=~s/CAA|CAG/Q/gi;
	$seq=~s/AAT|AAC/N/gi;
	$seq=~s/AAA|AAG/K/gi;
	$seq=~s/GAT|GAC/D/gi;
	$seq=~s/GAA|GAG/E/gi;
	$seq=~s/TGT|TGC/C/gi;
	$seq=~s/TTA|TTG|CTT|CTC|CTA|CTG/L/gi;
	$seq=~s/AGT|AGC|TCT|TCC|TCA|TCG/S/gi;
	$seq=~s/AGA|AGG|CGT|CGC|CGA|CGG/R/gi;
	$seq=~s/ATT|ATC|ATA/I/gi;
	$seq=~s/TGG/W/gi;
	$seq=~s/ATG/M/gi;
	$seq=~s/GTT|GTC|GTA|GTG/V/gi;
	$seq=~s/CCT|CCA|CCC|CCG/P/gi;
	$seq=~s/ACT|ACA|ACC|ACG/T/gi;
	$seq=~s/GCT|GCA|GCC|GCG/A/gi;
	$seq=~s/GGT|GGA|GGC|GGG/G/gi;
	$seq=~s/TAA|TAG|TGA//gi; #delete stop codon
	$seq=~s/ //g;
	return $seq;
}

sub guess_seq {
	#shift;
	my ($check_seq)=@_;
	#print"$check_seq \n";
	if($check_seq=~/[SYQNKRHDEVLIPWFM]/){return "Protein";}
	else{return "DNA_RNA";}
}

sub shortest_len_seqs {
	my $sequences_pnt=@_[0];
	#print "\$sequences_pnt=$sequences_pnt$newline_str";
	$shortest_len=0;
	my $seq_len,$i;
	for($i=0;$i<@$sequences_pnt;$i++){
		$$sequences_pnt[$i]=~s/\s//g;
		$seq_len=length ($$sequences_pnt[$i]);
		if($seq_len<$shortest_len||$shortest_len==0){$shortest_len=$seq_len;}
	}
	return $shortest_len;
}

sub region_for_draw{
	#my $reg_length, $reg_start;
	my ($FORM_pnt,$shortest_len,$sequences_pnt)=@_;
	if($reg_type eq "shortest"){
		$reg_length=$shortest_len;
		$reg_start=1;
	}elsif($reg_type eq"manual"){ 
		if($reg_start+$reg_length-1>$shortest_len){
			$warning_text.="The start position $reg_start plus sequences length $reg_length for drawing a logo is 
			bigger than $shortest_len (the shortest length of your inputed sequences). ";
			$reg_length=$shortest_len-$reg_start+1;
			$warning_text.="We have set the length to draw a logo = $reg_length (
			the shortest length of your inputed sequences minus The start position);$newline_str$newline_str";
		}
	}
	#print $$FORM_pnt{'reg_start'}.", ".$$FORM_pnt{'reg_length'}.$newline_str;
	#return($reg_start,$reg_length);
}

sub cal_bg_fre{ #calculate background nt frequency of DNA or RNA
	my ($FORM_pnt,$sequences_pnt)=@_;
	#$seq_type=$$FORM_pnt{'seq_type'};#NTs, AAs or codons
	my $i;
	if($bgfrqtype eq"manual"){
		my @bg_fre_input=split(/[^0-9.]+/,$$FORM_pnt{'bg_fre'});
		if(!$bg_fre_input[0]){shift @bg_fre_input;}
		%bgfre={};
		my $fre_summary=0;
		for($i=0;$i<@$seq_type;$i++){
			$bgfre{$$seq_type[$i]}=$bg_fre_input[$i];
			$fre_summary+=$bg_fre_input[$i];
		}
		if($seq_type eq "codons" && abs($fre_summary-3)>0.005){
			$bgfrqtype="auto";
			$error_text.="The summary of your inputed frequencies 
				is not 3 (3 is for all the three positions of codon, the sum of each position is 1). 
				Do not use % in the frequencies. The program
				 has calculated the background frequencies automatically; $newline_str";
		}elsif($seq_type ne "codons" && abs($fre_summary-1)>0.005){
			$bgfrqtype="auto";
			$error_text.="The summary of your inputed frequencies 
				is not 1. Do not use % in the frequencies. The program
				 has calculated the background frequencies automatically; $newline_str";
		}else{
			$bg_fre_str=join(', ',@bg_fre_input);
		}
	}
	if($bgfrqtype eq "auto"){ # eq"auto"
		my $temp_seq_str, $temp_seqlen;
		if($seq_type eq "codons"){
			$temp_seq_str="";
			$temp_seqlen=$reg_length-$reg_length%3; # codon should be 3x
			foreach (@$sequences_pnt){
				if(!$_){next;}
				$temp_seq_str.=substr($_,$reg_start-1,$temp_seqlen);
			}
			$temp_seq_str=~s/\w{3}/$& /g;
			$temp_seq_str=~s/\b\w/$&1/g;
			$temp_seq_str=~s/\B[a-zA-Z]\B/$&2/g;
			$temp_seq_str=~s/\w\b/$&3/g;
			$temp_seq_str=~s/ //g;
		}else{ # $seq_type eq NTs or AAs
			$temp_seq_str=join('',@$sequences_pnt);
		}
		my $total_len=length ($temp_seq_str);
		my $total_len_temp=$total_len,$total_len_true=0;
		my %total_num={},$current_num;
		%bgfre={};
		$bg_fre_str="";
		foreach (@$seq_type){
			#print"$temp_seq_str $newline_str ";
			$temp_seq_str=~s/$_//g;
			$current_num=length($temp_seq_str);
			#print"$current_num ";
			$total_num{$_}=$total_len_temp-$current_num;
			$total_num_true+=$total_num{$_};
			$bgfre{$_}=$total_num{$_}/$total_len;
			if($seq_type eq "codons"){$bgfre{$_}*=3;}
			$total_len_temp=$current_num;
			$bg_fre_str.=sprintf ("%6.4f, ", $bgfre{$_});
		}
		#print"$bg_fre_str$newline_str ";
		if($temp_seq_str){$error_text.="There are unexpected charactor(s) 
			 in your input sequences; $newline_str";
			$bg_fre_str="";
			foreach (@$seq_type){
				$bgfre{$_}=$total_num{$_}/$total_num_true;
				if($seq_type eq "codons"){$bgfre{$_}*=3;}
				$bg_fre_str.=$bgfre{$_}.", ";
			}
			$total_len=$total_num_true;
		}
	}
	$bg_fre_str="($seq_type) @$seq_type, ".$bg_fre_str;
}


sub draw_logo{
	my ($FORM_pnt,$sequences_pnt,$img_abs_dir)=@_;
	$logo_type_method="draw_logo_".$$FORM_pnt{'logo_type'};
	#print"\$logo_type_method=$logo_type_method$newline_str";
	draw_frame($FORM_pnt);
	my $pos;
	my $chi = new Statistics::ChisqIndep;
	for($pos=$reg_start-1;$pos<$reg_start+$reg_length-1;$pos++){
		my %posfre={},%test_p_value={};
		my %IC={},$IC_total=0,$IC_positive_total=0,$element,$element_bg;
		foreach (@$sequences_pnt){
			$element=substr($_,$pos,1);
			$posfre{$element}++;
		}
		foreach $element(keys %posfre) {
			$posfre{$element}/=$seq_number;
			if($posfre{$element}==0){next;}
			$element_bg=$element;
			if($seq_type eq "codons"){ #add codon position
				$element_bg.=($pos+1-$reg_start)%3+1;
				#print"$element_bg $newline_str ";
			}
			
			if($$FORM_pnt{'do_statis_test'}){
				$observed_num=int($seq_number*$posfre{$element}+0.5);
				$expected_num=int($seq_number*$bgfre{$element_bg}+0.5);
				if($seq_number>400 && $observed_num>20 && $expected_num>20){ #do chi square test
					my $chisq_test_array=[[$observed_num,$seq_number-$observed_num],[$expected_num,$seq_number-$expected_num]];
					$chi->load_data($chisq_test_array);
					$test_p_value{$element}=$chi->{p_value};
					#print "\$pos=$pos, $element p-value: ", $chi->{p_value}, "<br>";
				}else{# do fisher's exact test
					$test_p_value{$element} = calculateStatistic( n11=>$observed_num,
                                      n1p=>$observed_num+$expected_num,
                                      np1=>$seq_number,
                                      npp=>2*$seq_number);
				}
			}
			#print"$element: $posfre{$element}; $element_bg: $bgfre{$element_bg}"; 
			$IC{$element}=$posfre{$element}*(log($posfre{$element}/$bgfre{$element_bg})/log(2));
			#print" IC: $IC{$element}$newline_str ";
			$IC_total+=$IC{$element};
			if($IC{$element}>0){$IC_positive_total+=$IC{$element};}
		}
		#print"$seq_number $posfre{'A'}$newline_str";
		&$logo_type_method($pos,\%posfre,\%IC,$IC_total,$IC_positive_total,\%test_p_value,$FORM_pnt);
	}
	$image->Write($img_abs_dir.$img_name);
}

sub draw_logo_type1{
	my ($pos,$posfre_pnt,$IC_pnt,$IC_total,$IC_positive_total,$test_p_value_pnt,$FORM_pnt)=@_;
	my $logo_total_height=$middle_total_height*$IC_total/($y_max-$y_min);
	#print$logo_total_height.$newline_str;
	my $base_height=$middle_total_height-$logo_total_height;
	my $ele_num=keys %$posfre_pnt; #elements number
	#print$ele_num.$newline_str;
	my $i,$max=0,$element,$max_elem,$elem_height,$elem_height_scale,$elem_x,$elem_y,$color,$color_symbol;
	for($i=1;$i<$ele_num;$i++){
		$max=0;
		foreach $element (keys %$posfre_pnt){
			if($$posfre_pnt{$element}>$max){
				$max=$$posfre_pnt{$element};
				$max_elem=$element;
			}
		}
		$elem_height=$logo_total_height*$max;
		$elem_height_scale=$elem_height/$char_width*1.4;
		$elem_x=$margin_left+($pos-$reg_start+1)*$char_width;
		$base_height+=$elem_height;
		$elem_y=$margin_top+$base_height;
		delete($$posfre_pnt{$max_elem});
		if($$FORM_pnt{'do_statis_test'} && $$test_p_value_pnt{$max_elem}>$$FORM_pnt{'p_threshold'}){
			$color_symbol="Gray";
		}else{
			$color=$seq_type."_color";
			$color_symbol=$$color{$max_elem};
		}
		
		$image->Annotate(font=>$fonts[@fonts-4], pointsize=>$char_width,text=>$max_elem,fill=>$color_symbol,
			x=>$elem_x,y=>$elem_y,scale=>"1.4,".$elem_height_scale);
	}
}

sub draw_logo_type2{
	my ($pos,$posfre_pnt,$IC_pnt,$IC_total,$IC_positive_total,$test_p_value_pnt,$FORM_pnt)=@_;
	my $logo_positive_height=$middle_total_height*$IC_positive_total/($y_max-$y_min);
	#print$logo_total_height.$newline_str;
	my $base_height=$middle_total_height*$y_max/($y_max-$y_min)-$logo_positive_height;
	my $ele_num=keys %$IC_pnt; #elements number
	#print$ele_num.$newline_str;
	my $i,$max,$element,$max_elem,$elem_height,$elem_height_scale,$elem_x,$elem_y,$color,$color_symbol;
	for($i=1;$i<=$ele_num;$i++){
		$max=-50;
		foreach $element (keys %$IC_pnt){
			if($$IC_pnt{$element}>$max){
				$max=$$IC_pnt{$element};
				$max_elem=$element;
			}
		}
		#if($$IC_pnt{$element}==0){next;}
		#& $$IC_pnt{$element}!=0
		$elem_height=$middle_total_height*abs($max)/($y_max-$y_min);
		$elem_height_scale=$elem_height/$char_width*1.4;
		$elem_x=$margin_left+($pos-$reg_start+1)*$char_width;
		$base_height+=$elem_height;
		$elem_y=$margin_top+$base_height;
		delete($$IC_pnt{$max_elem});
		if($$test_p_value_pnt{$max_elem}>$$FORM_pnt{'p_threshold'}){
			$color_symbol="grey";
		}else{
			$color=$seq_type."_color";
			$color_symbol=$$color{$max_elem};
		}
		$image->Annotate(font=>$fonts[@fonts-4], pointsize=>$char_width,text=>$max_elem,fill=>$color_symbol,
			x=>$elem_x,y=>$elem_y,scale=>"1.4,".$elem_height_scale);
	}
}

sub draw_frame{
	my ($FORM_pnt)=@_;
	$char_width=$$FORM_pnt{'char_width'};
	$middle_total_height=$$FORM_pnt{'middle_total_height'};
	my $img_width=$margin_left+$reg_length*$char_width+1;
	my $img_height=$margin_top+$middle_total_height+$margin_bottom;
	my $img_size="${img_width}x${img_height}";
	my $temp,$i,$y;
	$image = Image::Magick->new(size=>$img_size,pointsize=>$char_width);
	$image->Read("xc:white");
	
  	@fonts = $image->QueryFont();
  	#foreach (@fonts){print;print"<br>";}
  	
	$temp=($margin_left-1).",".($margin_top-1).",".($img_width-1).",".($img_height-$margin_bottom-1);
	$image->Draw(primitive=>'rectangle', points=>$temp, fill=>"black");
	$temp=($margin_left).",".($margin_top).",".($img_width-2).",".($img_height-$margin_bottom-2);
	$image->Draw(primitive=>'rectangle', points=>$temp, fill=>'white');
	
	for($i=0;$i<=4;$i++){
		$y=$margin_top+$middle_total_height/4*$i;
		$temp=($margin_left-4).",".($y-1).",".($margin_left-1).",".($y-1);
		$image->Draw(primitive=>'line',points=>"$temp");
		$y_max=$$FORM_pnt{'y_max'};
		$y_min=$$FORM_pnt{'y_min'};
		$image->Annotate(font=>$fonts[@fonts-4], pointsize=>18,text=>$y_max-($y_max-$y_min)/4*$i,x=>$margin_left-3,y=>$y+6,align=>'right');
	}
	for($i=0;$i<$reg_length;$i+=3){
		$image->Annotate(font=>$fonts[@fonts-4], pointsize=>18,text=>$i+($$FORM_pnt{'first_num'}),x=>$margin_left+$i*$char_width+5,y=>$img_height-$margin_bottom+20);
	}
}

sub print_output{
	# information output:
	($FORM_pnt,$sequences_pnt,$imghttp_path)=@_;
	if($$FORM_pnt{'output_type'} eq "html"){
		#print http head#########
		print"Content-type:text/html\n\n";
		print"<meta http-equiv=\"Content-Type\" content=\"text/html; charset=\" />";
		print"<title>Blogo Result</title>";
		print"<body><h2>Blogo Result:</h2>";
	}
	
	if($$FORM_pnt{'output_type'} eq "remote_pipe"){
		print"Content-type:text/html\n\nBlogo:\n\n";
		#my $xml = XML::Simple->new();
		%out_hash=(
			'jobname'=>$$FORM_pnt{'jobname'},
			'Blogoimg'=>${imghttp_path}.$img_name,
			'logo_type'=>$$FORM_pnt{'logo_type'},
			'bg_fre_str'=>$bg_fre_str
		);
		foreach $out_name (keys %out_hash){
			print "$out_name=$out_hash{$out_name}\n";
		}
		#my $xml_content = $xml->XMLout($out_hash,NoAttr => 1);
		#print (Dumper $out_hash);
		#print $xml_content;
	}else{
		print($newline_str.("="x50).$newline_str);
		print"$newline_str Job name: $$FORM_pnt{'jobname'}$newline_str";
		print("<img src=${imghttp_path}$img_name>$newline_str");
		print("download Blogo image: <a target=_blank href=${imghttp_path}$img_name> $img_name</a>$newline_str");
		if($warning_text ne "$newline_str<b>Warning: </b>$newline_str"){print$warning_text;}
		if($error_text ne "$newline_str<b>Error: </b>$newline_str"){print$error_text;}
		print"$newline_str$newline_str Your input sequences are:$newline_str";
		my $num=0;
		foreach (@$sequences_pnt){
			if(!$_){next;}\
			$num++;
			if($num>200){
				print"......(more)......$newline_str";
				last;
			}
			print substr($_,$reg_start-1,$reg_length)." $newline_str";
		}
		print "$newline_str Start from $reg_start, Analysis length=$reg_length $newline_str";
		print "$newline_str Logo type = $$FORM_pnt{'logo_type'}$newline_str";
		print("$newline_str Background frequency:$newline_str$bg_fre_str$newline_str");
	}
	if($$FORM_pnt{'output_type'} eq "html"){
		print"</BODY>";
	}
}

sub run_Blogo{
	shift;
	($FORM_pnt,$seq_str_pnt,$imghttp_path)=@_;
	$newline_str=$newline_strs{$$FORM_pnt{'output_type'}};
	$warning_text="$newline_str<b>Warning: </b>$newline_str";
	$error_text="$newline_str<b>Error: </b>$newline_str";
	treat_raw_seq_str($FORM_pnt,$seq_str_pnt);
	$sequences_pnt=\@sequences;
	($seq_type,$reg_type,$reg_start,$reg_length,$bgfrqtype)=@$FORM_pnt{'seq_type','reg_type','reg_start','reg_length','bgfrqtype'};
	if($seq_type eq"trNT"){
		pre_check_seq($FORM_pnt,$sequences_pnt);
		if($reg_type eq"shortest"){
			shortest_len_seqs($sequences_pnt);
			region_for_draw($FORM_pnt,$shortest_len,$sequences_pnt);
		}
		foreach (@$sequences_pnt){
			if(!$_){next;}
			$_=substr($_,$reg_start-1,$reg_length);
		}
		translate_DNAs($sequences_pnt);
		$seq_type="AAs";
		$reg_type="manual";
		$reg_length=($reg_length-$reg_length%3)/3;
		$reg_start=1;
	}
	#pre-check (1) whether the input are protein sequences or not; 
	# 	   (2)whether the sequences number is enouph for a type2 logo:
	pre_check_seq($FORM_pnt,$sequences_pnt);
	shortest_len_seqs($sequences_pnt);
	# the region for drawing logos.##########
	region_for_draw($FORM_pnt,$shortest_len,$sequences_pnt);
	# background frequency of nt or amino acids:#############
	cal_bg_fre($FORM_pnt,$sequences_pnt);
	#calculate frequency on every position and draw a logo
	$img_name=$$FORM_pnt{'jobname'}."_".$seq_type."_B_".$$FORM_pnt{'logo_type'}.".".time().".".$$FORM_pnt{'image_format'};
	$img_abs_dir=$$FORM_pnt{'img_abs_dir'};
	#if(-e $img_abs_dir.$img_name){unlink $img_abs_dir.$img_name;}
	draw_logo($FORM_pnt,$sequences_pnt,$img_abs_dir);
	print_output($FORM_pnt,$sequences_pnt,$imghttp_path);
}

=LICENSE
<Blogo, display bias of biological sequences> Copyright (C) <2008> <Wencheng Li, lwcbio@yahoo.com.cn> 

This library is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License (http://www.gnu.org/licenses/lgpl.html) for more details.

Note: a copy of the GNU Lesser General Public License is available on the web
at <http://www.gnu.org/licenses/lgpl.html> and is included in this
distribution as GNU_Lesser_General_Public_License.txt.

=cut


1;
