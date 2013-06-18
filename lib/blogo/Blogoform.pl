#!/usr/bin/env perl

#print http head#########
print"Content-type:text/html\n\n";

require "setwebpath.pl";
if($sub_dir_under_webroot eq "/"){$sub_dir_under_webroot="";}
$seq_abs_dir=$ENV{'DOCUMENT_ROOT'}.$sub_dir_under_webroot."/blogo/pic/";
$page_HTTP_root=$sub_dir_under_webroot."/blogo/";

#default form data:
%formData=('jobname'=>"", 		'NTs'=>"checked",'AAs'=>"",'trNT'=>"",'codons'=>"",
	'plain_txt'=>"",'fasta'=>"checked",
	'sequences'=>"",
	'UPLOAD_FILE'=>"",		 	'auto'=>"checked",'bg_manual'=>"",
	'bg_fre'=>"0.25,0.25,0.25,0.25",
	'shortest'=>"checked",'re_manual'=>"",
	'reg_start'=>1,
	'reg_length'=>50,
	'first_num'=>1,
	'type1'=>"checked",'type2'=>"",
	'middle_total_height'=>200,
	'y_max'=>1,'y_min'=>-1,
	'char_width'=>25,
	'p_threshold'=>0.05
	);

if($ENV{"QUERY_STRING"}){
	$queryString = $ENV{"QUERY_STRING"};
	$queryString=~s/\+/ /g; 
	$queryString=~s/%(..)/pack("c",hex($1))/ge; 
	@tmpArray = split(/&/,$queryString); 
	foreach $curString(@tmpArray){ 
		($key,$value) = split(/=/,$curString); 
		$formData{$key}=$value; 
		if("seq_type seq_format logo_type"=~/$key/ || "auto shortest"=~/$value/){
			$formData{$value}="checked";
		}
		if($value eq "manual" && $key eq "bgfrqtype"){
			$formData{'bg_manual'}="checked";
		}
		if($value eq "manual" && $key eq "reg_type"){
			$formData{'re_manual'}="checked";
		}		
	}
	if($seq_file=$formData{'seq_file'}){ #sequences file
		open (SEQFILE,"$seq_abs_dir$seq_file");
		@seqs=<SEQFILE>;
		$formData{'sequences'}=join("",@seqs);
		close SEQFILE;
	}
}




print <<eof
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
<meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1" />
<title>Blogo-Run</title>
<link href="${page_HTTP_root}style.css" rel="stylesheet" type="text/css" />
</head>

<body>
<script language="javascript" src="${page_HTTP_root}form.js"></script>

<table width="97%" border="0" align="center" cellpadding="0" cellspacing="0">
  <tr>
    <td width="23%" bgcolor="#999999"><table width="115" border="0" cellspacing="0" cellpadding="0">
      <tr>
        <td width="26" bgcolor="#000000"><span class="B_style"><span class="STYLE4">B</span></span></td>
        <td width="89" bgcolor="#999999"><span class="STYLE9">L</span><span class="B_style"><span class="STYLE5">o</span></span><span class="STYLE9">g</span><span class="B_style"><span class="STYLE5">o</span></td>
      </tr>
    </table></td>
    <td width="77%" align="center" bgcolor="#999999"><span class="menu_style">&nbsp;&nbsp;&nbsp;<a href="${page_HTTP_root}index.htm" target="_self"> Home</a> | Run | <a href="${page_HTTP_root}example.htm" target="_self">Examples</a> | <a href="${page_HTTP_root}help.htm" target="_blank">Help</a> </span></td>
  </tr>
  <tr>
    <td colspan="2"><hr /></td>
  </tr>
  <tr>
    <td colspan="2">
	
	
    </td>
  </tr>
</table>
<table width="97%" border=0  cellspacing=0 align=center>

<form name=formlogo method="post" action="Blogo.pl" onSubmit="check_inputs();return false;" ENCTYPE="multipart/form-data">
<tr>
  <td width="146" valign="top"><a href="${page_HTTP_root}help.htm#job name" target="_blank">Job Name:</a></td>
  <td width="13">   
  <td width="581"><p>
  <input type=text name=jobname value=$formData{'jobname'}>
</p>
    <p>&nbsp; </p></tr><tr valign="top"></td>
<tr>
  <td valign="top"><span class="STYLE14"><a href="${page_HTTP_root}help.htm#Sequence Type" target="_blank">Your sequences:</a></span></td>
  <td>  
  <td>  <p>Type:
  <input type=radio name=seq_type onClick="treat_for_Nt()" value="NTs" $formData{'NTs'}>
  DNA or RNA &nbsp
  <input type=radio name=seq_type onClick="treat_for_AA()" value="AAs" $formData{'AAs'}>
  Protein &nbsp
  <input type=radio name=seq_type onClick="treat_for_AA()" value="trNT" $formData{'trNT'}>
  Translate to AA &nbsp
  <input type=radio name=seq_type onClick="treat_for_codon()" value="codons" $formData{'codons'}>
  Codon
  <br>
  Format: 
  <input type=radio name=seq_format value="fasta" $formData{'fasta'}>
  multi-fasta &nbsp
  <input type=radio name=seq_format value="plain_txt" $formData{'plain_txt'}>
  one sequence one line <br>
  <textarea rows=10 cols=80 name=sequences>$formData{'sequences'}</textarea>
    <br>or upload a file:<INPUT type="file" name="UPLOAD_FILE">
    
  </p>
    <p>        <br>
    </p></tr><tr valign="top"></td>


<tr>
  <td valign="top"><span class="STYLE14">  <a href="${page_HTTP_root}help.htm#Background frequencies" target="_blank">Background frequencies:</a><br>
</span></td>
  <td>  
  <td>	
	<p>
	  <input type=radio name=bgfrqtype value=auto $formData{'auto'}>
	  calculated from your sequences <br>
	    <input type=radio name=bgfrqtype value="manual" $formData{'bg_manual'}> 
	  input yourself:<br>
	  For DNA or RNA: (A,C,G,T)<br>
	  For proteins: (G,S,T,Y,C,Q,N,K,R,H,D,E,A,V,L,I,P,W,F,M)<br>
	  For codons: (A1 C1 G1 T1 A2 C2 G2 T2 A3 C3 G3 T3)<br>
	  <textarea rows=1 cols=80 name=bg_fre >$formData{'bg_fre'}</textarea>
	</p>
	<p>&nbsp;        </p></tr><tr valign="top"></td>


<tr><td valign="top"><span class="STYLE14"> <a href="${page_HTTP_root}help.htm#Choose region" target="_blank">Choose region:</a><br>
</span></td>
  <td>  
  <td>
	<p>
	  <input type=radio name=reg_type value="shortest" $formData{'shortest'}>
	  full length of your shortest sequence <br>
	    <input type=radio name=reg_type value="manual" $formData{'re_manual'}> 
	  input yourself:<br>
	  From
	  <input type=text onChange="set_manual_region()" name=reg_start size=5 value=$formData{'reg_start'}>
  &nbsp 
	  Length
  <input type=text onChange="set_manual_region()" size=5 name=reg_length value=$formData{'reg_length'}>
  <br>
  The number of the first symbol: <input type=text size=5 name=first_num value=$formData{'first_num'}>
	</p>
	<p>&nbsp;        </p></tr><tr valign="top"></td>


<tr><td valign="top"><span class="STYLE14"> <a href="${page_HTTP_root}help.htm#Logo properties" target="_blank">Logo properties:</a><br>
</span></td>
  <td>  
  <td>
	<p>Logo type 
	  <input type=radio name=logo_type onClick="treat_type1()" value=type1 $formData{'type1'}>
	  type 1 &nbsp
	  <input type=radio name=logo_type onClick="treat_type2()" value=type2 $formData{'type2'}>
	  type 2<br>
	  Middle height (Pixel) 
	  <input type=text name=middle_total_height size=5 value=200>
	  &nbsp<br>
	  Information Content: Y max 
	  <input type=text size=5 name=y_max value=$formData{'y_max'}>
	  &nbsp
	  Y min
	  <input type=text name=y_min size=5 value=$formData{'y_min'}>
	  <br>
	  Character width
	  <input type=text name=char_width size=5 value=25>
	  <br>

	  Image file format 
	  <select name=image_format>
	    <option selected value=jpg>jpg
	      <option value=gif>gif
	        <option value=tif>tif
          </select>
          <br>
          <input type=checkbox name=do_statis_test>
	  Statistical test: Only highlight symbols with P-Value less than: 
	  <input type=text name=p_threshold size=5 value=0.05>

    
	</p>
	<p>&nbsp;        </p></tr></td>

<tr><th colspan=3>
	<br><input type=submit value="submit"> 
	&nbsp <input type=reset value="reset" name="btn_reset">
	&nbsp <input type=button value="example" name="btn_example" onClick="treat_example()">
</th></tr>
</form>
</table>



</body>
</html>

eof

