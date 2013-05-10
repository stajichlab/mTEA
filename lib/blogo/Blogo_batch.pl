#!/usr/bin/env perl

=global variables:
%FORM
@seqfile_names
@sequences
=cut

sub read_config {
	# read config file:
	my $config_file;
	if(!$FORM{"conf_file"}){$FORM{"conf_file"}="Blogo.conf";}
	open (CONFIG_FILE,$FORM{"conf_file"}) || die "Can not find Blogo.conf file! \n";
	print @result;
	while (<CONFIG_FILE>){
		s/#.*//g;
		s/\s//g;
		chomp;
		if($_){
			@row=split(/=/,$_);
			print "$_\n";
			if(!$FORM{$row[0]}){$FORM{$row[0]}=$row[1];}
		}
	}
	$FORM{'output_type'}="text";
	@seqfile_names=split(/,/,$FORM{'file_names'});
	#print"\n@seqfile_names\n";
	close (CONFIG_FILE);
}

sub read_seq_file {
	use Blogo;
	$logo_obj=Blogo::new();
	#read sequense file:
	$seq_format=$FORM{'seq_format'};
	foreach $seqfile_name (@seqfile_names){
		$FORM{'jobname'}=$seqfile_name;
		my $input_file=$FORM{'file_path'}.$seqfile_name;
		open (INPUT_SEQ_FILE,$input_file) || die "Can not find $input_file! \n";
		@seq_content_array=<INPUT_SEQ_FILE>;
		#print@seq_content_array;
		my $sequences=join ("",@seq_content_array);
		$logo_obj->run_Blogo(\%FORM,\$sequences,$FORM{'img_abs_dir'});
		close (INPUT_SEQ_FILE);
	}
}

BEGIN {
	print(("="x50)."\n Hi! Welcome to Blogo!\n");
	if(@ARGV){
		my $blockstr,@block_value;
		foreach $blockstr (@ARGV){
			@block_value=split(/=/,$blockstr);
			$FORM{$block_value[0]}=$block_value[1];
		}
	}
	read_config;
	read_seq_file;
}


=LICENSE
<Blogo, display bias of biological sequences> Copyright (C) <2008> <Wencheng Li, lwcbio@yahoo.com.cn> 

This library is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License (http://www.gnu.org/licenses/lgpl.html) for more details.

Note: a copy of the GNU Lesser General Public License is available on the web
at <http://www.gnu.org/licenses/lgpl.html> and is included in this
distribution as GNU_Lesser_General_Public_License.txt.

=cut
