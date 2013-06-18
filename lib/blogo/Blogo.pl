#!/usr/bin/env perl

#read form papameters:
use CGI;
$req = new CGI;
%FORM = %{$req->Vars};

#foreach $key (keys %ENV){
#	print"$key=$ENV{$key}<br>";
#}
#set path:##########
require "setwebpath.pl";
if($sub_dir_under_webroot eq "/"){$sub_dir_under_webroot="";}
$img_abs_dir=$ENV{'DOCUMENT_ROOT'}."$sub_dir_under_webroot/blogo/img/";
$FORM{'img_abs_dir'}=$img_abs_dir;
$imghttp_path="http://".$ENV{'SERVER_NAME'}."$sub_dir_under_webroot/blogo/img/";
#$imghttp_path="$sub_dir_under_webroot/blogo/img/";
$cgi_abs_dir=$ENV{'SCRIPT_FILENAME'};
$cgi_abs_dir=~s|/[^/]*$|/|;
unshift (@INC , $cgi_abs_dir);
use Blogo;
$logo_obj=Blogo::new();
if(!$FORM{'output_type'}){$FORM{'output_type'}="html";}

my $file = $req->param("UPLOAD_FILE");
if($file){
	read_upload_file($file);
}else{
	$sequences=$FORM{'sequences'};
}

#run Blogo:##############
$logo_obj->run_Blogo(\%FORM,\$sequences,$imghttp_path);

#delete image file more than XXX hours
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = gmtime(time);
if($hour>=0 && $hour<=2 && $wday>4){
	#delete_old_temp_files($img_abs_dir,8);
}
delete_old_temp_files($img_abs_dir,2);

sub delete_old_temp_files{
	($abs_dir,$max_hours)=@_;
	opendir (DIR, $abs_dir);
	my $filename;
	while ($filename=readdir DIR){
		if (-f $abs_dir.$filename){
			my ($dev,$ino,$mode,$nlink,$uid,$gid,$rdev,$size,$filenametime,$mtime,$ctime,
			$blksize,$blocks)=stat $abs_dir.$filename;
			if((time()-$ctime)/3600>$max_hours){
				unlink $abs_dir.$filename;
			}
		}
	}
	close DIR;
}

sub read_upload_file{
	my ($file)=@_;
	my $fileName = $file;
	$fileName =~ s/^.*(\\|\/)//; 
	open (UPLOADFILE,">$img_abs_dir$fileName");
	binmode(UPLOADFILE);
	while (my $bytesread = read($file, my $buffer, 1024)) {
		print UPLOADFILE $buffer;
	}
	close (UPLOADFILE);

	open(UPLOADFILE,"$img_abs_dir$fileName");
	my $line;
	while($line=<UPLOADFILE>){
		$sequences.=$line;
	}
	close (UPLOADFILE);
}

=LICENSE
<Blogo, display bias of biological sequences> Copyright (C) <2008> <Wencheng Li, lwcbio@yahoo.com.cn> 

This library is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License (http://www.gnu.org/licenses/lgpl.html) for more details.

Note: a copy of the GNU Lesser General Public License is available on the web
at <http://www.gnu.org/licenses/lgpl.html> and is included in this
distribution as GNU_Lesser_General_Public_License.txt.


=cut
