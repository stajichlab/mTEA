#!/usr/bin/env perl

use LWP::Simple;
use Data::Dumper;
use HTTP::Request::Common;
require LWP::UserAgent;

%query_str=(
	'jobname'=>"temp25",
	'seq_type'=>"AAs",
	'seq_format'=>"plain_txt",
	'sequences'=>"",
	'UPLOAD_FILE'=>["E:/lwc/website/Blogo/seq/test1.seq"],
	'bgfrqtype'=>"auto",
	'bg_fre'=>"0.25,0.25,0.25,0.25",
	'reg_type'=>"shortest",
	'reg_start'=>1,
	'reg_length'=>50,
	'first_num'=>1,
	'logo_type'=>"type2",
	'middle_total_height'=>200,
	'y_max'=>1,
	'y_min'=>-1,
	'char_width'=>25,
	'image_format'=>"jpg",
	'do_statis_test'=>"on",
	'p_threshold'=>0.05,
	'output_type'=>"remote_pipe");

#UPLOAD_FILE should be changed
#output_type=remote_pipe(default), html or text

# Set the path of Blogo server, use any one of them:
$website="http://acephpx.cropdb.org/cgi-bin/Blogo/Blogo.pl";
#$website="http://www.bioinformatics.org/blogo/cgi-bin/Blogo/Blogo.pl";


my $ua = LWP::UserAgent->new;
$ua->timeout(25);
#$ua->env_proxy; # set proxy:
#$ua->proxy(['http'], 'http://123.123.123.123:80/');
#$ua->proxy(['ftp'], 'http://123.123.123.123:21/');

my $response=$ua->request(POST $website,'Content_Type'=>"multipart/form-data", 'Content'=>\%query_str);

if($response->is_success){
	$content=$response->content;
}else{
	die $response->status_line;
}
print	$content;
