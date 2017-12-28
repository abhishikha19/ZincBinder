#!/usr/bin/perl 
############################### Header Information ##############################               
require 'cgi.perl';
use CGI;;
$query = new CGI;
&ReadParse;
print &PrintHeader;
################################ Reads Input Data ##############################                
$seq = $query->param('seq');
$file = $query->param('file');
$svm_th = $query->param('svm_fpr');
#################Validation Of Input Sequence Data (file upload) ###################################        
if($file ne '' && $seq eq '')
{
    $file=~m/^.*(\\|\/)(.*)/;
    while(<$file>)
    {
        $seqfi .= $_;
    }
}
elsif($seq ne '' && $file eq '')
{

   $seqfi="$seq";
}
##############ACTUAL PROCESS BEGINS FROM HERE#######################                            
$ran= int(rand 10000);
$dir = "/webservers/cgidocs/mkumar/temp/shikha/znbinder/zinc_$ran";
system "mkdir $dir";
system "chmod 777 $dir";
$nam = "seq".".input";
open(FP1,">$dir/$nam");
if($seqfi !~ m/\>/)
{
    print FP1 ">seq\n";
}
print FP1 "$seqfi\n";
close FP1;
open(FH_FINAL,">$dir/seq.final") or die "$!";
open(FH_SEQ,"$dir/seq.input") or die "$!";
while($line_seq=<FH_SEQ>)
{
    chomp($line_seq);
    $line_seq=~ s/[^_>a-zA-Z0-9]//g;
    print FH_FINAL "$line_seq\n";
}
close FH_SEQ;
close FH_FINAL;
open(FP9,">>/webservers/cgidocs/mkumar/temp/shikha/znbinder/zinc.que");                       
print FP9 "$dir:$svm_th";
close(FP9);
if(!defined($pid = fork)){
    die;
}
elsif(!$pid){
    close(STDIN);close(STDOUT);close(STDERR);
    sleep(30);
    system "/usr/bin/perl /webservers/cgi-bin/znbinder/zn_pred2.pl $dir";
    exit;
}
else{
    print  "</head><body>\n";
    print  "<table ALIGN = \"CENTER\" border=\"0\" cellpadding=\"0\" cellspacing=\"0\" class=\"tbl1\" width=\"100%\"><tr><td colspan=\"4\"></td></tr>\n";
    print  "<td width=\"92%\"><br><h2 ALIGN = \"CENTER\"> Sequence Submitted for Prediction</h2><HR ALIGN =\"CENTER\"> </HR>\n";
    print  "<font size=3 color=black ><b><center>Thanks for using zincbinder Web-server</center></b></font>";
    print  "<p>";
    print  "<font size=3 color=black><b><center>Your request is being processed. Please wait for a few minutes till the process completes. ";
    print  "If you have any problem or suggestions please contact Dr. Manish Kumar <a href='mailto:manish\@south.du.ac.in'><font size=4 color red><b>[manish\@south.du.ac.in\]</b></font></a>.<p>Your job number is <font color=red>$ran</font>.<br></p></b>";
    print  "<p>";
    print  "<meta http-equiv=\"refresh\" content=\"20;url=/cgi-bin/znbinder/chkres?c=$ran\"></center>\n";
    print  "</table><h2>&nbsp;</h2></td></tr></table>\n";
    print  "</body>";
}



