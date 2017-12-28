#!/usr/bin/perl 
###############################Read input data from zn_pred1.pl ##############################
$random=$ARGV[0];
print"$random\n";
chomp($random);
open(FT,"/webservers/cgidocs/mkumar/temp/shikha/znbinder/zinc.que");
while($liner=<FT>)
{
    chomp($liner);
    @cod=split(/:/,$liner);
    $fnum="$cod[0]";
    @directory = split('/',$cod[0]);
    $usr1="usr1";
    $fnum=~s/$usr1//g;
    $fnum=~s/[^0-9]//g;
    if($fnum==$random)
    {
        $dir = "$cod[0]";
        $thr = "$cod[1]";
    }
}
$num_of_seq= `grep -c '>' $random/seq.input`;
system "/usr/bin/perl /webservers/cgi-bin/znbinder/splitfasta.pl $random/seq.final";
system "mv /webservers/cgi-bin/znbinder/*.fa $random";
###############################actual process begins from here ##############################
if($num_of_seq <=5)
{
    system "/usr/bin/perl /webservers/cgi-bin/znbinder/fasta.pl $random/seq.final > $random/seq_2line";
}
if ($num_of_seq >5)
{
    system "/usr/bin/perl /webservers/cgi-bin/znbinder/fasta.pl  $random/seq.final > $random/seq_2line1";
    system "head -10 $random/seq_2line1 > $random/seq_2line";
}
system "/usr/bin/perl  /webservers/cgi-bin/znbinder/CHED_extract_1.pl   $random/seq_2line  > $random/CHED_result1.txt";
system " sort -n $random/CHED_result1.txt >$random/CHED_result1_sort.txt";
system "/usr/bin/perl  /webservers/cgi-bin/znbinder/CHED_extract_2.pl   $random/seq_2line  >$random/CHED_result2.txt";
system "/usr/bin/perl  /webservers/cgi-bin/znbinder/patgen.pl -i $random/CHED_result2.txt  -w  19  >$random/CHED_win19.pat";
system "mkdir $random/PSSM";
system "chmod 777 $random/PSSM";
$a=0;
open(FILE, "$random/seq_2line")or die;
while ($line1=<FILE>)
{
    chomp($line1);
    if ($line1=~m/^>/)
    {
        $nextline=<FILE>;
        chomp ($nextline);
	$app = "$line1"."\n"."$nextline";
	$pssm[$a]=$app;
	$line1=~ s/^>//;
	$pssm_out="$line1".'_pssm';
	$blast_out="$line1".'.out';
	open(TMP,">$random/temp") or die "$!";
	print TMP "$pssm[$a]\n";
	close TMP;
	system "cat $random/temp";
	$a++;
	system "/usr/local/bin/blast-2.2.18/bin/blastpgp  -d /webservers/cgi-bin/znbinder/nr90/nr90 -i $random/temp -j 3 -h 0.001 -m 8 -Q $random/$pssm_out -o $random/$blast_out";
	system "mv $random/$pssm_out $random/PSSM/.";
    }
}
system "/usr/bin/perl /webservers/cgi-bin/znbinder/extract_pssm.pl $random/CHED_win19.pat  $random/seq_2line $random/pssm_temp $random/PSSM  > $random/CHED_win19.mtx";
system "/usr/bin/cut -d '#' -f2 $random/CHED_win19.mtx > $random/CHED_final_win19.mtx";
system "/usr/local/bin/svm_classify   $random/CHED_final_win19.mtx  /webservers/cgi-bin/znbinder/level1_model  $random/svm_score";
system "/usr/bin/cut -d ':' -f1 $random/CHED_win19.pat > $random/CHED_id";
system "/usr/bin/cut -d ':' -f3 $random/CHED_win19.pat > $random/CHED_pattern";
system "/usr/bin/paste -d ':' $random/CHED_id  $random/CHED_pattern  >$random/CHED_id_pattern";
system "/usr/bin/paste $random/CHED_id_pattern  $random/CHED_result1_sort.txt  $random/svm_score > $random/total_result";
###############################Final result processing################################################
system "cat 1.html >>$random/result.html";
open(FH2,">>$random/result.html") or die "$!";
print FH2 "<h4 align ='center'><font color ='#990033'>Zincbinder Prediction Result</font></h4>\n";
print FH2 "<table align='center' border = '0'>";
print FH2 "<tr hight='100'><td height='20' colspan='4' align='center'><font color ='#990033'><h5>Query Search Detail</h5></font></td></tr>";
print FH2 "<tr bgcolor='#F3E5AB'><td><font color ='#990033'>JOB-ID</font></td><td align ='center' colspan='3' >$directory[7]</td></tr>";
print FH2 "<tr bgcolor='#F3E5AB'><td><font color ='#990033'>Number of Query Sequences</font></td><td align ='center' colspan='3'>$num_of_seq</td></tr>";
print FH2 "<tr bgcolor='#F3E5AB'><td><font color ='#990033'>Predicted on</font></td><td align ='center' colspan='3'><iframe src='http://free.timeanddate.com/clock/i57l4a99/n176/tcf3e5ab' frameborder='0' width='114' height='30'></iframe></td></tr>";
###################################print HTML result#######################################
print FH2 "<tr hight='100'><td height='20' colspan='4' align='center'><font color ='#990033'><h5>Prediction Result</h5></font></td></tr>";
print FH2 "<tr>";
print FH2 "<th bgcolor=\"\#5AD7D7\">Protein-ID</th>";
print FH2"<th  bgcolor=\"\#5AD7D7\">zinc binding residue</th>";
print FH2 "<th bgcolor=\"\#5AD7D7\">position</th>";
print FH2 "<th bgcolor=\"\#5AD7D7\">znbinding score</th>";
print FH2 "</tr>";
open (FH, "$random/total_result") or die "$!";
$count=0;
while ($line=<FH>)
{
    chomp($line);
    @value=split(/\t/,$line);
    @aa=split(/:/,$value[0]);
    $aa[0]=~ s/\>//g;
    @CHED=split(//,$aa[1]);
    $vv = sprintf("%.1f",$value[2]);
    if ($vv >= $svm_th)
    {
	$count++;
	if($count%2==0)
      {
	  print FH2 "<tr hight='100' bgcolor='#CCFFFF'>";
	  print FH2 "<td align='center'>$aa[0]</td><td align='center'>$CHED[9]</td><td align='center'>$value[1]</td><td align='center'>$vv</td>\n";
	  print FH2 "</tr>";
      }
    elsif($count%2==1)
    {
        print FH2 "<tr hight='40' bgcolor='#FFFFCC'>";
	print FH2 "<td align='center'>$aa[0]</td><td align='center'>$CHED[9]</td><td align='center'>$value[1]</td><td align='center'>$vv</td>\n";
	print FH2 "</tr>";
    }
 	
}
    
 }
close FH;
print FH2 "</table></p>";
print FH2"</div>";
print FH2 "<div id='page-bottom'><div id='page-bottom-content'><p>ZincBinder is maintained by Department of Biophysics, University of Delhi South Campus. We acknowledge your comments orcontribution to this resources. Please send your sugestions to Dr. Manish Kumar (Email:manish at south dot du dot ac dot in)</p></div><br class='clearfix'/></div></div>";
print FH2 "</body>";
print FH2 "</html>";



	    
