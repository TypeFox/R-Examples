autoemail<-function(eadd, sfile, hnote="Exam Results")
{

  ######### eadd   = email address
  #######   sfile = name of file to be sent
  #######   hnote="Exam Results"  (subject line for email)
###########   this program creates the perl script used for sending files to students
  if(missing(hnote)) {  hnote="Exam Results"  }

  
  efile = paste(sep="", sfile, "efile.prl")
  if(file.exists( efile) ) { system(paste(sep=" ", "rm ", efile)) }

  cat(file=efile, "#!/usr/bin/perl", sep="\n", append=TRUE)


  U1 = unlist(strsplit(split="@", as.character(eadd) ))


  toadd = paste(sep="", U1[1] , "\\@", U1[2])
  scont = scan(file=sfile, what="", sep="\n")
  mess = paste(scont, collapse="\\n")
  ##  mess =  "This is a test of the email function."


  sendline = paste(sep="", "sendEmail(\"", toadd, "\",", "\"jonathan.lees\\@unc.edu\", \"",  hnote ,"\"," , "\"", mess ,  "\");")

  cat(file=efile, sendline, sep="\n", append=TRUE)

  cat(file=efile, "# Simple Email Function", sep="\n", append=TRUE)
  cat(file=efile, "# ($to, $from, $subject, $message)", sep="\n", append=TRUE)
  cat(file=efile, "sub sendEmail", sep="\n", append=TRUE)
  cat(file=efile, "{", sep="\n", append=TRUE)
  cat(file=efile, "my ($to, $from, $subject, $message) = @_;", sep="\n", append=TRUE)
  cat(file=efile, "my $sendmail = '/usr/lib/sendmail';", sep="\n", append=TRUE)
  cat(file=efile, "open(MAIL, \"|$sendmail -oi -t\");", sep="\n", append=TRUE)
  cat(file=efile, "print MAIL \"From: $from\\n\";", sep="\n", append=TRUE)
  cat(file=efile, "print MAIL \"To: $to\\n\";", sep="\n", append=TRUE)
  cat(file=efile, "print MAIL \"Subject: $subject\\n\\n\";", sep="\n", append=TRUE)
  cat(file=efile, "print MAIL \"$message\\n\";", sep="\n", append=TRUE)
  cat(file=efile, "close(MAIL);", sep="\n", append=TRUE)
  cat(file=efile, "}", sep="\n", append=TRUE)


  return(efile)

}
