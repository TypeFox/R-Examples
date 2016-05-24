`outBLAST` <-
function(RID,oot=oot)
{
  t <- 5
  blastout <- scan(paste("http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?RID=",RID,"&SHOW_OVERVIEW=no&FORMAT_TYPE=XML&ALIGNMENTS=0&NCBI_GI=yes&CMD=Get",sep=""),what="",sep="\n",quiet=TRUE)
  while (blastout[1] != "<?xml version=\"1.0\"?>" & t <= oot)
  {
    Sys.sleep(t)
    blastout <- scan(paste("http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?RID=",RID,"&SHOW_OVERVIEW=no&FORMAT_TYPE=XML&ALIGNMENTS=0&NCBI_GI=yes&CMD=Get",sep=""),what="",sep="\n",quiet=TRUE)
    t <- t + 5
  }
  if(blastout[1] == "<?xml version=\"1.0\"?>") blastout <- xml(c("<root>",blastout[3:length(blastout)],"</root>")) else blastout <- "out of time"
  blastout
}

