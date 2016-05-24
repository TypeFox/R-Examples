`blast` <-
function(X,program="blastn",database="nr",entrezquery="none",nb=5,oot=35)
{
  Y <- NULL
  listRID <- NULL
  
  if(class(X)=="DNAbin") Y <- lapply(as.character(X),paste,collapse="")

  if(class(X)!="DNAbin")
  {
    for (i in 1:length(X))
    {
	GB <- getGB(X[i])
	mol_type <- GB['INSDSet/INSDSeq/INSDSeq_moltype/#']
	if((mol_type=="DNA"|mol_type=="RNA")&(program=="tblastn"|program=="blastp")) stop("change the \"program\" parameter for nucleotide query")
	if((mol_type=="AA")&(program=="blastn"|program=="blastx"|program=="tblastx")) stop("change the \"program\" parameter for protein query")
	Y <- c(Y,X[i])
    }
  }
  
  for(i in 1:length(Y)) listRID <- c(listRID,submitBLAST(Y[i],database=database,entrezquery=entrezquery,program=program,nb=nb))	

  if(class(X)!="DNAbin") names(listRID) <- Y else  names(listRID) <- names(X)

  hit <- list(NULL)
  for (i in 1:length(listRID))
  {
    Sys.sleep(5)
    #out <- outBLAST(listRID[i])

    t <- 5
    blastout <- scan(paste("http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?RID=",listRID[i],"&SHOW_OVERVIEW=no&FORMAT_TYPE=XML&ALIGNMENTS=0&NCBI_GI=yes&CMD=Get",sep=""),what="",sep="\n",quiet=TRUE)
    while (blastout[1] != "<?xml version=\"1.0\"?>" & t <= oot)
    {
      Sys.sleep(t)
      blastout <- scan(paste("http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?RID=",listRID[i],"&SHOW_OVERVIEW=no&FORMAT_TYPE=XML&ALIGNMENTS=0&NCBI_GI=yes&CMD=Get",sep=""),what="",sep="\n",quiet=TRUE)
      t <- t + 5
    }
    if(blastout[1] == "<?xml version=\"1.0\"?>") blastout <- xml(c("<root>",blastout[3:length(blastout)],"</root>")) else blastout <- "out of time"

    blasthit <- sortBLAST(blastout)

    hit[[i]] <- blasthit
  }
names(hit) <- names(listRID)
hit
}

