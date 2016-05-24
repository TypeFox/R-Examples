`submitBLAST` <-
function(ID,database,program,entrezquery,nb)
{
	if(entrezquery!="none") blastin <- scan(paste("http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?QUERY=",ID,"&DATABASE=",database,"&HITLIST_SIZE=",nb,"&FILTER=L&EXPECT=10&PROGRAM=",program,"&CLIENT=web&SERVICE=plain&NCBI_GI=on&PAGE=Nucleotides&ENTREZ_QUERY=",entrezquery,"&CMD=Put",sep=""),what="raw")
	if(entrezquery=="none") blastin <- scan(paste("http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?QUERY=",ID,"&DATABASE=",database,"&HITLIST_SIZE=",nb,"&FILTER=L&EXPECT=10&PROGRAM=",program,"&CLIENT=web&SERVICE=plain&NCBI_GI=on&PAGE=Nucleotides&CMD=Put",sep=""),what="raw")
	RID <- blastin[rev(grep("RID",blastin))[1] + 2]
	RID
}

