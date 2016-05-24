`getGB` <-
function(seqID)
{
	URL <- paste("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=",seqID,"&rettype=gbc",sep = "")
	GBFile <- scan(file = URL, what = "", sep = "\n", quiet = TRUE)
	GB <- xml(c("<root>",GBFile[3:length(GBFile)],"</root>"))
	GB
}

