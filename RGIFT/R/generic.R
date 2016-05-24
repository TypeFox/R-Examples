#
#Some generic functions
#

GIFTcomment<-function(txt)
{
	cat("//", txt, sep="")
}


GIFTcategory<-function(catname)
{
	cat("\n\n$CATEGORY: ", catname, "\n\n")
}

GIFTQName<-function(qname)
{
	cat("::", qname, "::", sep="")
}

GIFTparse<-function(txt)
{
#Special characters are: ~ / $ { }

spchar<-c("~", "/", "$", "{", "}") 

for(c in spchar)
{
	txt<-gsub(paste("\\", c, sep=""),paste("\\\\", c, sep=""), txt)
}

	return(txt)
}
