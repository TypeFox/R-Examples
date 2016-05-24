print.RMark.version <- function()
{ library(help=RMark)$info[[1]] -> version
	version <- version[pmatch("Version",version)]
	if(!is.null(version))
	{
		um <- strsplit(version," ")[[1]]
  	    version <- um[nchar(um)>0][2]
	}
	hello <- paste("This is RMark ",version,"\n",sep="")
	packageStartupMessage(hello)
}

.onAttach <- function(...) { 
	print.RMark.version()
	checkForMark()
}

create_markpath=function()
{
	markpath=Sys.which("mark.exe")
	if(markpath!="")
	{
		markpath="mark.exe"
        return(markpath)
	}
	if(!exists("MarkPath"))
	{
		markpath=c("c:/Program Files/Mark","c:/Program Files (x86)/Mark")
		markpath=markpath[file.exists(markpath)]
	}else
	{
		if(substr(MarkPath,nchar(MarkPath),nchar(MarkPath))%in%c("\\","/")) MarkPath=substr(MarkPath,1,(nchar(MarkPath)-1))
		markpath=MarkPath
	}
	markstrings=c("mark.exe","mark32.exe","mark64.exe")
	if(length(markpath)!=0)
	{
		markpath=as.vector(sapply(markpath,function(x)paste(x,markstrings,sep="/")))
		which.exists=file.exists(markpath)
	}else
		which.exists=rep(FALSE,3)
	if(any(which.exists)) 
	{
		if(which.exists[1])
		   markpath=shQuote(markpath[1])
	    else
	    {
		   if(R.Version()$arch=="x86_64")
		   {
			   if(which.exists[3])
				   markpath=shQuote(markpath[3])
			   else
			   {
				   warning("\n Warning:mark64.exe does not exist. Using mark32.exe\n")
				   markpath=shQuote(markpath[2])
			   }
		   } else
		   {
			   if(which.exists[2])
				   markpath=shQuote(markpath[2])
			   else
			   {
				   warning("\n Warning:mark32.exe does not exist. Using mark64.exe\n")
				   markpath=shQuote(markpath[3])
			   }
	       }
	   }
    } else
	{
		if(exists("MarkPath"))message("no mark.exe found in specified MarkPath location. Looking for an exe in operating system Path.\n")
		if(!exists("markpath") || length(markpath)>1)
		{
			inPath=Sys.which(markstrings)!=""
			if(inPath[1])
				markpath=shQuote(markstrings[1])
			else
			if(inPath[3]&R.Version()$arch=="x86_64")
				markpath=shQuote(markstrings[3])
			else
			if(inPath[2])
				markpath=shQuote(markstrings[2])
			else
				markpath=NULL
		}
   }
return(markpath)
}



checkForMark<-function()
{
	if(R.Version()$os=="mingw32")
	{
	   markpath=create_markpath()
	   if(is.null(markpath))
	   {
		   cat("Warning: Software mark.exe,mark32.exe or mark64.exe not found in path or in c:/Program Files/mark or c:/Program Files (x86)/mark\n. It is available at http://www.cnr.colostate.edu/~gwhite/mark/mark.htm\n")
	       cat('         If you have mark.exe, you will need to set MarkPath object to its location (e.g. MarkPath="C:/Users/Jeff Laake/Desktop"')
       }
   }else
	   if(exists("MarkPath")) 
       {
	      isep="/"
	      if(substr(MarkPath,nchar(MarkPath),nchar(MarkPath))%in%c("\\","/")) isep=""
  		  MarkPath=paste(MarkPath,"mark",sep=isep)
	      if(!file.exists(MarkPath)) 
		     cat(paste("mark executable cannot be found at specified MarkPath location:",MarkPath,"\n"))		
        } else
        {
 	       if(Sys.which("mark")=="")
	       {
		       cat("Warning: Software mark not found in path.\n")
		       cat('         If you have mark executable, you will need to set MarkPath object to its location (e.g. MarkPath="C:/Users/Jeff Laake/Desktop"')
	       }  
	    }
	invisible()
}


