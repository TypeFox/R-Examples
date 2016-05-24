# copyright (C) 2004 Witold Eryk Wolski.

##==================================================================
## multipartpost
##==================================================================

# changed to post binary files (including \0) by Andreas Westfeld, 2/2/2006
# added cookies by Andreas Westfeld, 13/8/2007
# 
postToHost <- function(host, path, data.to.send, referer, port=80, ua, accept,
	accept.language, accept.encoding, accept.charset, contenttype, cookie)
{
	if(missing(path))
		path <- "/"
	if(missing(ua))
		ua <- "Mozilla/5.0 (X11; U; Linux i686; en-US; rv:1.8.1.5) Gecko/20070719 Iceweasel/2.0.0.5 (Debian-2.0.0.5-2)"
	if(missing(referer))
		referer <- NULL
	if(missing(accept))
		accept <- "text/xml,application/xml,application/xhtml+xml,text/html;q=0.9,text/plain;q=0.8,image/png,*/*;q=0.5"
	if(missing(accept.language))
		accept.language <- "de,de-de;q=0.8,en-us;q=0.5,en;q=0.3"
	if(missing(accept.encoding))
		accept.encoding <- "gzip,deflate"
	if(missing(accept.charset))
		accept.charset <- "ISO-8859-1,utf-8;q=0.7,*;q=0.7"
	if(missing(contenttype))
		contenttype <- "application/octet-stream"
	if(missing(cookie))
		cookie <- NULL
	if(missing(data.to.send))
		stop("No data to send provided")
	if(!inherits(data.to.send,"list"))
		stop("Data to send have to be a list")

	dc <- 0; #counter for strings
	#make border
	xx <- as.integer(runif(29, min=0, max=9))
	bo <- paste(xx, collapse="")
	bo <- paste(paste(rep("-", 29), collapse=""), bo, sep="")

	header <- NULL
	header <- c(header,paste("POST ", path, " HTTP/1.1\r\n", sep=""))
	if (port==80)
		header <- c(header,paste("Host: ", host, "\r\n", sep=""))
	else
		header <- c(header,paste("Host: ", host, ":", port,
			"\r\n", sep=""))
	header <- c(header,paste("User-Agent: ", ua, "\r\n", sep=""))
	header <- c(header,paste("Accept: ", accept, "\r\n", sep=""))
	header <- c(header,paste("Accept-Language: ", accept.language,
		"\r\n", sep=""))
	header <- c(header,paste("Accept-Encoding: ", accept.encoding,
		"\r\n", sep=""))
	header <- c(header,paste("Accept-Charset: ", accept.charset,
		"\r\n", sep=""))

	header <- c(header,"Connection: close\r\n")
	if (!is.null(referer))
		header <- c(header,paste("Referer: ", referer, "\r\n", sep=""))
	if (!is.null(cookie))
		header <- c(header,paste("Cookie: ", cookie, "\r\n", sep=""))
	header <- c(header,paste("Content-Type: multipart/form-data; boundary=",substring(bo,3),"\r\n",sep=""))

	mcontent <- NULL # keeps the content.

	for(x in 1:length(data.to.send)) {
		val <- data.to.send[[x]]
		key <- names(data.to.send)[x]
		if (typeof(val)=="list") {
			ds <- c(charToRaw(sprintf("%s\r\nContent-Disposition: form-data; name=\"%s\"; filename=\"%s\"\r\nContent-Type: %s\r\n\r\n", bo, key, val$filename, contenttype)), val$object, charToRaw("\r\n"))
		} else {
			ds <- charToRaw(sprintf("%s\r\nContent-Disposition: form-data; name=\"%s\"\r\n\r\n%s\r\n", bo,as.character(key),as.character(val)))
		}
		dc <- dc + length(ds)
		mcontent <- c(mcontent,ds)
	}

	dc <- dc + length(charToRaw(bo))+4;
	header <- c(header,paste("Content-Length: ",dc,"\r\n\r\n",sep=""))
	mypost <- c(charToRaw(paste(header, collapse="")),mcontent,
		charToRaw(paste(bo,"--\r\n",sep="")))
	rm(header,mcontent)

	scon <- socketConnection(host=host,port=port,open="a+b",blocking=TRUE)
	writeBin(mypost, scon, size=1)

	output <- character(0)
	#start <- proc.time()[3]
	repeat{
		ss <- rawToChar(readBin(scon, "raw", 2048))
		output <- paste(output,ss,sep="")
		if(regexpr("\r\n0\r\n\r\n",ss)>-1) break
		if(ss == "") break
	}
	close(scon)
	return(output)
}

##================================================================
## GET
##================================================================

getToHost <- function(host,path,referer,port=80)
{
  if(missing(path))
    path<-"/"
  if(missing(referer))
    referer<-""
  
  fp <- make.socket(host=host, port=port,server=FALSE)
  header <- character(0)
  header <- c(header,paste("GET ",path," HTTP/1.1\n",sep=""))
  header <- c(header,paste("Host: ",host,"\n",sep=""))
  header <- c(header,"Connection: close\n")
  header <- c(header,paste("Referer: ",referer,"\n",sep=""))
  header <- c(header,"User-Agent: Mozilla/4.05C-SGI [en] (X11; I; IRIX 6.5 IP22)\n")
  header <- c(header,"Accept: */*\n\n")
  tmp<-paste(header,collapse="")

  write.socket(fp,tmp)
  output <- character(0)
  #int <- 1
  repeat{
    ss <- read.socket(fp,loop=FALSE)
    output <- paste(output,ss,sep="")
    if(regexpr("\r\n0\r\n\r\n",ss)>-1) break
    if (ss == "") break
  }
  close.socket(fp)
  return(output)
}

getToHost2 <- function(host,path,referer,port=80)
{
	if(missing(path))
		path <- "/"
	#if(missing(ua))
		ua <- "User-Agent: Mozilla/4.05C-SGI [en] (X11; I; IRIX 6.5 IP22)\n"
	if(missing(referer))
		referer <- ""
	header <- NULL
	header <- c(header,paste("GET ",path," HTTP/1.1\n",sep=""))
	header <- c(header,paste("Host: ",host,"\n",sep=""))
	header <- c(header,"Connection: close\n")
	header <- c(header,paste("Referer: ",referer,"\n",sep=""))
	header <- c(header,ua)
	header <- c(header,"Accept: */*\n\n")

	scon <- socketConnection(host=host,port=port,open="a+b",blocking=TRUE)
	writeBin(charToRaw(paste(header,collapse="")), scon, size=1)

	output <- character(0)
	repeat{
		ss <- rawToChar(readBin(scon, "raw", 2048))
		output <- paste(output,ss,sep="")
		if(regexpr("\r\n0\r\n\r\n",ss)>-1) break
		if(ss == "") break
	}
	close(scon)
	return(output)
}

##================================================================
## Simple Post
##================================================================

simplePostToHost <- function(host, path, datatosend, referer, contenttype,
				port=80, maxlen=131063L) {
	if(!missing(datatosend))
		lengthdatatosend <- length(charToRaw(datatosend))
	else {
		datatosend<-character(0)
		lengthdatatosend <- 0
	}
	if(missing(path)) path<-"/"
	if(missing(referer)) referer <- ""
	if(missing(contenttype))
		contenttype <- "application/x-www-form-urlencoded"
	#make the header
	header<-character(0)
	header<-c(header, paste("POST ", path, " HTTP/1.1\n", sep=""))
	header<-c(header, paste("Host: ", host, "\n", sep=""))
	header<-c(header, paste("Referer: ", referer, "\n", sep=""))
	header<-c(header, paste("Content-Type: ", contenttype, "\n", sep=""))
	header<-c(header, paste("Content-Length: ",
				lengthdatatosend, "\n", sep=""))
	header<-c(header, "Connection: Keep-Alive\n\n")
	#add the data.
	header <- paste(c(header,datatosend,"\n"), collapse="")
	##establish the connection.
	fp <- make.socket(host=host, port=port,server=FALSE)
	write.socket(fp,header)
	output <- character(0)
	# read as long as maxlen bytes have been returned
	repeat {
		ss <- read.socket(fp, maxlen, loop=FALSE)
		output <- paste(output, ss, sep="")
		if (length(charToRaw(ss)) < maxlen) break
		# if(regexpr("\r\n0\r\n\r\n", ss)>-1) break
	}
	close.socket(fp)
	return(output)
}
