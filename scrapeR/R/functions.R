# ================
# HELPER FUNCTIONS
# ================

# A function to parse headers from HTML/XML source code
headerParse<-function(source) {
	headers<-sapply(sapply(source,substr,start=1,stop=10000),gsub,pattern="(.+)\r\n\r\n+<.+",replacement="\\1")
	headers<-strsplit(headers,"\r\n")
	
	nixThese<-lapply(headers,function(x){which(x=="")})
	
	for(i in 1:length(nixThese)) {
		if(length(nixThese[[i]])!=0) {
			headers[[i]]<-headers[[i]][-nixThese[[i]]]
		}
	}
	#headers<-lapply(headers,unlist,use.names=FALSE)
	
	# names
	valNames<-lapply(headers,gsub,pattern="^(.+):\\s+.*",replacement="\\1")
	valNames<-lapply(valNames,gsub,pattern="HTTP/[0-9]{1}\\.[0-9]+\\s+[0-9]{3}\\s+.+",replacement="statusCode")
	
	# values
	vals<-lapply(headers,gsub,pattern="^.+:\\s+(.*)",replacement="\\1")
	vals<-lapply(vals,gsub,pattern="HTTP/[0-9]{1}\\.[0-9]+\\s+([0-9]{3})\\s+.+",replacement="\\1")
	
	for(i in 1:length(headers)) {
		names(vals[[i]])<-valNames[[i]]
	}
	vals
}

# A function to retrieve urls contained in any "Location" headers
urlRedirects<-function(headers) {
	contentLocationHeaders<-sapply(lapply(headers,names),grep,pattern="Content-Location")
	locationHeaders<-sapply(lapply(headers,names),grep,pattern="^Location$")
	if(any(sapply(contentLocationHeaders,length)>0)) {
		these<-which(sapply(contentLocationHeaders,length)>0)
		for(i in these) {
			nixThese<-match(contentLocationHeaders[[i]],locationHeaders[[i]])
			locationHeaders[[i]]<-locationHeaders[[i]][-nixThese]
		}
	}
	return<-vector("list",length=length(headers))
	names(return)<-names(headers)
	for(i in 1:length(headers)) {
		if(length(locationHeaders[[i]])>0) {
			urls<-headers[[i]][locationHeaders[[i]]]
			names(urls)<-NULL
			return[[i]]<-urls
		}
	}
	return
}

# A function to strip headers from HTML/XML source code
stripHeaders<-function(source,.encoding) {
	if(identical(.encoding,integer())) {
		startStopPositions<-gregexpr(pattern="^(.+?\r\n\r\n+)<",source)
	} else {
		startStopPositions<-gregexpr(pattern="^(.+?\r\n\r\n+)<",source,useBytes=T)
	}
	
	start<-unlist(lapply(startStopPositions,attr,which="match.length"))
	end<-sapply(source,nchar)
	mapply(substr,x=source,start=start,stop=end)
}

# This function takes a vector and cuts it into chunks of size x or smaller
sliceIt<-function(vec,x) {
	if(length(vec)<=x) {
		list<-vector("list")
		list[[1]]<-vec
	} else {
		result<-cut(1:length(vec),ceiling(length(vec)/x),labels=FALSE)
		list<-vector("list",length=length(unique(result)))
		for(i in unique(result)) {
			list[[match(i,unique(result))]]<-vec[which(result==i)]
		}
	}
	list
}

# =================
# PRIMARY FUNCTIONS
# =================
scrape<-function(url=NULL,object=NULL,file=NULL,chunkSize=50,maxSleep=5,
	userAgent=unlist(options("HTTPUserAgent")),follow=FALSE,headers=TRUE,
	parse=TRUE,isXML=FALSE,.encoding=integer(),verbose=FALSE) {
	
	# Require libraries
	require(XML)
	require(RCurl)
	
	# Perform a bunch of validity checks up front...
	if(is.null(url) && is.null(file) && is.null(object)) {
		stop("You must supply either a vector of URLs, an R object name, or a vector of file paths.  Function call terminated.")
	} else if(!is.null(url) && !is.null(object)) {
		stop("You must supply either a vector of URLs or an R object name, but not both.  Function call terminated.")
	} else if(!is.null(url) && !is.null(file)) {
		stop("You must supply either a vector of URLs or a vector of file paths, but not both.  Function call terminated.")
	} else if(!is.null(object) && !is.null(file)) {
		stop("You must supply either an R object name or a vector of file paths, but not both.  Function call terminated.")
	}
	
	if(!is.null(url)) {
		if(!all(sapply(url,is.character))) {
			stop("At least one of the values in the supplied URL vector is not a character string.  Function call terminated.")
		}
	} else if(!all(sapply(file,is.null))) {
		if(!all(sapply(file,is.character))) {
			stop("At least one of the values in the supplied vector of file paths is not a character string.  Function call terminated.")
		}
	}
	
	# Decide between whether using URL, object name, or file.  First, if URL...
	if(!is.null(url) && is.null(object) && is.null(file)) {
		
		# Check the chunkSize parameter.  It must be >0
		chunkSize<-as.numeric(chunkSize)
		if(chunkSize<=0) {
			chunkSize<-1
		}
		
		# Check the maxSleep parameter.  It must be >0
		maxSleep<-as.numeric(maxSleep)
		if(maxSleep<=0) {
			maxSleep<-1
		}
		
		# Chunk the URLs
		chunks<-sliceIt(url,chunkSize)
		if(verbose) {
			cat(paste("The vector of URLs has been split into",length(chunks),"chunk(s).\n",sep=" "))
		}
		
		# Prep the lists that will hold the header values, redirection URLs, and the page source code
		hdrs<-rep(list(NULL),length(url))
		names(hdrs)<-url

		sourceCode<-as.list(rep(list(NULL),length(url)))
		names(sourceCode)<-url
		
		redirectURL<-rep(list(NULL),length(url))

		if(follow && headers) {	# follow==T and headers==T
		
			# Set up the special CURL request parameters
			#curl=getCurlHandle()
			#curlSetOpt( .opts = list(header=TRUE,useragent=userAgent,followlocation=TRUE),curl=curl)
				
			for(i in 1:length(chunks)) {	# for each chunk...
				if(verbose) {
					cat(paste("\tProcessing chunk ",i," of ",length(chunks),"...\n",sep=""))
				}
				theseSlots<-match(chunks[[i]],url)	# match these URLs in the primary url vector
				
				# Grab source code, follow URL redirects, and include headers
				sourceCode[theseSlots]<-getURLAsynchronous(chunks[[i]],header=TRUE,useragent=userAgent,followlocation=TRUE,.encoding=.encoding)
				
				# Progress reports, if desired, and a nap
				if(verbose) {
					cat(paste("\tProcessing chunk ",i," of ",length(chunks)," -- ",round(i/length(chunks)*100,2),"% complete.\n",sep=""))
				}
				if(i==length(chunks)) {
					if(verbose) {
						cat(paste("\t------------------------------------------------------\n",sep=""))
					}
				} else {
					sleepThis<-ceiling(runif(1,min=0,max=maxSleep))
					if(verbose) {
						cat(paste("\tNow sleeping for ",sleepThis," seconds.\n",sep=""))
						cat(paste("\t------------------------------------------------------\n",sep=""))
					}
					Sys.sleep(sleepThis)
				}
			}
			# Update the source code, header, and URL redirect lists
			hdrs<-headerParse(sourceCode)
			redirectURL<-urlRedirects(hdrs)
			sourceCode<-stripHeaders(sourceCode,.encoding)

		} else if(!follow && headers) {	# follow==F and headers==T
			for(i in 1:length(chunks)) {	# for each chunk...
				if(verbose) {
					cat(paste("\tProcessing chunk ",i," of ",length(chunks),"...\n",sep=""))
				}
				theseSlots<-match(chunks[[i]],url)	# match these URLs in the primary url vector
				
				# Grab source code and include headers
				sourceCode[theseSlots]<-getURLAsynchronous(chunks[[i]],useragent=userAgent,header=TRUE,.encoding=.encoding)
				
				# Progress reports, if desired, and a nap
				if(verbose) {
					cat(paste("\tProcessing chunk ",i," of ",length(chunks)," -- ",round(i/length(chunks)*100,2),"% complete.\n",sep=""))
				}
				if(i==length(chunks)) {
					if(verbose) {
						cat(paste("\t------------------------------------------------------\n",sep=""))
					}
				} else {
					sleepThis<-ceiling(runif(1,min=0,max=maxSleep))
					if(verbose) {
						cat(paste("\tNow sleeping for ",sleepThis," seconds.\n",sep=""))
						cat(paste("\t------------------------------------------------------\n",sep=""))
					}
					Sys.sleep(sleepThis)
				}
			}
			
			# Update the source code, and header lists
			hdrs<-headerParse(sourceCode)
			sourceCode<-stripHeaders(sourceCode,.encoding)
			
		} else if(follow && !headers) {	# follow==T and headers==F
			# First, build a temporary headers holder, since it's ultimately going to be discarded
			hdrsTemp<-rep(list(NULL),length(url))
			
			for(i in 1:length(chunks)) {	# for each chunk...
				if(verbose) {
					cat(paste("\tProcessing chunk ",i," of ",length(chunks),"...\n",sep=""))
				}
				theseSlots<-match(chunks[[i]],url)	# match these URLs in the primary url vector
				
				# Grab source code, follow URL redirects, and include headers
				sourceCode[theseSlots]<-getURLAsynchronous(chunks[[i]],useragent=userAgent,followlocation=TRUE,header=TRUE,.encoding=.encoding)
				
				# Progress reports, if desired, and a nap
				if(verbose) {
					cat(paste("\tProcessing chunk ",i," of ",length(chunks)," -- ",round(i/length(chunks)*100,2),"% complete.\n",sep=""))
				}
				if(i==length(chunks)) {
					if(verbose) {
						cat(paste("\t------------------------------------------------------\n",sep=""))
					}
				} else {
					sleepThis<-ceiling(runif(1,min=0,max=maxSleep))
					if(verbose) {
						cat(paste("\tNow sleeping for ",sleepThis," seconds.\n",sep=""))
						cat(paste("\t------------------------------------------------------\n",sep=""))
					}
					Sys.sleep(sleepThis)
				}
			}
			# Update the source code, temporary headers, and URL redirect lists
			hdrsTemp<-headerParse(sourceCode)
			redirectURL<-urlRedirects(hdrsTemp)
			sourceCode<-stripHeaders(sourceCode,.encoding)
			rm(hdrsTemp)	# Nix the temporary headers holder
		} else if(!follow && !headers){	# follow==F and headers==F
			for(i in 1:length(chunks)) {	# for each chunk...
				if(verbose) {
					cat(paste("\tProcessing chunk ",i," of ",length(chunks),"...\n",sep=""))
				}
				theseSlots<-match(chunks[[i]],url)	# match these URLs in the primary url vector
				
				# Grab source code
				sourceCode[theseSlots]<-getURLAsynchronous(chunks[[i]],useragent=userAgent,.encoding=.encoding)
				
				# Progress reports, if desired, and a nap
				if(verbose) {
					cat(paste("\tProcessing chunk ",i," of ",length(chunks)," -- ",round(i/length(chunks)*100,2),"% complete.\n",sep=""))
				}
				if(i==length(chunks)) {
					if(verbose) {
						cat(paste("\t------------------------------------------------------\n",sep=""))
					}
				} else {
					sleepThis<-ceiling(runif(1,min=0,max=maxSleep))
					if(verbose) {
						cat(paste("\tNow sleeping for ",sleepThis," seconds.\n",sep=""))
						cat(paste("\t------------------------------------------------------\n",sep=""))
					}
					Sys.sleep(sleepThis)
				}
			}
		}
		
		if(parse) {	# are we parsing these results?
			if(isXML) {	# are these results well-formatted XML?
				discard<-capture.output(returnThis<-sapply(sourceCode,xmlParse,asText=TRUE))
			} else {
				discard<-capture.output(returnThis<-sapply(sourceCode,htmlParse,asText=TRUE))
			}
		} else {
			if(!is.list(sourceCode)) {
				returnThis<-as.list(sourceCode)
			}
		}
		rm(sourceCode)
	}
	
	# Otherwise, if object name...
	if(is.null(url) && !is.null(object) && is.null(file)) {
		obj<-get(object)	# fetch object
		if(length(obj)==1) {
			temp<-vector("list",length=1)
			temp[[1]]<-obj
			obj<-temp
			rm(temp)
		}
		redirectURL<-lapply(obj,attr,which="redirect.url")	# fetch any redirect.url attributes
		hdrs<-lapply(obj,attr,which="headers")	# fetch any headers attributes
		
		if(parse) {	# are we parsing these results?
			if(isXML) {	# are these results well-formatted XML?
				discard<-capture.output(returnThis<-sapply(obj,xmlParse,asText=TRUE))
			} else {
				discard<-capture.output(returnThis<-sapply(obj,htmlParse,asText=TRUE))
			}
		} else {
			returnThis<-obj
		}
		rm(obj)
	}
	
	# Otherwise, if file path...
	if(is.null(url) && is.null(object) && !is.null(file)) {
		
		# The redirect.url and headers attributes do not apply here
		redirectURL<-NULL
		hdrs<-NULL
		
		if(parse) {	# are we parsing these results?
			if(isXML) {	# are these results well-formatted XML?
				discard<-capture.output(returnThis<-sapply(file,xmlParse))
			} else {
				discard<-capture.output(returnThis<-sapply(file,htmlParse))
			}
		} else {
			stop("If you are supplying a vector of local file paths, you must set parse=TRUE or else nothing will happen.  Function call terminated.")
		}
	}
	
	# If there are headers to work with, then attach them as attributes to the return value
	if(any(!sapply(hdrs,is.null))) {
		these<-which(!sapply(hdrs,is.null))
		for(i in these) {
			attr(returnThis[[i]],which="headers")<-hdrs[[i]]
		}
	}
	
	# If there are redirection URLs to work with, then attach them as attributes to the return value
	if(any(!sapply(redirectURL,is.null))) {
		these<-which(!sapply(redirectURL,is.null))
		for(i in these) {
			attr(returnThis[[i]],which="redirect.url")<-redirectURL[[i]]
		}
	}
	
	returnThis
}
