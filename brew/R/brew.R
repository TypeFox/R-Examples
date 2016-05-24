#
#  brew - text templating for R (www.r-project.org)
#
#  Copyright (C) 2007 Jeffrey Horner
#
#  brew is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

BRTEXT <- 1
BRCODE <- 2
BRCOMMENT <- 3
BRCATCODE <- 4
BRTEMPLATE <- 5
DELIM <- list()
DELIM[[BRTEXT]] <- c('','')
DELIM[[BRCODE]] <- c('<%','%>')
DELIM[[BRCOMMENT]] <- c('<%#','%>')
DELIM[[BRCATCODE]] <- c('<%=','%>')
DELIM[[BRTEMPLATE]] <- c('<%%','%%>')

.bufLen <- 0
.cache <- NULL

# Experimental function for setting the size of internal buffers.
# There's anecdotal evidence that suggests larger buffer sizes
# means faster parsing. 
setBufLen <- function(len=0){
	.bufLen <<- len
	invisible(NULL)
}

brewCache     <- function(envir=NULL) {
	if (missing(envir)) return(.cache)
	.cache <<- envir
	invisible(NULL)
}
brewCacheOn  <- function() brewCache(new.env(hash=TRUE,parent=globalenv()))
brewCacheOff <- function() brewCache(NULL)

`.brew.cached` <- function(output=stdout(),envir=parent.frame()){
	# Only sink if caller passed an argument
	sunk <- FALSE
	if (!missing(output)) {
		sunk <- TRUE
		sink(output)
	}

	text <- get('text')
	brew.cat <- function(from,to) cat(text[from:to],sep='',collapse='')
	.prev.brew.cat <- NULL
	if (exists('.brew.cat',envir=envir)){
		.prev.brew.cat <- get('.brew.cat',pos=envir)
	}
	assign('.brew.cat',brew.cat, envir=envir)

	code <- get('code')
	ret <- try(eval(code,envir=envir))

	# sink() will warn if trying to end the real stdout diversion
	if (sunk && unclass(output) != 1) sink()

	if(!is.null(.prev.brew.cat)){
		assign('.brew.cat',.prev.brew.cat,envir=envir)
	} else {
		rm('.brew.cat',envir=envir)
	}

	invisible(ret)
}

`brew` <-
function(file=stdin(),output=stdout(),text=NULL,envir=parent.frame(),run=TRUE,parseCode=TRUE,tplParser=NULL,chdir=FALSE){

	file.mtime <- canCache <- isFile <- closeIcon <- FALSE
	filekey <- file # we modify file when chdir=TRUE, so keep same cache key

	# Error check input
	if (is.character(file) && file.exists(file)){
		isFile <- closeIcon <- TRUE
		if (chdir || isTRUE(getOption('brew.chdir'))){
			isURL <- length(grep("^(ftp|http|file)://", file)) > 0L
			if(isURL)
				warning("'chdir = TRUE' makes no sense for a URL")
			if(!isURL && (path <- dirname(file)) != ".") {
				owd <- getwd()
				if(is.null(owd))
					warning("cannot 'chdir' as current directory is unknown")
				else on.exit(setwd(owd), add=TRUE)
				setwd(path)
				file <- basename(file)
			}
		}

	} else if (is.character(text) && nchar(text[1]) > 0){
		closeIcon <- TRUE
		icon <- textConnection(text[1])
	} else if (inherits(file,'connection') && summary(file)$"can read" == 'yes') {
		icon <- file
	} else {
		stop('No valid input')
		return(invisible(NULL))
	}

	# Error check output
	if (inherits(output,'connection')){
		if (summary(output)$"can write" != 'yes')
			stop('output connection is not writeable')
	} else if ( !is.character(output) )
		stop('No valid output')

	# Error check env
	if (!is.environment(envir)){
		warning('envir is not a valid environment')
		return(invisible(NULL))
	}

	# Can we use the cache
	if (!is.null(.cache) && isFile && run && is.null(tplParser)){
		canCache <- TRUE
		if (exists(filekey,.cache)){
			file.cache <- get(filekey,.cache)
			file.mtime <- file.info(file)$mtime
			if (file.cache$mtime >= file.mtime){
				brew.cached <- .brew.cached
				environment(brew.cached) <- file.cache$env
				if (!missing(output)) {
					return(brew.cached(output,envir))
				} else {
					return(brew.cached(envir=envir))
				}
			}
		}
	}

	# Not using cache, open input file if needed
	if (isFile) icon <- file(file,open="rt")

	state <- BRTEXT
	text <- code <- tpl <- character(.bufLen)
	textLen <- codeLen <- as.integer(0)
	textStart <- as.integer(1)
	line <- ''
	
	while(TRUE){
		if (!nchar(line)){
			line <- readLines(icon,1)
		   	if (length(line) != 1) break
			line <- paste(line,"\n",sep='')
		}
		if (state == BRTEXT){

			spl <- strsplit(line,DELIM[[BRCODE]],fixed=TRUE)[[1]]

			# Beginning markup found
			if (length(spl) > 1){

				if (nchar(spl[1])) {
					text[textLen+1] <- spl[1]
					textLen <- textLen + 1
				}
				line <- paste(spl[-1],collapse='<%')

				# We know we've found this so far, so go ahead and set up state.
				state <- BRCODE

				# Now let's search for additional markup.
				if (regexpr('^=',spl[2]) > 0){
					state <- BRCATCODE
					line <- sub('^=','',line)
				} else if (regexpr('^#',spl[2]) > 0){
					state <- BRCOMMENT
				} else if (regexpr('^%',spl[2]) > 0){
					if (is.null(tplParser)){
						text[textLen+1] <- '<%'
						textLen <- textLen + 1
					}
					line <- sub('^%','',line)
					state <- BRTEMPLATE
					next
				}

				if (textStart <= textLen) {
					code[codeLen+1] <- paste('.brew.cat(',textStart,',',textLen,')',sep='')
					codeLen <- codeLen + 1
					textStart <- textLen + 1
				}
			} else {
				text[textLen+1] <- line
				textLen <- textLen + 1
				line <- ''
			}
		} else {
			if (regexpr("%%>",line,perl=TRUE) > 0){
				if (state != BRTEMPLATE)
					stop("Oops! Someone forgot to close a tag. We saw: ",DELIM[[state]][1],' and we need ',DELIM[[state]][2])
				spl <- strsplit(line,"%%>",fixed=TRUE)[[1]]
				if (!is.null(tplParser)){
					tpl[length(tpl)+1] <- spl[1]
					# call template parser
					tplBufList <- tplParser(tpl)
					if (length(tplBufList)){
						textBegin <- textLen + 1;
						textEnd <- textBegin + length(tplBufList) - 1
						textLen <- textEnd
						text[textBegin:textEnd] <- tplBufList
					}
					tpl <- character()
				} else {
					text[textLen+1] <- paste(spl[1],'%>',sep='')
					textLen <- textLen + 1
				}
				line <- paste(spl[-1],collapse='%%>')
				state <- BRTEXT
				next
			}
			if (regexpr("%>",line,perl=TRUE) > 0){
				spl <- strsplit(line,"%>",fixed=TRUE)[[1]]
				line <- paste(spl[-1],collapse='%>')

				n <- nchar(spl[1])
				# test  for '-' immediately preceding %> will strip trailing newline from line
				if (n > 0) {
					if (substr(spl[1],n,n) == '-') {
						line <- substr(line,1,nchar(line)-1)
						spl[1] <- substr(spl[1],1,n-1)
					}
					text[textLen+1] <- spl[1]
					textLen <- textLen + 1
				}

				# We've found the end of a brew section, but we only care if the
				# section is a BRCODE or BRCATCODE. We just implicitly drop BRCOMMENT sections
				if (state == BRCODE){
					code[codeLen+1] <- paste(text[textStart:textLen],collapse='')
					codeLen <- codeLen + 1
				} else if (state == BRCATCODE){
					code[codeLen+1] <- paste('cat(',paste(text[textStart:textLen],collapse=''),')',sep='')
					codeLen <- codeLen + 1
				}
				textStart <- textLen + 1
				state <- BRTEXT
			} else if (regexpr("<%",line,perl=TRUE) > 0){
				stop("Oops! Someone forgot to close a tag. We saw: ",DELIM[[state]][1],' and we need ',DELIM[[state]][2])
			} else {
				if (state == BRTEMPLATE && !is.null(tplParser))
					tpl[length(tpl)+1] <- line
				else {
					text[textLen+1] <- line
					textLen <- textLen + 1
				}
				line <- ''
			}
		}
	}
	if (state == BRTEXT){
		if (textStart <= textLen) {
			code[codeLen+1] <- paste('.brew.cat(',textStart,',',textLen,')',sep='')
			codeLen <- codeLen + 1
			textStart <- textLen + 1
		}
	} else {
		stop("Oops! Someone forgot to close a tag. We saw: ",DELIM[[state]][1],' and we need ',DELIM[[state]][2])
	}

	if (closeIcon) close(icon)

	if (run){

		brew.env <- new.env(parent=globalenv())
		assign('text',text,brew.env)
		assign('code',parse(text=code,srcfile=NULL),brew.env)
		brew.cached <- .brew.cached
		environment(brew.cached) <- brew.env

		if (canCache){
			if (file.mtime == FALSE) file.mtime <- file.info(file)$mtime
			assign(filekey,list(mtime=file.mtime,env=brew.env),.cache)
		}

		if (!missing(output)) {
			return(brew.cached(output,envir))
		} else {
			return(brew.cached(envir=envir))
		}
	} else if (parseCode){
		brew.env <- new.env(parent=globalenv())
		assign('text',text,brew.env)
		assign('code',parse(text=code,srcfile=NULL),brew.env)
		brew.cached <- .brew.cached
		environment(brew.cached) <- brew.env
		invisible(brew.cached)
	} else {
		invisible(list(text=text,code=code))
	}
}

.onAttach <- function(library, pkg)
{
	unlockBinding('.bufLen',asNamespace('brew'))
	unlockBinding('.cache',asNamespace('brew'))
}
