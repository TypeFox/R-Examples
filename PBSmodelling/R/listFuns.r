#writeList------------------------------2009-02-05
# Writes a list to a file in "D" or "P" format.
# Arguments:
#  x        - list to save
#  fname    - file to write list to
#  format   - write list in "D" or "P" format
#  comments - include string as comment at the top of file
#-------------------------------------------ACB/RH
writeList <- function(x, fname="", format="D", comments="") {
	NoComments<-missing(comments)
	comments <- sub("^", "#", comments)
	if (format=="D") {
		dput(x, fname)
		if (file.exists(fname) && !NoComments) {
			output <- scan(fname, what=character(0), sep="\n")
			output <- paste(output, collapse="\n")

			sink(fname)
                        on.exit(expr=sink());

                        #add comments
                        if (any(comments!="#")) {
                                cat(paste(comments, collapse="\n"));
                                cat("\n");
                        }
                        #spit out original output from dput
                        cat(output)
                        cat("\n")
		}
		return(fname)
	}
	if (format=="P") {
		if (!is.list(x) || length(x) == 0)
			stop("x must be a non-empty list.")
		.writeList.P(x, fname, comments)
		return(fname)
	}
	stop(paste("format \"",format,"\" not recognized."))
}
#----------------------------------------writeList


#.writeList.P---------------------------2009-02-10
# Saves list x to disk using "P" format
#-------------------------------------------ACB/RH
.writeList.P <- function( x, fname="", comments, prefix="") {
	if (fname!="") {
                sink(fname)
                on.exit(expr=sink())
        }
	if (!missing(comments)) {
		cat(paste(comments,collapse="\n")); cat("\n")
        }

        # check for list elements that are missing names
	xNames=names(x)
        if (is.null (xNames))
                # names(x) returns NULL if no names are set
                xNames <- names(x) <- paste("X",1:length(x),sep="")
        else {
                # assign X# names for NA/"" cases
                z <- is.na(xNames) | is.element(xNames,"")
                if (any(z))
                        xNames[z] <- names(x)[z] <- paste("X",1:sum(z),sep="")
        }

	#check for errors
	for(i in 1:length(x)) {
		dat=x[[i]]
                # note if a list (A) contains a list (B), is.vector will return
                # true for (B)
		if (!is.vector(dat) && !is.matrix(dat) && class(dat)!="data.frame") next
		#prepare character strings with quotes if spaces exist
		if (is.character(dat) && length(dat)>0) {
			for(j in 1:length(dat)) {
				#only strings with spaces need quotes
				if (typeof(dat[j])=="character") {
					x[[i]][j] <- .addslashes(dat[j])
				}
			}
		}
	}

	#start cat-ing keys and values
	for(i in 1:length(x)) {
		if (prefix != "")
			nam=paste(prefix, xNames[i], sep="$")
		else
			nam=xNames[i]
		dat=x[[i]]
		if (!is.vector(dat) && !is.matrix(dat) && class(dat)!="data.frame" && !is.array(dat))  next
		#print varName
		cat(paste("$", nam, "\n", sep=""))
		# we must handle lists first because is.vector returns true for a list
		if (class(dat) == "list") {
			#list within the list
			.writeList.P (dat, fname="", prefix=nam)
		}
		else if (is.vector(dat)) {
			#print names
			vecNames<-names(dat)
			if (is.null(vecNames)) vecNames=""
			vecNames <- .addslashes(vecNames)
			cat(paste("$$vector mode=\"", typeof(dat), "\" names=", vecNames, "\n", sep=""))
			cat(dat); cat("\n")
		}
		else if( is.matrix( dat ) ) {
			#print colnames
			matColNames <- colnames( dat )
			matRowNames <- rownames( dat )
			if (is.null(matColNames))
				matColNames <- ""
			if (is.null(matRowNames))
				matRowNames <- ""
			matColNames <- .addslashes(matColNames)
			matRowNames <- .addslashes(matRowNames)
			cat(paste("$$matrix mode=\"", typeof( dat ), "\" rownames=", matRowNames, " colnames=", matColNames, " ncol=", ncol( dat ), "\n", sep=""))

			for(j in 1:dim(x[[i]])[1]) {
				cat(dat[j,]); cat("\n")
			}
			
		}
		else if ( is.array(dat) ) {
			d=dim(dat); nr=d[1]; nc=d[2]; nd=length(d) # dimension info
			#get the dimensional names
			dim_names_flat <- c() #single dimension vector
			dim_names <- dimnames( dat ) #note: dimnames may also have a names attribute
			if( is.null( dim_names ) ) {
				dim_names_flat <- ""
			} else {
				#ensure a name exists, assign an index value otherwise
				if( is.null( names( dim_names ) ) )
					names( dim_names ) <- 1:length(dim_names)
				#flatten into a vector c( 1st_dimname, 1st_dim_element_1, 2, ... n, 2nd_dimname, ... )
				for( i in 1:length(dim_names) ) 
					dim_names_flat <- c( dim_names_flat, names(dim_names)[i], dim_names[[i]] )
			}
			dim_names_flat <- .addslashes( dim_names_flat )

			cat(paste("$$array mode=\"", typeof(dat), "\" dim=",.addslashes(d)," byright=FALSE",
				" byrow=TRUE dimnames=", dim_names_flat, "\n", sep="")) 
			dhi=d[3:nd]; nhi=length(dhi) # extra dimensions above 2
			idx=1:nhi; index=letters[10+(idx)]
			ex1=paste(paste("for(",rev(index)," in 1:",dhi[rev(idx)],"){",sep=""),collapse=" ")
			ex2="for(j in 1:nr){"
			ex3=paste("cat(dat[j,,",paste(index,collapse=","),"]); cat(\"\\n\")",sep="")
			ex4=paste(rep("}",nhi+1),collapse="")
			expr=paste(ex1,ex2,ex3,ex4,collapse=" ")
			eval(parse(text=expr)) 
		}
		else if (class(dat)=="data.frame") {
			cat("$$data "); 
			#ncol
			cat("ncol="); cat(dim(dat)[2]); cat(" ");
			#modes
			cat("modes=\"")
			for (j in 1:length(dat)) {
				if (j>1) cat(" ")
					cat(typeof(dat[[j]])) }
			cat("\" ")
			#rownames
			cat("rownames="); cat(.addslashes(rownames(dat))); cat(" ")
			#colnames
			cat("colnames="); cat(.addslashes(colnames(dat))); cat(" ")
			#byrow
			cat("byrow=TRUE"); cat("\n")
			for(j in 1:dim(dat)[1]) {
				for(k in 1:dim(dat)[2]) {
					cat(dat[j,k]); cat(" ") }
				cat("\n") }
		}
	}
}
#-------------------------------------.writeList.P


#readList-------------------------------2014-10-21
# Returns a list object from files originally 
#  formatted in one of the following ways:
#  "D" = created by R functions `dput' or `dump'
#  "R" = R list object as ASCII (e.g., Windows History file)
#  "P" = PBS-formatted file (see `writeList')
#  "C" = Comment-delimited file (e.g., Awatea/ADMB input files)
# Arguments:
#  fname - file to read
#  Change (Anisa Egeli): There is a check to see if the file exists. 
#  This allows try(readList(...), silent=TRUE) to catch the error.
#-------------------------------------ACB/AE/NB/RH
readList <- function(fname) {
	if(!file.exists(fname))
		stop(paste("File", fname, "does not exist."))
	#detect file type
	ff <- scan(fname, what=character(), sep="\n", quiet=TRUE)

	if (any(grepl("^[ \t]*structure", ff))) fileformat <- "D"  # file created by `dput' or `dump'
	else if (any(grepl("^[ \t]*list", ff))) fileformat <- "R"  # R list object
	else if (any(grepl("^[ \t]*\\$", ff)))  fileformat <- "P"  # file in PBS format
	else if (any(grepl("^[ \t]*#", ff)))    fileformat <- "C"  # attempt to interpret comment-delimited data (adapted from PBSawatea)
	else stop("Unknown file format detected (maybe not a list object).")
#browser();return()
	#for(i in 1:length(ff)) {
	#	if (!any(grep("^[ \t]*#", ff[i]))) {
	#		if (any(grep("^[ \t]*structure", ff[i]))) fileformat <- "D"
	#		else if (any(grep("^[ \t]*list", ff[i]))) fileformat <- "R"
	#		else if (any(grep("^[ \t]*\\$", ff[i])))  fileformat <- "P"
	#		else stop("unknown fileformat detected.")
	#		break;
	#	}
	#}
	if (fileformat == "R" || fileformat == "D") 
		return(eval(parse(fname)))
	if (fileformat == "P") 
		return(.readList.P(fname))
	if (fileformat == "C") 
		return(.readList.C(fname))
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~readList

#.readList.P----------------------------2009-02-05
# Read list in "P" format.
#----------------------------------------ACB/AE/RH
.readList.P <- function(fname) {
	#srcfile will be modified, orgfile is untouched and only used for user debug error messages
	srcfile=orgfile=scan(fname, what=character(), sep="\n", quiet=TRUE, blank.lines.skip=FALSE)
	data=list(); j=0; halt=FALSE; str=""
	extendLine <- FALSE    #used for extending a single line into lines with \
	extendLineNumber <- 0  #where a new widget starts - used for error messages

	if (!length(srcfile)) {stop("Input file is empty\n")}
	#print("loop start"); print(date());
	#if comments were striped out earlier, we would lose the line count.
	for(i in 1:length(srcfile)) {
		if (!any(grep("^[[:space:]]*(#.*)?$", srcfile[i]))) {
			srcfile[i] <- .stripComments(srcfile[i])
			#append last string onto new string if applicable
			if (extendLine == TRUE)
				str <- paste(str, srcfile[i], sep=" ")
			else {
				str <- srcfile[i]
				extendLineNumber <- i
			}
			#determine if this string is extended by a \ at the end.
			tmp <- sub('\\\\$', '', str)
			if (tmp==str) #no sub took place
				extendLine = FALSE
			else
				extendLine = TRUE
			str <- tmp
			#parse the line once it is complete (no \)
			if (extendLine == FALSE) {
				j <- j + 1
				data[[j]]<-list(str=str, line.start=extendLineNumber, line.end=i)
			}
		}
	}

	#convert the "data" list into a real list
	varName=varOptions=NULL
	varData=retData=list()    #list to return
	for(i in 1:length(data)) {
		str <- data[[i]]$str

		#varOptions (optional)
		if (substr(str,1,2)=="$$") {
			if (!is.null(varOptions))
				stop("extra $$ line found")
			if (is.null(varName))
				stop("$$ line found before $ line")
			varOptions <-data[[i]]
			varOptions$str = substr(varOptions$str, 3, nchar(varOptions$str)) #remove $$
		}
		#varName
		else if (substr(str,1,1)=="$") {
			if (!is.null(varName)) {
				#save data into the retData list
				listelem <- paste("retData", paste("[[\"", paste(strsplit (varName, "\\$")[[1]], collapse="\"]][[\""), "\"]]", sep=""), sep="")
				savedata <- paste(listelem, " <- .readList.P.convertData(varOptions, varData, fname, orgfile)", sep="")
				eval(parse(text=savedata))
				if (is.null(eval(parse(text=listelem)))) halt<-TRUE
				varName <- varOptions <- NULL
				varData <- list()
			}
			varName <- .trimWhiteSpace(substr(str, 2, nchar(str)))
			if (!any(grep("^[a-zA-Z0-9_.$ ]+$", varName))) {
				.catError(
					paste("Variable name \"", varName,"\" is not valid", sep=""), fname, 
					data[[i]]$line.start, data[[i]]$line.end, 
					orgfile, "readList error"
					)
				halt<-TRUE
			}
			line.start <- data[[i]]$line.start
		}
		else {
			varData[[length(varData)+1]] <-data[[i]]
		}
	}

	#save anything from after
	if (!is.null(varName)) {
		listelem <- paste("retData", paste("[[\"", paste(strsplit(varName, "\\$")[[1]], collapse="\"]][[\""), "\"]]", sep=""), sep="")
		savedata <- paste(listelem, " <- .readList.P.convertData(varOptions, varData, fname, orgfile)", sep="")
		eval(parse(text=savedata))
		if (is.null(eval(parse(text=listelem)))) halt<-TRUE
	}

	if (halt) {stop("Errors were found in the file. Unable to continue\n")}

	return(retData)
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~.readList.P

#.readList.P.convertData----------------2008-03-05
# Helper function to convert data into proper mode.
#-------------------------------------------ACB/RH
.readList.P.convertData <- function(varOptions, varData, fname="", sourcefile=list()) {
	if (is.null(varOptions)) {
		#simple format with no options
		if (length(varData)==0) return(varData)
		else if (length(varData)==1) {
			#just a vector
			return(.autoConvertMode(.convertParamStrToVector(varData[[1]]$str, fname, varData[[1]]$line.start)))
		}
		else {
			#some sort of matrix to parse
			dimSize <- c(length(varData),0) #num of rows
			matData <- c()                  #vector to hold values byrow
			for(i in 1:length(varData)) {
				tmp <- .convertParamStrToVector(varData[[i]]$str, fname, varData[[i]]$line.start)
				if (dimSize[2]==0)
					dimSize[2] <- length(tmp)
				else if (length(tmp)!=dimSize[2]) {
					.catError(paste("Matrix row (line ",varData[[i]]$line.start,") lenght should match first row (line ",varData[[1]]$line.start,") length of ", dimSize[2], sep=""), fname, 
					varData[[i]]$line.start, varData[[i]]$line.end, 
					sourcefile, "readList error")
					return(NULL)
				}
				matData <- append(matData, tmp)
			}
			matData <- .autoConvertMode(matData)
			return(matrix(matData, dimSize[1], dimSize[2], byrow=TRUE))
		}
	}
	#otherwise varOptions was given (in string format)
	#convert it into a list first
	opts <-.getParamFromStr(varOptions$str, fname, varOptions$line.start, varOptions$line.end, sourcefile, .pFormatDefs)
	if (is.null(opts)) stop("Errors were detected")
	if (length(varData)==0) return(eval(parse(text=paste(opts$mode,"()",sep="")))) # return empty typeof

	#flatten all data into a vector (of characters)
	x <- c()
	for(i in 1:length(varData)) {
		#weird things happen if its x[i] <- as.vector(.convert...)
		x <- c(x, .convertParamStrToVector(varData[[i]]$str, fname, varData[[i]]$line.start))
	}
	if(opts$type=="vector") {
		x <- .forceMode(x, opts$mode)
		if (any(opts$names!="")) {
			names(x)<-opts$names
		}
		return(x)
	}
	else if(opts$type=="matrix") {
		x <- .forceMode(x, opts$mode)
		#calculate dims
		nrow <- length(x)/opts$ncol
		if (as.integer(nrow)!=nrow) {
			.catError(paste("Matrix data length [", length(x), "] is not a sub-multiple of ncol [", opts$ncol, "]", sep=""), fname, 
			varOptions$line.start, varData[[length(varData)]]$line.end, 
			sourcefile, "readList error")
			return(NULL)
		}
		#convert to matrix
		mat <- matrix(x, nrow, opts$ncol, byrow=opts$byrow)
		#add colnames
		if (any(opts$colnames!="")) {
			if (length(opts$colnames)!=opts$ncol) {
				.catError(paste("Matrix colnames length [", length(opts$colnames), "] is not equal to ncol [", opts$ncol, "]", sep=""), fname, 
				varOptions$line.start, varData[[length(varData)]]$line.end, 
				sourcefile, "readList error")
				return(NULL)
			}
			colnames(mat)<-opts$colnames
		}
		#add rownames
		if (any(opts$rownames!="")) {
			if (length(opts$rownames)!=nrow) {
				.catError(paste("Matrix rownames length [", length(opts$rownames), "] is not equal to nrow [", nrow, "]", sep=""), fname, 
				varOptions$line.start, varData[[length(varData)]]$line.end, 
				sourcefile, "readList error")
				return(NULL)
			}
			rownames(mat)<-opts$rownames
		}
		return(mat)
	}
	else if(opts$type=="array") {
		x <- .forceMode(x, opts$mode)
		opts$dim <- .convertMode(opts$dim, "numeric")
		if (any(is.na(opts$dim))) {
			.catError("dim values must be numeric", fname, 
			          varOptions$line.start, varData[[length(varData)]]$line.end, 
			          sourcefile, "readList error")
			return(NULL)
		}
		#check dims works
		if (length(x)!=prod(opts$dim)) {
			.catError(paste("dim [product ",prod(opts$dim),"] do not match the length of object [",length(x),"]", sep=""), fname, 
			          varOptions$line.start, varData[[length(varData)]]$line.end, 
			          sourcefile, "readList error")
			return(NULL)
		}
		x=.convertVecToArray(x,opts$dim,byright=opts$byright,byrow=opts$byrow)

		if( all( opts$dimnames == "" ) )
			return( x )

		#restore dimnames
		# example: dimnames(Titanic) -> there are 4 dimensions (class, sex, age, survived)
		# these dimensions have different number of names: i.e. Sex has two: male, female
		# dimnames contains first the name of the element "sex", followed by the labels for the each dimension
		# "male", "female". dim(Titanic)[2] tells us sex only has two labels, so the next label is a dimension name
		# ex: dimnames="Class 1st 2nd 3rd Crew Sex Male Female Age Child Adult Survived No Yes" dim="4 2 2 2"
		dim_name_dimensions <- length( opts$dim )
		dim_names <- list()

		for( i in 1:dim_name_dimensions ) {
			#j points to name of the dimension
			if( i == 1 )
				j <- 1
			else
				j <- sum( opts$dim[ 1:(i-1) ] + 1 ) + 1
			#dim_name_elements are the element names for a particular dimension
			dim_name_elements <- ( j + 1 ) : ( j + opts$dim[ i ] )
			dim_names[[ i ]] <- opts$dimnames[ dim_name_elements ]
			names( dim_names )[ i ] <- opts$dimnames[ j ]
		}
		dimnames( x ) <- dim_names
		return(x)
	}
	else if(opts$type=="data") {
		#check ncol works
		if (length(x)%%opts$ncol>0) {
			.catError(paste("dataframe data length [", length(x), "] is not a sub-multiple of ncol [", opts$ncol, "]", sep=""), fname, 
			varOptions$line.start, varData[[length(varData)]]$line.end, 
			sourcefile, "readList error")
			return(NULL)
		}
		if (opts$ncol != length(opts$colnames)) {
			.catError(paste("Data colnames length [", length(opts$colnames), "] is not equal to ncol [", opts$ncol, "]", sep=""), fname, 
			varOptions$line.start, varData[[length(varData)]]$line.end, 
			sourcefile, "readList error")
			return(NULL)
		}
		if (opts$ncol != length(opts$modes)) {
			.catError(paste("Data modes length [", length(opts$modes), "] is not equal to ncol [", opts$ncol, "]", sep=""), fname, 
			varOptions$line.start, varData[[length(varData)]]$line.end, 
			sourcefile, "readList error")
			return(NULL)
		}
		#calculate nrow
		nrow <- length(x)/opts$ncol
		#break up data into a vector of a list, such that each element represents a column
		dataCols <- list()
		if (opts$byrow) {
			for(i in 1:length(x)) {
				j <- i%%opts$ncol
				if (j==0)
					j <- opts$ncol
				if (length(dataCols)<j)
					dataCols[[j]] <- x[i]
				else
					dataCols[[j]] <- c(dataCols[[j]], x[i])
			}
		}
		else {
			for(i in 1:length(x)) {
				j <- as.integer((i-1)/(length(x)/opts$ncol))+1

				if (length(dataCols)<j)
					dataCols[[j]] <- x[i]
				else
					dataCols[[j]] <- c(dataCols[[j]], x[i])
			}
		}
		#create data.frame and use colnames to refer to each colum
		#the data.frame will be stored as 'ret'
		txt <- "ret <- data.frame("
		for(i in 1:opts$ncol) { #for each column
			#convert into propper mode
			dataCols[[i]] <- .convertMode(dataCols[[i]], opts$modes[i])
			if (i>1)
				txt <- paste(txt, ", ", sep="")
			name <- opts$colnames[i]
			txt <- paste(txt, name, "=dataCols[[", i, "]]", sep="")
		}
		txt <- paste(txt, ")", sep="")
		eval(parse(text=txt))
		#add rownames if any exist
		if (any(opts$rownames!="")) {
			if (length(opts$rownames)!=nrow) {
				.catError(paste("Data rownames length [", length(opts$rownames), "] is not equal to nrow [", nrow, "]", sep=""), fname, 
				varOptions$line.start, varData[[length(varData)]]$line.end, 
				sourcefile, "readList error")
				return(NULL)
			}
			rownames(ret)<-opts$rownames
		}
		return(ret)
	}
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~.readList.P.convertData

#.readList.C----------------------------2014-10-21
# Read ADMB-like input file and create a list.
# Adapted from `readAD' in PBSawatea.
#-----------------------------------------------RH
.readList.C = function(fname) {
	Otxt = readLines(fname) # original text
	otxt = .trimWhiteSpace(Otxt)
	utxt = otxt[!is.element(substring(otxt,1,3),"###")] # use text (remove data not comments)
	xlst = strsplit(utxt,"")
	xlst = xlst[sapply(xlst,length)>0]
	ntxt = sapply(xlst,function(x){paste(clipVector(x,clip="\t",end=2),collapse="")})
	#ntxt = gsub("\\\t\\\t","\t",ntxt)   # internal cleanup
	ntxt = gsub("[[:space:]]+"," ",ntxt)   # internal cleanup
	vlst = list(); vcom=NULL; acom=NULL

	for (i in 1:length(ntxt)) {
		if (substring(ntxt[i],1,1)=="#") {
			expr = paste("vlst[[\"L",pad0(i,3),": Comment\"]] = ",deparse(ntxt[i],width.cutoff=500),sep="")
			comm = .trimWhiteSpace(substring(ntxt[i],2,)) ; names(comm)=i
			acom = c(acom,comm) }
		else {
			vcom = c(vcom,comm)
#if (i==4) {browser();return()}
			if ( all(sapply(strsplit(ntxt[i],""),function(x){grepl("[-0-9. ]",x)}))) # all numeric components
				expr = paste0("vlst[[\"L",pad0(i,3),": Data {",comm,"}\"]] = c(",gsub("[[:space:]]+",",",ntxt[i]),")")
			else
				expr = paste0("vlst[[\"L",pad0(i,3),": Data {",comm,"}\"]] = \"",gsub("\"","\\\\\"",ntxt[i]),"\"")
		}
		eval(parse(text=expr))
	}
	
	# description of variables with inputs
	vdesc = unique(vcom)
	gcomm = acom[!is.element(acom,vdesc)] # general comments
	vcomm = acom[is.element(acom,vdesc)]  # variable comments (may be duplicated)
	nvars = length(vdesc)
	names(vdesc) = paste("v",pad0(1:nvars,3),sep="")
	vars = as.list(vdesc)

	dnam = names(vlst); names(dnam)=1:length(dnam)
	dnam = dnam[grep("Data",dnam)]
	dnam = substring(dnam,13,nchar(dnam)-1)
	
	for (i in names(vdesc)) {
		ivar = vlst[as.numeric(names(dnam[is.element(dnam,vdesc[i])]))]
		ilen = length(ivar)
		if (ilen==1) vars[[i]] = ivar[[1]]
		else if (ilen>=2) {
			jvar=NULL
			for (j in 1:length(ivar))
				jvar=rbind(jvar,ivar[[j]])
			vars[[i]] = jvar
		}
		else stop("Conversion to a list object failed.")
	}
	names(vars)=vdesc[names(vars)]
#browser();return()
	return(vars)
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~readList.C


#unpackList-----------------------------2012-12-06
# Make local/global variables from the components of a named list.
#-------------------------------------------ACB/RH
unpackList <- function(x, scope="L") {
	namx <- names(x); nx <- length(namx);
	if (nx > 0) for (i in 1:nx) {
		if (namx[i] != "") {
			if (scope=="L")
				assign(namx[i], x[[i]], pos=parent.frame(1))
			else if (scope=="P")
				assign(namx[i], x[[i]], envir = .PBSmodEnv)
			else if (scope=="G")
				eval(parse(text="assign(namx[i], x[[i]], envir = .GlobalEnv)"))
		}
	}
	namx[namx != ""]
}
#---------------------------------------unpackList


#packList-------------------------------2013-07-16
# Pack a list in target environment with existing
# objects from parent or user-specified environment.
# NOTE:-------------
# New 'packList' takes advantage of the accessor functions:
#  'tget', 'tcall', and 'tput'.
# Less complicated than former 'packList' (temporarily 
#  available as '.packList.deprecated'), and should be way faster.
#-----------------------------------------------RH
packList=function(stuff, target="PBSlist", value, penv=NULL, tenv=.PBSmodEnv)
{
	if (is.null(penv)) penv = parent.frame() # for a parent envir, need to call this inside the function NOT as an argument
	if (!is.vector(stuff) || !is.character(stuff))
		showAlert("Provide a vector of names denoting objects")
	target = as.character(substitute(target))
	if (target %in% lisp(envir=tenv))
		eval(parse(text=paste("tget(",target,",tenv=tenv)",sep="")))
	else
		eval(parse(text=paste(target,"=list()")))
	if (!missing(value)) {# use explicit value instead of objects
		eval(parse(text=paste(target,"[[\"",stuff,"\"]] = ",paste(deparse(value),collapse="\n"),sep="")))
	} else {
		for (i in stuff)
			eval(parse(text=paste(target,"[[\"",i,"\"]] = tcall(",i,",tenv=penv)",sep="")))
	}
	eval(parse(text=paste("tput(",target,",tenv=tenv)",sep="")))
	invisible()
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~packList


#lisp-----------------------------------2012-12-12
lisp = function(name, pos=.PBSmodEnv, envir=as.environment(pos), all.names=TRUE, pattern){
	ls(name,pos,envir,all.names,pattern)
}
#---------------------------------------------lisp


#.packList.deprecated-------------------2012-12-03
# Pack a list with (i) existing objects using their
# names or (ii) one explicit value.
#-----------------------------------------------RH
.packList.deprecated=function(stuff, target="PBSlist", value, tenv=.PBSmodEnv)
{
	penv = parent.frame() # need to call this inside the function NOT as an argument
	# Deparse bad objects: those that break code (see function 'deparse')
	deparseBO = function(x){
		if (is.list(x) && !is.data.frame(x)) {
			sapply(x,function(x){deparseBO(x)},simplify=FALSE) # recursion through lists within lists
		} else {
			xclass=class(x)
			if (mode(x) %in% c("call","expression","(","function","NULL"))
				x=paste(deparse(x),collapse="")
			class(x)=xclass
			return(x)
		}
	}
	if (!is.vector(stuff) || !is.character(stuff))
		showAlert("Provide a vector of names denoting objects")
	target=deparse(substitute(target))
	target=gsub("^\"","",gsub("\"$","",target)) # strip leading and ending escaped quotes
	endpos=regexpr("[\\[$]",target)-1; if (endpos<=0) endpos=nchar(target)
	base=substring(target,1,endpos)
	if (!exists(base,envir=tenv)) 
		assign(base,list(),envir=tenv)
	if (!missing(value)) { # use explicit value instead of objects
		objet=paste(deparse(value),collapse="\n")
		eval(parse(text=paste(target,"[[\"",stuff,"\"]]=",objet,sep="")),envir=tenv) } #pack explicit value into the list
	else {
		for (s in stuff) {
			if (!exists(s,envir=penv)) next
			eval(parse(text=paste("objet=get(\"",s,"\",envir=penv)",sep=""))) #grab the local object
			if (is.list(objet) && !is.data.frame(objet)) {
				atts=attributes(objet)
				objet=deparseBO(objet)
				natts=setdiff(names(atts),names(attributes(objet)))
				# retain additional attributes of the original list
				if (length(natts)>0) { 
					for (i in natts) attr(objet,i)=atts[[i]] }
				# Reminder: applying original class can cause a display error if 
				# underlying objects (e.g., functions, calls) have been converted to strings.
				#lclass=sapply(objet,class,simplify=FALSE)
				#for(i in 1:length(objet)) attr(objet[[i]],"class")=lclass[[i]]
			}
			objet=paste(deparse(objet),collapse="\n")
			eval(parse(text=paste(target,"[[\"",s,"\"]]=",objet,sep="")),envir=tenv) } #pack into the list
	}
	invisible() }
#-----------------------------.packList.deprecated

#===== THE END ===================================

