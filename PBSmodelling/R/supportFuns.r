#================================================|
#                 PBS Modelling                  |
# Support functions (in alphabetic order)--------|
# Authors:                                       |
#  Jon T. Schnute <Jon.SchnuteJ@dfo-mpo.gc.ca>   |
#  Alex Couture-Beil <alex@mofo.ca>              |
#  Rowan Haigh <Rowan.HaighR@dfo-mpo.gc.ca>      |
#  Anisa Egeli <>                                |
#================================================|


#calcFib--------------------------------2006-08-28
# Calculate a vector containing fibonacci numbers
# Arguments:
#  len    - return the last "len" calculated numbers 
#  n      - calculate the nth number
#  method - use .C, .Call, R code, or closed form
#----------------------------------------------JTS
calcFib <- function(n, len=1, method="C") {
	if (n<0)
		return(NA)
	if (len>(n+1))
		len <- (n+1)

	switch(casefold(method),
	       c=.fibC(n,len),
	       call=.fibCall(n,len),
	       r=.fibR(n,len),
	       closed=.fibClosedForm(n,len)
	       )
}

.fibCall <- function(n, len=1) {
	retArr <- numeric(len)
	out <- .Call("fibonacci2", as.integer(n), as.integer(len), PACKAGE="PBSmodelling")
	return(out) }

.fibC <- function(n, len=1) {
	retArr <- numeric(len)
	out <- .C("fibonacci", as.integer(n), as.integer(len), as.numeric(retArr), PACKAGE="PBSmodelling")
	x <- out[[3]]
	return(x) }

.fibR <- function(n, len=1) {
	retArr <- numeric(len)
	xa <- 0; xb <- 1;
	for(i in 0:n) {
		#init conds: fib(0)=0, fib(1)=1
		if (i <= 1) { xn <- i }
		#fib(n)=fib(n-1)+fib(n-2)
		else {
			xn <- xa+xb; xa <- xb; xb <- xn }
		## save results if iteration i is within the 
		## range from n-len to n
		j <- i - n + len;
		if (j>0) retArr[j] <- xn
	}
	return(retArr) }

.fibClosedForm <- function(n, len=1) {
	n <- (n-(len-1)):n
	phi <- (1+sqrt(5))/2
	return(round((phi^n - (1-phi)^n)/sqrt(5)))
}
#------------------------------------------calcFib


#calcGM---------------------------------2006-08-28
# Calculates the geometric mean of a vector of numbers
# Return the geometric mean of a vector of numbers
# Arguments:
#  x      - Vector of numbers
#  offset - Added value to validate zeroes
#  exzero - If TRUE, exclude zeroes
#----------------------------------------------JTS
calcGM <- function (x, offset = 0, exzero = TRUE) {
	x <- x[!is.na(x)]
	if (exzero) 
		x <- x[x > 0 & !is.na(x)]
	n <- length(x)
	if (n == 0) 
		return(0)
	x <- x + offset
	g <- exp(mean(log(x)))
	return(g)
}
#-------------------------------------------calcGM


#clearAll-------------------------------2012-12-06
# Remove all data in the global environment
# Arguments:
#  hidden  - if T remove all variables including dot variables
#  verbose - list all removed variables
#  PBSsave - if TRUE, do not remove .PBSmod
#----------------------------------------ACB/AE/RH
clearAll <- function(hidden=TRUE, verbose=TRUE, PBSsave=TRUE, pos=.PBSmodEnv) {
	objs <- ls(all.names = TRUE, pos = pos) #.GlobalEnv")
	if (verbose && length(objs))
		print(objs)
	rmlist <- ls(all.names = hidden, pos = pos) #.GlobalEnv")
	if(PBSsave)
		rmlist=rmlist[rmlist!=".PBSmod"]
	rm(list = rmlist, pos = pos) #.GlobalEnv")
	if (verbose) {
		cat("Removed:\n")
		if (length(rmlist))
			print(rmlist)
		else
			cat("\n")

		cat("Remaining:\n")
		if (length(ls(all.names = TRUE, pos = pos))) #.GlobalEnv")))
			print(ls(all.names = TRUE, pos = pos))    #.GlobalEnv"))
		else
			cat("\n")
	}
	invisible() }
#-----------------------------------------clearAll


#clearPBSext----------------------------2013-07-16
# Disassociate any number of file extensions from commands
#  previously save with setPBSext.
# Argument:
#  ext - optional character vector of file extensions to
#        clear; if unspecified, all associations are removed
#--------------------------------------------AE/RH
clearPBSext=function(ext){
  .initPBSoptions()
  tget(.PBSmod)
  if(missing(ext))
    .PBSmod$.options$openfile <- list() #packList("openfile",".PBSmod$.options",list()) #.PBSmod$.options$openfile<<-list()
  else{
    oldLen=length(.PBSmod$.options$openfile)
    #eval(parse(text=".PBSmod$.options$openfile <<- .removeFromList(.PBSmod$.options$openfile, ext)"))
    .PBSmod$.options$openfile <- .removeFromList(.PBSmod$.options$openfile, ext)
    if(oldLen!=length(.PBSmod$.options$openfile))
      #packList(".optionsChanged",".PBSmod$.options",TRUE) #.PBSmod$.options$.optionsChanged<<-TRUE
      .PBSmod$.options$.optionsChanged <- TRUE
  }
  tput(.PBSmod)
  invisible()
}
#--------------------------------------clearPBSext


#clearRcon------------------------------2009-05-14
# Clear the R console display.
#-----------------------------------------------NO
clearRcon <- function(os=.Platform$OS.type) {
	if (os!="windows") {
		err="'clearRcon' needs Windows OS to use Windows Scripting"
		cat(err,"\n"); return(invisible(err)) }
	tdir <- tempdir()
	fname <- paste(tdir, "\\clearRcon.vbs", sep="")
	cat('Dim pShell\n', file=fname)
	cat('Set pShell = CreateObject("WScript.Shell")\n', file=fname, append=TRUE)
	cat('pShell.AppActivate "R Console"\n', file=fname, append=TRUE)
	cat('pShell.SendKeys "^L"\n', file=fname, append=TRUE)
	cat('Set pShell = Nothing\n', file=fname, append=TRUE)
	system(paste("cscript //NoLogo", fname), minimized=TRUE)
	invisible(fname) }
#----------------------------------------clearRcon


#clipVector-----------------------------2007-06-09
# Clip a vector at one or both ends of a 
# specified value or string
#-----------------------------------------------RH
clipVector <- function (vec,clip,end=0) {
	#clip=as.character(substitute(clip))
	if (is.null(names(vec))) names(vec) <- 1:length(vec)
	if (any(end==c(0,1))) { # clip vector from the front
		znot = !vec%in%clip
		zzz  = cumsum(znot)
		z1   = match(1,zzz)
		vec  = vec[z1:length(vec)] }
	if (any(end==c(0,2))) { # clip vector from the end
		vrev = rev(vec)
		znot = !vrev%in%clip
		zzz  = cumsum(znot)
		z1   = match(1,zzz)
		vrev = vrev[z1:length(vrev)]
		vec  = rev(vrev) }
	return(vec)
}
#---------------------------------------clipVector


#compileDescription---------------------2006-08-28
#  Convert a GUI description file into a complete GUI desc List
#  which can be passed directly to createWin
# Arguments:
#  descFile - filename of GUI description file
#  outFile  - filename to save list to. 
#             WARNING: this will overwrite the file, if it currently exists
#----------------------------------------------ACB
compileDescription <- function(descFile, outFile="") {
	if (outFile!="")
		sink(outFile)
	x<-parseWinFile(descFile)
	cat(paste(
	"#This file was automaticaly generated from window description file \"",
	descFile, "\"\n#with the compileDescription function.\n\n", sep=""))
	cat("#This list can then be passed directly to createWin()\n\n")
	cat(paste("#To assign this list to a variable: GUIdesc<-eval(parse(\"", outFile, "\"))\n", sep=""))
	writeList(x)
	cat("\n")
	if (outFile!="")
		sink()
}
#-------------------------------compileDescription


#convSlashes----------------------------2009-02-16
# Convert unix "/" to R's "\\" if OS is windows.
#-----------------------------------------------RH
convSlashes=function(expr, os=.Platform$OS.type, addQuotes=FALSE){
	if (os=="windows") 
		expr=gsub("/","\\\\",expr)
	else 
		expr=gsub("\\\\","/",expr)
	if (addQuotes) expr=paste("\"",expr,"\"",sep="")
	return(expr) }
#--------------------------------------convSlashes


#createVector---------------------------2012-12-11
# Create a GUI with a vector widget and button
# Arguments:
#  vec:          a vector of widget variable names
#                if vec is named, then the names are used as widget variable 
#                names and the values are used as the default value
#  vectorLabels: if supplied, this vector of labels are printed above each entry box
#                There should be one label for every variable defined in vec
#                i.e. length(vectorLabels)==length(vec)
#  func:         function name as a string
#                If given, this function will be called whenever data is entered
#                i.e. Enter pressed, or submit button clicked. This user function
#                would then most likely use getWinVal()
#  windowname:   windowname to use for this GUI
# Output: If no user defined function is given (see func paramater), then global variables 
#         matching the variable name is set with the value of the widget
#         whenever text focus is in a widget and enter is pressed, or when submit is pushed.
#         Otherwise, func will be called and it is the user's responsibility to  make use of getWinVal
#-------------------------------------------ACB/RH
createVector <- function (vec, vectorLabels=NULL, func="", windowname="vectorwindow", env=NULL) {
	if(is.null(env)) env = parent.frame()
	if (is.null(names(vec))) {
		namesVal <- vec
		valuesVal <- ""
	}
	else {
		namesVal <- names(vec)
		valuesVal <- vec
	}
	if (!is.character(func))
		stop("func must be a character string")
	namesVal <- as.character(namesVal)
	if (is.null(vectorLabels)) 
		vecLabels <- names(vec)
	else {
		if (length(vectorLabels) != length(vec)) 
			stop("length of parameters vec and vectorLabels should be the same length")
		vectorLabels <- as.character(vectorLabels)
		vecLabels <- vectorLabels
	}
	winList <- list(list(title = "Vector", windowname = windowname, vertical = TRUE, 
		onclose = "", .widgets = list(list(type = "vector", names = namesVal, 
		length = length(vec), labels = vecLabels, values = valuesVal, 
		font = "", vertical = FALSE, "function" = func, enter = TRUE, 
		action = "", mode = "numeric", width = 6, sticky = "", 
		padx = 0, pady = 0), list(type = "button", "function" = func, 
		text = "Go", padx = 0, pady = 0)), .menus = list()))
	createWin(winList,env=env)
	invisible(winList)
}
#-------------------------------------createVector


#evalCall-------------------------------2014-03-04
# Evaluates a function call, resolving conflicting arguments.
#-----------------------------------------------RH
evalCall=function(fn, argu=list(), ..., envir=parent.frame(), checkdef=FALSE, checkpar=FALSE){
	fnam=as.character(substitute(fn))
	fnam.def=paste(fnam,"default",sep=".")

	base=formals(fnam)
	if (checkdef && exists(fnam.def)) {
		defs=formals(fnam.def)
		udefs=defs[setdiff(names(defs),names(base))]
		base=c(base,udefs) }
	base=base[setdiff(names(base),"...")]
	if (checkpar) {
		pars=par(); pars=pars[setdiff(names(pars),names(base))] # use only pars not in base
		forms=c(base,pars) } # all possible formal arguments
	else forms=base

	dots=list(...)
	# check to see if user has manipulated dots before passing to `evalCall`
	if (length(dots)==1 && !is.null(names(dots)) && names(dots)=="dots")
		dots = dots[["dots"]]

	argus=argu[setdiff(names(argu),names(dots))]
	given=c(argus,dots)
	allow=given[intersect(names(given),names(forms))]
	strargs=sapply(allow,deparse,width.cutoff=500,simplify=FALSE)
	strargs=sapply(strargs,paste,collapse="",simplify=FALSE) # collapse multiple strings
	argspec=character(0)

	for (i in names(strargs)) argspec=c(argspec,paste(i,"=",strargs[[i]]))
	expr=paste(fnam,"(",paste(argspec,collapse=","),")",sep="")
	eval(parse(text=expr)) 
	invisible(expr) }
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~evalCall


#findPat--------------------------------2006-08-28
# Searches all patterns in pat from vec, and returns the
# matched elements in vec.
# Arguments:
#  pat - character vector of patterns to match in vec.
#  vec - character vector where matches are sought.
#-----------------------------------------------RH
findPat <- function (pat, vec) {
	n <- length(vec)
	z <- NULL
	for (xstr in pat) {
		a <- regexpr(xstr, vec)
		b <- (1:n)[a > 0]
		z <- union(z, b) }
	found <- vec[z]
	return(found)
}
#------------------------------------------findPat


#findProgram----------------------------2006-08-28
# Given string name of a program - return the 
# complete path to program.  This was made without
# knowing Sys.which() existed - we may want to 
# deprecate this at some point.
#----------------------------------------------ACB
findProgram <- function( name, includename = FALSE )
{
	tmp <- Sys.which( name )
	if( includename == FALSE )
		tmp <- dirname( tmp )
	return( tmp )
}
#--------------------------------------findProgram


#focusRgui------------------------------2011-09-09
# Set focus to the RGui window.
#--------------------------------------------NO/RH
focusRgui = function (os = .Platform$OS.type) {
	if (os != "windows") {
		err = "'focusRgui' needs Windows OS to use Windows Scripting"
		cat(err, "\n")
		return(invisible(err))
	}
	tdir <- tempdir()
	fname <- paste(tdir, "\\focusRgui.vbs", sep = "")
	cmd <- "Set w = CreateObject(\"WScript.Shell\"): w.AppActivate(\"RGui\") : w.SendKeys(\"% x\") "
	cat(cmd, file = fname)
	system(paste("cscript //NoLogo", fname), minimized = TRUE)
	invisible(fname) 
}
#----------------------------------------focusRgui


#genMatrix------------------------------2006-08-28
#  Generate a test matrix for use in plotBubbles
# Arguments:
#  m     - number of rows
#  n     - number of columns
#  mu    - mean value of distribution
#  sigma - std deviation of distribution
#----------------------------------------------JTS
genMatrix <- function (m,n,mu=0,sigma=1) {
	matrix(rnorm(m*n,mean=mu,sd=sigma), m, n)
}
#----------------------------------------genMatrix


#getPBSext------------------------------2012-12-04
# Retrieve previously saved command.  Argument:
#  ext - file extension
#-------------------------------------------ACB/RH
getPBSext <- function(ext) {
	if (!exists(".PBSmod",envir=.PBSmodEnv))
		stop(".PBSmod was not found")
	tget(.PBSmod)
	if (missing(ext))
		return(.PBSmod$.options$openfile)
	if ( is.null( .PBSmod$.options$openfile[[ ext ]] ) )
		return(NULL)
	return(.PBSmod$.options$openfile[[ext]])
}
#----------------------------------------getPBSext


#getPBSoptions--------------------------2012-12-04
# Retrieve a user option.  Argument:
#   option - name of option to retrieve
#-------------------------------------------ACB/RH
getPBSoptions <- function(option) {
	tget(.PBSmod)
	if (missing(option))
		return(.PBSmod$.options)
	return(.PBSmod$.options[[option]])
}
#------------------------------------getPBSoptions


#isWhat---------------------------------2008-09-08
# Prints the class, mode, type, and attributes
# of the given object x.
#----------------------------------------------JTS
isWhat <- function(x) {
  cat("class: "); print(class(x));
  cat("mode:  "); print(mode(x));
  cat("type:  "); print(typeof(x));
  att=attributes(x)
  cat(paste("attributes:",ifelse(is.null(att)," NULL\n","\n"),sep=""))
  if (!is.null(att)) print(att)
  invisible() }
#-------------------------------------------isWhat


#openFile-------------------------------2013-04-05
# Opens a file for viewing based on System file
# extension association or .PBSmod$.options$openfile
#----------------------------------------ACB/RH/NB
openFile <- function(fname="", package=NULL)
{
	if (!exists(".PBSmod",envir=.PBSmodEnv)) .initPBSoptions()

	.openFile=function(fname, package)
	{
		if (fname=="")  fname=getWinAct()[1] # does this work?
		if(!is.null(package)) {
			# relative within package
			pkg_file <- system.file(fname, package=package)
			if( pkg_file == "" )
				stop(paste("File \"", fname,"\" does not exist in package \"", package, "\"", sep=""))
			fname <- pkg_file
		} else {
			# relative/absolute within file system
			fname <- normalizePath (fname);
		}
		# system.file and normalizePath both fail if the file
		# doesn't exist; if we've gotten this far, we
		# should be okay to open it
		# remove everything that precedes extension
		ext <- sub("^.*\\.", "", fname)
		tget(.PBSmod)
		if ( is.null( .PBSmod$.options$openfile[[ ext ]] ) ) {
			# nothing previously set with setPBSext
			if ( exists("shell.exec", mode="function" ) )  {
				# Windows
				shell.exec(fname)
				return(fname)
			} else if( file.exists( "/usr/bin/open" ) ) {
				# Mac OS X
				system( paste( "open", fname ) )
				return(fname)
			} else if( file.exists( "/usr/bin/xdg-open" ) ) {
				# Linux (xdg-open is a desktop-independent tool)
				system( paste( "xdg-open", fname ) )
				return(fname)
			} else if( file.exists( "/usr/bin/gnome-open" ) ) {
				# Linux (Gnome desktop)
				system( paste( "gnome-open", fname ) )
				return(fname)
			}
			stop(paste("There is no program associated with the ",
				"extension '", ext, "'\n",
				"Please set an association with the ",
				"setPBSext command\n"))
		} else {
			# matches extension previously set with setPBSext
			cmd <- getPBSext(ext)
			# RH: gsub needs to expand WinOS "\\" delimiters introduced by `normalizePath` above.
			cmd <- gsub("%f", gsub("\\\\","\\\\\\\\",fname), cmd)
			if (.Platform$OS.type=="windows")
				shell(cmd,wait=FALSE)
			else
				system(cmd,wait=FALSE)
			return(cmd)
		}
	}
	ops=sapply(fname,.openFile, package=package)
	invisible(ops) }
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~-openFile


#openUG---------------------------------2009-12-04
# Open package User Guide 'pkg-UG.pdf' if it exists.
# Essentially a wrapper for 'openFile'.
#-----------------------------------------------RH
openUG = function(pkg="PBSmodelling"){
	pkgnam = as.character(substitute(pkg))
	openFile(paste("/doc/",pkgnam,"-UG.pdf",sep=""),pkgnam)
}
#-------------------------------------------openUG


#pad0-----------------------------------2011-12-12
# Takes numbers (or text coerced to numeric), 
#  converts them to integers then text, and 
#  pads them with leading zeroes.
# Arguments:
#    x - Vector of numbers
#    n - Length of padded integer
#    f - Factor of 10 to expand x by
# Note: For meaningful results, n should be at least as
#       large as the order of factored x (x * 10^f).
#-----------------------------------------------RH
pad0 <- function (x, n, f = 0) {
	xin <- x
	fin <- f
	flist <- as.list(f); names(flist) = f
	xnum <- suppressWarnings(as.numeric(xin))
	zstr <- is.na(xnum); znum = !zstr
	for (f in fin) {
		xout <- xin
		if (any(znum)) {
			x <- xnum[znum]
			xord <- max(ceiling(log10(abs(x * 10^f))), na.rm = TRUE);
			if (any(max(abs(x * 10^f)) == 10^(-10:10))) xord <- xord + 1;
			if (n < xord) n <- xord  # No padding occurs if n<=xord
			x <- round(x, f) * 10^f
			xneg <- x < 0;
			x <- abs(x);  x <- format(x, scientific=FALSE)
			xout[znum] = x }
		if (any(zstr)) 
			xout[zstr] <- paste(xout[zstr],paste(rep(0,f),collapse=""),sep="")
		xout <- gsub(" ", "", xout); nx <- length(xout);
		base0 <- rep(paste(rep(0, n), collapse = ""), nx);
		nchr <- nchar(xout); ndiff <- n - nchr;
		add0 <- substring(base0, 1, ndiff);
		xout <- paste(add0, xout, sep = "");
		if (any(znum))
			xout[znum][xneg] <- paste("-", xout[znum][xneg], sep = "");
		attr(xout, "input") <- xin
		flist[[as.character(f)]] <- xout
	}
	if (length(fin)==1)
		padout = flist[[as.character(fin)]]
	else if (length(xin)==1) {
		padout = sapply(flist,function(x){x[1]})
		attr(padout, "input") <- xin }
	else padout <- flist
	return(padout) }
#---------------------------------------------pad0


#pause----------------------------------2006-08-28
# Pause, typically between graphics displays
# Arguments:
#  s  - string to display to user
#----------------------------------------------JTS
pause <- function (s = "Press <Enter> to continue") {
	cat(s)
	readline()
	invisible()
}
#--------------------------------------------pause


#readPBSoptions-------------------------2012-12-04
# Load PBS options from a text file. The loaded options will
# overwrite existing ones in memory; however, an existing
# option in memory will not be cleared if this option does
# not exist in the options file.
# Input:
#  fname - name of options file (or path to this file)
# Output:
#  returns FALSE if file did not exist or if read failed
#  otherwise returns TRUE
#-------------------------------------------ACB/RH
readPBSoptions=function(fname="PBSoptions.txt"){
	.initPBSoptions()
	optList=try(readList(fname), silent=TRUE)
	if(class(optList)=="try-error")
		return(FALSE)
	tget(.PBSmod)
	#eval(parse(text=".PBSmod$.options <<- .mergeLists(.PBSmod$.options, optList)"))
	.PBSmod$.options <- .mergeLists(.PBSmod$.options, optList)
	if(fname!="PBSoptions.txt")
		#packList(".optionsFile",".PBSmod$.options",fname) #.PBSmod$.options$.optionsFile<<-fname
		.PBSmod$.options$.optionsFile <- fname
	#packList(".optionsChanged",".PBSmod$.options",NULL) #.PBSmod$.options$.optionsChanged<<-NULL
	.PBSmod$.options$.optionsChanged <- NULL
	tput(.PBSmod)
}
#-----------------------------------readPBSoptions


#runDemos-------------------------------2013-06-27
# Display a GUI to display something equivalent to R's demo()
#-------------------------------------------ACB/RH
runDemos <- function (package) {
	#if (!exists(".dwd",where=1)) assign(".dwd",getwd(),envir=.GlobalEnv)
	#if (!exists(".dls",where=1)) assign(".dls",c(".dls",ls(pos = 1, all.names=TRUE)),envir=.GlobalEnv)
	if (!exists(".dwd",envir=.PBSmodEnv)) assign(".dwd",getwd(),envir=.PBSmodEnv) #.GlobalEnv)
	if (!exists(".dls",envir=.PBSmodEnv)) assign(".dls",c(".dls",ls(pos = .PBSmodEnv, all.names=TRUE)),envir=.PBSmodEnv) #.GlobalEnv)
	try(closeWin(),silent=TRUE)
	x <- demo(package = .packages(all.available = TRUE))
	if (missing(package)) {
		#display a list of packages to choose from
		pkgDemo <- unique(x$results[,"Package"])
		wintext <- c( "window title=\"R Demos\" name=pbsdemo onclose=.dClose",
			"label text=\"R Demos                                                  \" font=\"bold underline\" fg=blue padx=10 sticky=w",
			"label text=\"Select a package to view available demos:\" sticky=W padx=12 font=\"bold 10\"" )

		#create droplist labels with counts
		pkg_labels <- c()
		for(pkg in pkgDemo) {
			len <- length(x$results[,"Package"][x$results[,"Package"]==pkg])
			if (len==1)
				items <- "(1 demo)"
			else
				items <- paste("(",len," demos)", sep="")
			pkg_labels <- c( pkg_labels, paste( pkg, items ) )
		}

		wintext[ length(wintext) + 1 ] <- paste( "droplist name=pkg values=", .addslashes( pkgDemo, asvector=TRUE ), " labels=", .addslashes( pkg_labels, asvector = TRUE ), " add=FALSE mode=character", sep="" )

		wintext[ length(wintext) + 1 ] <- "button function=.viewPkgDemo action=pkg text=\"View Demos\" sticky=w padx=12 bg=aliceblue"

		createWin(wintext, astext=TRUE, env = parent.env( environment() ) )
		return(invisible(NULL))
	}
	#display demos from a certain package
	eval(parse(text=paste("OK=require(",package,",quietly=TRUE)",sep="")))
	if (!OK) {
		mess=paste(package,"package is not available")
		showAlert(mess); stop(mess) }
	x <- x$results[x$results[,"Package"]==package,]

	# It is critical that the label widget named package *only* contain the package name + whitespace -> it is used by getWinVal() later
	wintext <- c( paste( "window title=\"R Demos:", package, "\" name=pbsdemo onclose=.dClose", sep="" ),
		paste( "label name=package text=", .addslashes( package ), " font=\"bold underline\" fg=red3 padx=10 sticky=w", sep="" ),
		"label text=\"Select a demo to view:\" sticky=W padx=12 font=\"bold 10\"" )

	if (is.null(dim(x))) {
		tmp<-names(x)
		dim(x)<-c(1,4)
		colnames(x)<-tmp
	}

	droplist_data <- c()
	for(j in 1:length(x[,1])) {
		demoDir <- file.path(x[j,"LibPath"], package, "demo")
		#====== DEPRECATED functions in `tools` =======
		#path <- tools::list_files_with_type(demoDir, "demo")
		#path <- path[x[j,"Item"]==tools::file_path_sans_ext(basename(path))]
		#==============================================
		dext = paste(c("R", "r"),"$",sep="")
		dfile = findPat(dext,list.files(path=demoDir)) #list.files(path=demoDir,pattern=dext)
		dpref = sub("([^.]+)\\.[[:alnum:]]+$", "\\1", dfile)
		disin = is.element(x[,"Item"],dpref)
		path  = paste(demoDir,dfile[disin],sep="/")
#browser();return()
		droplist_data[ j ] <- path [j]
	}
	titles <- x[,"Title"]
	title_cut_off <- 50 #cut off titles longer than this
	labels <- paste( x[,"Item"], " ::: ", substring( titles, 1, title_cut_off ), ifelse( nchar(titles) > title_cut_off, "...", "" ), sep="" )

	wintext[ length( wintext ) + 1 ] <- paste( "droplist name=demo values=", .addslashes( droplist_data, asvector=TRUE ), " labels=", .addslashes( labels, asvector=TRUE ), " add=FALSE mode=character sticky=W padx=20 width=55 function=.dUpdateDesc", sep="" )
	wintext[ length( wintext ) + 1 ] <- paste( "label name=demo_desc text=", .addslashes( x[1,"Title"] ), " sticky=W wraplength=500 padx=20 pady=\"20 0\"", sep="" )

	wintext[ length( wintext ) + 1 ] <- "grid 1 3 sticky=w pady=3"
	wintext[ length( wintext ) + 1 ] <- "button function=.viewPkgDemo action=demo text=\"Run Demo\" sticky=w padx=12 bg=greenyellow"
	wintext[ length( wintext ) + 1 ] <- "button function=.viewPkgDemo action=source text=\"View Source\" sticky=w padx=12"
	wintext[ length( wintext ) + 1 ] <- "button function=runDemos action=\"\" text=\"All Packages\" sticky=w padx=12"

	createWin( wintext, astext=TRUE, env=parent.env( environment() ) )
	return(invisible(NULL))
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~runDemos

#.viewPkgDemo---------------------------2009-03-04
# Display a GUI to display something equivalent to R's demo()
#-------------------------------------------ACB/RH
.viewPkgDemo <- function() {
	act <- getWinAct()[1]
	if (act=="pkg") {
		package=getWinVal("pkg")$pkg
		eval(parse(text=paste("OK=require(",package,",quietly=TRUE)",sep="")))
		if (!OK) {
			mess=paste(package,"package is not available")
			showAlert(mess); stop(mess) }
		return(runDemos(package)) }
	if (act=="demo") {
		demo <- getWinVal("demo")$demo
		source(demo, echo=TRUE, max.deparse.length=100)
		return(invisible(NULL))
	}
	if (act=="source") {
		demo <- getWinVal("demo")$demo
		openFile(demo)
		return(invisible(NULL))
	}
}
#.dClose--------------------------------2013-07-02
# Function to execute on closing runDemos().
#----------------------------------------------ACB
.dClose <- function() {
	act <- getWinAct()[1];
	closeWin();
	setwd(tcall(.dwd))
	#if (is.null(act) || act=="demo") {
		remove(list = setdiff(ls(pos=.PBSmodEnv, all.names=TRUE), tcall(.dls)), pos=.PBSmodEnv)
		#remove(list = c(".dwd", ".dls"), pos=.PBSmodEnv)  # final good-bye
	#}
	return()
}
#.dUpdateDesc---------------------------2009-03-04
# Updates demo window with demo descriptions.
#----------------------------------------------ACB
.dUpdateDesc <- function() {
	vals <- getWinVal()
	demo.id <- vals$demo.id
	package <- .trimWhiteSpace( vals$package )
	x <- demo(package = .packages(all.available = TRUE))
	x <- x$results[x$results[,"Package"]==package,]
	if (is.null(dim(x))) {
		tmp<-names(x)
		dim(x)<-c(1,4)
		colnames(x)<-tmp
	}
	setWinVal( list( demo_desc=x[demo.id,"Title"] ) )
}
#=====================================END runDemos


#runExample-----------------------------2012-12-10
# Display a single GUI example.
#-----------------------------------------------RH
runExample <- function (ex, pkg="PBSmodelling") {

	.runExHelperQuit <- function() {
		#closeWin(name=allWin)
		setwd(tget(.cwd))
		remove(list = setdiff(ls(pos = .PBSmodEnv), tget(.cls)), pos = .PBSmodEnv)
		junk=file.remove(tget(.cfl)) # assign only to suppress printing T/F
		return() }
	assign(".cls",ls(pos=.PBSmodEnv, all.names=TRUE),envir=.PBSmodEnv)    # current list in .PBSmodEnv
	assign(".cwd",getwd(),envir=.PBSmodEnv)                               # current working system directory
	assign(".cfl",list.files(tempdir(),full.names=TRUE),envir=.PBSmodEnv) # current file list in temporary system directory
	assign(".runExHelperQuit",.runExHelperQuit,envir=.PBSmodEnv) #.GlobalEnv)

	pdir <- system.file(package = pkg)                # package directory
	edir <- paste(pdir, "/", "examples", sep = "")    # examples directory
	fnam <- paste(edir, list.files(edir), sep = "/")  # file names in examples directory
	bnam <- basename(fnam)                            # basenames

	if (missing(ex) || !any(is.element(paste(ex,"Win.txt",sep=""),setdiff(bnam,"runExamplesWin.txt")))) {
		mess = paste("Your example does not exist.\nChoose from:\n     ",
			paste(setdiff(findPrefix("Win.txt",edir),"runExamples"),collapse="\n     "),sep="")
		showAlert(mess)
		stop ("Choose another example") }
	rtmp <- tempdir()                                 # R's temporary directory
	file.copy(fnam, rtmp, overwrite = TRUE)
	setwd(rtmp)
	wnam <- paste(ex,"Win.txt",sep="")  # window description file or XML talk file
	wtxt <- readLines(wnam)
	wtxt[1] = paste(wtxt[1], "onClose=.win.runExHelperQuit")
	writeLines(wtxt,wnam)
	rnam <- paste(ex,".r",sep="")                     # R code file
	if (is.element(rnam,bnam)) source(rnam,local=.PBSmodEnv) # seems to see the function in .PBSmodEnv
	invisible() }
#---------------------------------------runExample


#runExamples----------------------------2012-12-17
# Display a master GUI to display examples
#-------------------------------------------RH/ACB
runExamples <- function () {
	allWin=c("runE","window","widWin","testW","choisir","presentwin") # all potential windows open
	.runExHelper <- function() {
		getWinVal(scope = "L")
		act <- getWinAct()[1]
		if (!exists("act") || !exists("eN")) return()
		if (act == "quit") { 
			tget(.runExHelperQuit)()
		} else if (act=="clear") {
			wtxt <- "No examples chosen"
			closeWin(name=setdiff(allWin,"runE"))
		} else if (act == "__USE_EDIT__") {
			if (wtxt == "\n" || wtxt == "No examples chosen\n") return()
			winDesc <- strsplit(wtxt, "\n")[[1]]
			createWin(winDesc, astext = TRUE)
			return()
		} else if (act=="swisstalk") {
			#TODO refactor this - should be more generic
			closeWin(name=setdiff(allWin,"runE"))
			presentTalk( "swisstalk.xml" )
		} else {
			if (act!="TestFuns")
				closeWin(name=setdiff(allWin,c("runE","window")))
			f <- paste(act, ".r", sep = "" )
			#assert that files match case
			stopifnot( any( dir() == f ) == TRUE )
			source( f ,local=.PBSmodEnv )
			wnam <- paste(act, "Win.txt", sep = "") # window description file
			wtxt <- paste(readLines(wnam), collapse = "\n")
		}
		setWinVal(list(wtxt=wtxt), winName="runE")
	}
	.runExHelperQuit <- function() {
		closeWin(name=allWin)
		setwd(tget(.cwd))
		remove(list = setdiff(ls(pos = .PBSmodEnv), tget(.cls)), pos = .PBSmodEnv)
		junk=file.remove(tget(.cfl)) # assign only to suppress printing T/F
		return()
	}
	assign(".runExHelper",.runExHelper,envir=.PBSmodEnv) #.GlobalEnv)
	assign(".runExHelperQuit",.runExHelperQuit,envir=.PBSmodEnv) #.GlobalEnv)
	assign(".cls",ls(pos=.PBSmodEnv, all.names=TRUE),envir=.PBSmodEnv)    # current list in .PBSmodEnv
	assign(".cwd",getwd(),envir=.PBSmodEnv)                               # current working system directory
	assign(".cfl",list.files(tempdir(),full.names=TRUE),envir=.PBSmodEnv) # current file list in temporary system directory
	pckg <- "PBSmodelling"
	pdir <- system.file(package = pckg)               # package directory
	edir <- paste(pdir, "/", "examples", sep = "")    # examples directory
	fnam <- paste(edir, list.files(edir), sep = "/")  # file names in examples directory
	rtmp <- tempdir()                                 # R's temporary directory
	file.copy(fnam, rtmp, overwrite = TRUE)
	setwd(rtmp)
	createWin("runExamplesWin.txt")
	msg <- paste("Examples are running in ", rtmp, sep = "")
	setWinVal(list(wtxt = msg), winName = "runE") 
	invisible() }
#--------------------------------------runExamples


#selectDir------------------------------2006-06-06
# Prompts user to select a directory - and returns it.
#----------------------------------------------ACB
selectDir <- function( initialdir = getwd(), mustexist = TRUE, title = "", usewidget = NULL )
{
	#get val from widget
	if( is.null( usewidget ) == FALSE )
		initialdir <- getWinVal()[[ usewidget ]]

	d <- tclvalue( tkchooseDirectory( initialdir = initialdir, mustexist = mustexist, title = title ) )

	#set val to widget
	if( is.null( usewidget ) == FALSE && d != "" ) {
		val <- list()
		val[[ usewidget ]] <- d
		setWinVal( val )
	}
		
	return( d )
}
#----------------------------------------selectDir


#selectFile-----------------------------2006-06-06
# Prompts a user to open/save a file(s).
#----------------------------------------------ACB
selectFile <- function(
	initialfile="", 
	initialdir=getwd(), 
	filetype=list(c("*", "All Files")), 
	mode="open", 
	multiple = FALSE, 
	title="",
	defaultextension="",
	usewidget = NULL
	)
{
	mode <- tolower( mode )
	stopifnot( mode == "open" || mode == "save" )

	if( is.null( usewidget ) == FALSE ) {
		#value in widget should be used as file
		initialfile <- getWinVal()[[ usewidget ]]
		if( multiple == TRUE )
			stop( "multiple and usewidget can not be used together" )
	}
	
	#prepare filetypes to correct format for tk
	filetypes <- ""
	for(i in 1:length(filetype)) {
		filetype[[i]]
		if (is.na(filetype[[i]][2]))
			filetype[[i]][2] <- filetype[[i]][1]
		if (filetype[[i]][1] != "*" && substr(filetype[[i]][1],1,1)!=".")
			filetype[[i]][1] <- paste(".", filetype[[i]][1], sep="")
		filetypes <- paste(filetypes, " {{", filetype[[i]][2], "} {", filetype[[i]][1], "}}", sep="")
	}
	filetypes <- .trimWhiteSpace(filetypes)
	
	#prepare common args
	args <- list(
		initialdir=initialdir,
		initialfile=initialfile,
		filetypes=filetypes,
		title=title,
		defaultextension=defaultextension
	)
	
	if( mode == "open" ) {
		args[[ "multiple" ]] <- multiple
		files <- tclvalue(do.call( tkgetOpenFile, args ))
	} else if( mode == "save" ) {
		if( multiple == TRUE ) stop( "multiple=TRUE is not supported by \"save\" mode" )
		files <- tclvalue(tkgetSaveFile( initialdir=initialdir, initialfile=initialfile, filetypes=filetypes ))
	} else {
		stop("mode not supported")
	}
	
	#split up string
	if( files == "" )
		return( c() )

	if( is.null( usewidget ) == FALSE ) {
		#store value in widget
		stopifnot( multiple == FALSE )
		val <- list()
		val[[ usewidget ]] <- files
		setWinVal( val )
		return( files )
	}
	if( multiple == FALSE )
		return( files )
	return( .tclArrayToVector( files ) )
}
#---------------------------------------selectFile


#setPBSext------------------------------2012-12-04
# Associate a new command with file types;
#  use "%f" in cmd to designate where the filename will be placed.
# AE: Added the setting of 'optionsChanged'
# Arguments:
#  ext - file extension
#  cmd - cmd to open these types of files
#----------------------------------------ACB/AE/RH
setPBSext <- function(ext, cmd) {
	if (!exists(".PBSmod",envir=.PBSmodEnv)) 
		stop(".PBSmod was not found")
	tget(.PBSmod)
	if (!any(grep("%f", cmd)))
		stop(paste("No %f was found in supplied command \"", cmd, 
		           "\".\n%f must be used to indicate where the filename will ",
		           "be inserted by openfile().\n",
		           "Did you mean \"", cmd, " %f\"?", sep=""))
		           
	if(is.null(.PBSmod$.options$openfile[[ext]]) || .PBSmod$.options$openfile[[ext]]!=cmd)
		#packList(".optionsChanged",".PBSmod$.options",TRUE) #.PBSmod$.options.optionsChanged<<-TRUE
		.PBSmod$.options.optionsChanged <- TRUE
	#eval(parse(text=".PBSmod$.options$openfile[[ext]] <<- cmd"))
	.PBSmod$.options$openfile[[ext]] <- cmd
	tput(.PBSmod)
}
#----------------------------------------setPBSext


#setPBSoptions--------------------------2012-12-04
# Change user options. Arguments:
#   option - name of option to change
#   value  - new value of option
# AE: Now a value '.PBSmod$.options$.optionsChanged' is set to TRUE when an option is changed,
#   so that the user doesn't always have to be prompted to save the options file.
#   By default, '.PBSmod$.options$.optionsChanged' is not set or NULL.
#   Also, if an option is set to "" or NULL then it is removed.
#   '.initPBSoptions()' is now called first (options starting with a dot "." do not set '.optionsChanged').
# RH: if the value is a sublist of an option, it can be changed individually using 'sublist=TRUE'.
#----------------------------------------ACB/AE/RH
setPBSoptions <- function(option, value, sublist=FALSE) {
	.initPBSoptions()
	if(!is.null(value) && length(value)==1 && value=="") value=NULL
	tget(.PBSmod)
	if(substr(option, 1, 1)!="." && !identical(.PBSmod$.options[[option]], value))
		#packList(".optionsChanged",".PBSmod$.options",TRUE) #.PBSmod$.options$.optionsChanged<<-TRUE
		.PBSmod$.options$.optionsChanged <- TRUE
	if(is.null(value) && !sublist)
		#eval(parse(text=".PBSmod$.options <<- .removeFromList(.PBSmod$.options, option)"))
		.PBSmod$.options <- .removeFromList(.PBSmod$.options, option)
	else{
		if(is.list(value) && sublist){
			for (i in 1:length(value)){
				ii=names(value[i]); if (ii=="") next
				ival=value[[i]]
				#txt=paste(".PBSmod$.options$",option,ifelse(ii=="","","$"),ii," <<- ival",sep="")
				txt=paste(".PBSmod$.options$",option,ifelse(ii=="","","$"),ii," <- ival",sep="")
				eval(parse(text=txt))
			}
		}
		else
			#eval(parse(text=".PBSmod$.options[[option]] <<- value"))
			.PBSmod$.options[[option]] <- value
	}
	tput(.PBSmod)
}
#------------------------------------setPBSoptions


#show0----------------------------------2011-09-08
# Shows decimal places including zeroes (string)
# Return character representation of number with
# specified decimal places.
# Arguments:
#  x  - Number as scalar or vector
#  n  - Number of decimal places to show, include zeroes
#  add2int - If TRUE, add zeroes on the end of integers
#  round2n - if TRUE, round x first to n decimal places
#-----------------------------------------------RH
show0 <- function (x, n, add2int=FALSE, round2n=FALSE) {
	if (round2n) x = round(x,n)
	x <- as.character(x)
	oldx <- x
	pnt <- regexpr("\\.", x)
	z <- grep(-1, pnt)
	x[z] <- paste(x[z], ".", sep = "")
	pnt[z] <- nchar(x)[z]
	int <- substring(x, 1, pnt)
	end <- substring(x, pnt + 1)
	nx <- length(end)
	base0 <- rep(paste(rep(0, n), collapse = ""), nx)
	nchr <- nchar(end)
	ndiff <- n - nchr
	add0 <- substring(base0, 1, ndiff)
	newx <- paste(int, end, add0, sep = "")
	if (!add2int) 
		newx[z] <- oldx[z]
	return(newx) }
#--------------------------------------------show0


#showArgs-------------------------------2009-02-23
#  show arguments of a widget definition
# Arguments:
#  widget - only show information about supplied widget
#-------------------------------------------ACB/RH
showArgs <- function(widget, width=70, showargs=FALSE) {
	x <- .widgetDefs
	if (missing(widget)) widget=sort(names(x))
	x=x[widget]; xnam=names(x)
	out=character(0) # output file

	for(i in 1:length(x)) {
		#print widget name, and underline it
		cat(xnam[i])
		cat(paste("\n",paste(rep("-",nchar(xnam[i])),collapse=""),"\n",sep=""))
		
		expr=xnam[i]
		for(j in 2:length(x[[i]])) { 
			argu=x[[i]][[j]]$param
			if (x[[i]][[j]]$required==TRUE) { }
			else if ( any( names( x[[i]][[j]] ) == "default" ) ) {
				if (x[[i]][[j]]$class=="character" || x[[i]][[j]]$class=="characterVector")
					delim="\""
				else
					delim=""				
				default <- ifelse( is.null( x[[i]][[j]]$default ), "NULL", paste(delim,x[[i]][[j]]$default,delim,sep="") )
				argu=c(argu,"=")
				argu=c(argu,default)
			}
			else {
				cat("\n\n")
				stop(paste(xnam[i],"::",x[[i]][[j]]$param, "is not required, but has no default."))
			}
			expr=paste(expr,paste(argu,collapse=""),sep=" ")
		}
		wexp=strwrap(expr,width=width,exdent=3)
		cat(paste(wexp,collapse="\n")); cat("\n\n")
		out=c(out,expr) # save for output

		if (showargs) {
			for(j in 1:length(x[[i]])) { 
				cat(x[[i]][[j]]$param)
				if (x[[i]][[j]]$required==TRUE) {
					cat("\t"); cat("(required)") }
				cat("\n") }
			cat("\n\n") }
	}
	invisible(out) }
#-----------------------------------------showArgs


#showHelp-------------------------------2015-01-20
# Show HTML help files for package contents.
#-----------------------------------------------RH
showHelp <- function(pattern="methods", ...) {
	oldopts = options(); on.exit(options(oldopts))
	options(stringsAsFactors=FALSE)
	options(warn = -1)
	Apacks = .packages(all.available = TRUE) # all packages
	Spacks = findPat(pattern,Apacks)         # show packages that match the pattern
	npacks = length(Spacks)
	if (npacks==0) { print("No such package"); return() }
	unpackList(list(...))
	# `home' code from `help.start'
	home <- 
	if (!is.element("remote",ls(envir=sys.frame(sys.nframe()))) || is.null(remote)) {
		if (as.numeric(R.version$svn) >= 67548)
			httpdPort = tools::startDynamicHelp(NA)
		else {
			showAlert("This function only available for R (>= 3.2.0)")
			stop("This function only available for R (>= 3.2.0)", call. = FALSE)
		}
		if (httpdPort > 0L) {
			if ("update" %in% ls(envir=sys.frame(sys.nframe())) && update)
				make.packages.html(temp = TRUE)
			paste0("http://127.0.0.1:", httpdPort) 
		}
		else stop("showHelp() requires the HTTP server to be running", call. = FALSE)
	}
	else remote
	getURL = function(x) {
		url=paste(home,"/library/",x,"/html/00Index.html",sep="")
		browseURL(url)
		return(url) }
	URLs=sapply(Spacks,getURL)
	invisible(list(Apacks=Apacks,Spacks=Spacks,URLs=URLs))
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~showHelp


#showPacks------------------------------2009-02-17
# Show packages that need to be installed.
#-----------------------------------------------RH
showPacks = function(packs=c("PBSmodelling","PBSmapping","PBSddesolve",
            "rgl","deSolve","akima","deldir","sp","maptools","KernSmooth")) {
	warn <- options()$warn
	options(warn = -1)
	Apacks = .packages(all.available = TRUE)   # all packages
	Ipacks = sort(findPat(packs,Apacks))       # installed packages that match those required
	Mpacks = sort(setdiff(packs,Ipacks))       # missing packages
	if (length(Mpacks)==0)
		showAlert("Congratulations! All specified pakages installed.","Package Check","info")
	else {
		mess=paste("You need to install these packages:\n     ",paste(Mpacks,collapse="\n     "),sep="")
		showAlert(mess,"Package Check","warning")
	}
	options(warn = warn)
	invisible(list(Apacks=Apacks,Ipacks=Ipacks,Mpacks=Mpacks)) }
#----------------------------------------showPacks


#showRes--------------------------------2008-04-29
# Show results of the calculation in string x
#----------------------------------------------JTS
showRes <- function(x, cr=TRUE, pau=TRUE) {
  cat(">",x," "); if(cr==TRUE) cat("\n");
  if(pau) pause("...") else cat("...\n");
  xres <- eval(parse(text=x));
  print(xres); 
  cat("************** Results shown above **************\n\n");
  if(pau) pause();
  invisible(xres); };
#------------------------------------------showRes


#showVignettes--------------------------2013-06-27
# Display a GUI to display something equivalent to R's vignette()
# Arguments: package = string specifying a package name.
#--------------------------------------------AE/RH
showVignettes = function(package){
	if (!exists(".dwd",envir=.PBSmodEnv)) assign(".dwd",getwd(),envir=.PBSmodEnv)
	if (!exists(".dls",envir=.PBSmodEnv)) assign(".dls",c(".dls",ls(pos = .PBSmodEnv, all.names=TRUE)),envir=.PBSmodEnv)
	closeWin();
	x <- vignette()
	xres = as.data.frame(x[["results"]])
	if (missing(package)) {
		nvig = sapply(split(xres$Item,xres$Package),length)
		tvig = paste(names(nvig)," (",nvig," vignette",ifelse(nvig>1,"s",""),")",sep="") # radio button labels
		rvig = paste("radio name=pkg value=",names(nvig),
			" text=\"",tvig,"\" mode=character sticky=W padx=12",sep="")
		win = c("window title=\"R Vignettes\" name=pbsvignette onclose=\".dClose\"",
			"grid 3 1",
			"label text=\"Select a package to view\\navailable vignettes:\" sticky=W padx=12 font=\"bold 11\"",
			paste("grid ",length(rvig)," 1",sep=""),
			rvig,
			"button function=.viewPkgVignette action=pkg text=\"View Vignettes\" sticky=E padx=12 bg=greenyellow")
		assign("xxy",win,envir=.PBSmodEnv)
		createWin(win,astext=TRUE)
		return(invisible(NULL))
	}
	#display vignettes from a certain package
	xpac <- xres[is.element(xres$Package,package),]
	if (nrow(xpac)==0) stop("No such package on your system")
	if (is.null(dim(xpac))) {
		tmp <- names(xres); dim(xpac) <- c(1,dim(xres)[2]); colnames(xpac) <- tmp }
	vdir <- file.path(xpac[1,"LibPath"], package, "doc") # always the same for a single package
	#====== DEPRECATED functions in `tools` =======
	#path <- tools::list_files_with_type(vdir, "vignette")
	#path <- path[xpac[,"Item"]==tools::file_path_sans_ext(basename(path))]
	#==============================================
	vext = paste(c(outer(c("R", "r", "S", "s"), c("nw", "tex"), paste, sep = ""), "Rmd"),"$",sep="")
	vfile = findPat(vext,list.files(path=vdir)) #list.files(path=vdir,pattern=vext)
	vpref = sub("([^.]+)\\.[[:alnum:]]+$", "\\1", vfile)
	visin = is.element(xpac[,"Item"],vpref)
	path  = paste(vdir,vfile[visin],sep="/")
#browser();return()
	nvig = sapply(split(xpac$Item,xpac$Package),length)
	rvig = paste("radio name=vignette value=\"",path,"\" text=\"",xpac[,"Item"],"\" mode=character font=\"underline 10\" sticky=W padx=12",sep="")
	lvig = paste("label text=\"",xpac[,"Title"],"\" sticky=W  wraplength=500 padx=20",sep="")
	win = c("window title=\"R Vignettes\" name=pbsvignette onclose=\".dClose\"",
		"grid 3 1",
		"label text=\"Select a vignette:\" sticky=W padx=12 font=\"bold 11\"",
		paste("grid ",length(rvig)," 2 byrow=FALSE",sep=""),
		rvig,lvig,
		"grid 1 3",
		"button function=.viewPkgVignette action=vignette text=\"View Vignette\" sticky=W padx=12 bg=lightsteelblue1",
		"button function=.viewPkgVignette action=source text=\"View Source\" sticky=W padx=12 bg=lightsteelblue1",
		"button function=showVignettes action=\"\" text=\"All Packages\" sticky=W padx=12 bg=greenyellow")
	assign("xx",win,envir=.PBSmodEnv)
	createWin(win,astext=TRUE)
	return(invisible(NULL))
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~showVignettes

#.viewPkgVignettes----------------------2013-07-02
# Display a GUI to display something equivalent to R's vignette()
#--------------------------------------------AE/RH
.viewPkgVignette <- function() {
	act <- getWinAct()[1]
	if (act=="pkg")
		return(showVignettes(getWinVal("pkg")$pkg))
	vignette <- getWinVal("vignette")$vignette
	if (act=="vignette")
		vignette <- paste(sub("([^.]+)\\.[[:alnum:]]+$", "\\1", vignette), ".pdf", sep="")
		#vignette <- paste(tools::file_path_sans_ext(vignette), ".pdf", sep="")  ### DEPRECATED
	openFile(vignette)
	return(invisible(NULL))
}
#================================END showVignettes


#testWidgets----------------------------2006-08-28
# Display a "master" GUI that displays other sample GUIs
#-----------------------------------------------RH
testWidgets <- function () {
	.testWidHelper <- function() {
		getWinVal(scope="L");
		if (getWinAct()[1]=="__USE_EDIT__") {
			if (wtxt=="\n" || wtxt=="No widgets displayed\n")
				return()
			winDesc <- strsplit(wtxt, "\n")[[1]]
			createWin(winDesc, astext=TRUE)
			return()
		}
		if (wN==0) {
			wtxt <- "No widgets displayed";
			closeWin(name="widWin");
		}
		else {
			pckg <- "PBSmodelling"; dnam <- "testWidgets";
			act  <- getWinAct()[1];

			#get window descript
			wtmp <- paste(dnam,"/",act,".txt",sep="");
			wnam <- system.file(wtmp,package=pckg)
			wtxt <- paste(readLines(wnam),collapse="\n");

			#source r file if appropriate (can I source it under a different environment to protect myself?)
			rtmp <- paste(dnam,"/",act,".r",sep="");
			rnam <- system.file(rtmp,package=pckg)

			if( file.exists( rnam ) )
				source( rnam )

			createWin(wnam);

		}
		setWinVal(list(wtxt=wtxt), winName="testW");
	}
	assign(".testWidHelper",.testWidHelper,envir=.PBSmodEnv) #.GlobalEnv)
	pckg <- "PBSmodelling"; dnam <- "testWidgets";
	wtmp <- paste(dnam,"/","testWidgetWin.txt",sep="");
	wnam <- system.file(wtmp,package=pckg)
	createWin(wnam);
}
#--------------------------------------testWidgets


#tget/tcall/trpint/tput-----------------2012-12-06
# Functions to get, put, and print objects into the .PBSmodEnv.
# CRAN packages can no longer modify user's working environment.
#-----------------------------------------------RH
tget = function (x, penv=NULL, tenv=.PBSmodEnv) {
	if (is.null(penv)) penv = parent.frame() # need to call this inside the function NOT as an argument
	xnam = as.character(substitute(x))
	if (exists(xnam,envir=tenv)) {
		eval(parse(text=paste("tgot=get(\"",xnam,"\",envir=tenv)",sep="")))
		eval(parse(text=paste("assign(\"",xnam,"\",tgot,envir=penv)",sep="")))
		return(invisible(tgot)) # useful for calling remote functions
	}
	invisible()
}
tcall = function (x, penv=NULL, tenv=.PBSmodEnv) {
	if (is.null(penv)) penv = parent.frame() # need to call this inside the function NOT as an argument
	xnam = as.character(substitute(x))
	if (exists(xnam,envir=tenv)) {
		eval(parse(text=paste("tgot=get(\"",xnam,"\",envir=tenv)",sep="")))
		#if (class(tgot)=="function")
		return(invisible(tgot)) # useful for calling remote functions
	}
	invisible()
}
tprint = function (x, penv=NULL, tenv=.PBSmodEnv) {
	if (is.null(penv)) penv = parent.frame() # need to call this inside the function NOT as an argument
	xnam = as.character(substitute(x))
	if (exists(xnam,envir=tenv)) {
		eval(parse(text=paste("tgot=get(\"",xnam,"\",envir=tenv)",sep="")))
		print(tgot) # useful for jsust seeing objects in .PBSmodEnv
	}
	invisible()
}
tput = function (x, penv=NULL, tenv=.PBSmodEnv) {
	if (is.null(penv)) penv = parent.frame() # need to call this inside the function NOT as an argument
	xnam = as.character(substitute(x))
	if (exists(xnam,envir=penv))
		eval(parse(text=paste("assign(\"",xnam,"\",get(\"",xnam,"\",envir=penv),envir=tenv)",sep="")))
	invisible()
}
# versions of tget/tprint for window calls
.win.tget = function(){
	act = getWinAct()[1]
	eval(parse(text=paste("tget(",act,")()",sep="")))
}
.win.tcall = function(){
	act = getWinAct()[1]
	eval(parse(text=paste("tcall(",act,")()",sep="")))
}
.win.tprint = function(){
	act = getWinAct()[1]
	eval(parse(text=paste("tprint(",act,")()",sep="")))
}
# functions called from window description files
.win.runExHelperQuit = function(){ tcall(.runExHelperQuit)() }
.win.closeALL        = function(){ tcall(closeALL)() }
.win.closeSDE        = function(){ tcall(closeSDE)() }
.win.closeChoice     = function(){ tcall(.closeChoice)() }
.win.makeChoice      = function(){ tcall(.makeChoice)() }
.win.chFile          = function(){ tcall(chFile)() }
.win.chTest          = function(){ tcall(chTest)() }
.win.restoreCWD      = function(){
	if(exists("cwd",envir=.PBSmodEnv) && tcall(cwd)!=getwd())
		setwd(tcall(cwd))
}
#---------------------------tget/tcall/trpint/tput


#view-----------------------------------2011-10-31
# View first/last/random n element/rows of an object.
#-----------------------------------------------RH
view <- function (obj, n=5, last=FALSE, random=FALSE, print.console=TRUE, ...) {
	getn=function(n,N,last,random,...) {
		n=min(n,N)
		if (random) return(sample(1:N,n,...)) 
		n1=ifelse(last,N-n+1,1); n2=ifelse(last,N,n)
		return(n1:n2) }
	showVec=function(obj,n,last,random,...){
		N=length(obj); if (N==0) return("empty vector")
		v.vec=obj[getn(n,N,last,random,...)]
		return(v.vec) }
	showTab=function(obj,n,last,random,...){
		N=nrow(obj); if (N==0) return("empty table")
		mess = paste(c("v.tab=obj[getn(n,N,last,random,...)",rep(",",length(dim(obj))-1),"]"),collapse="")
		eval(parse(text=mess))
		return(v.tab) }
	showLis=function(obj,n,last,random,print.console,...){
		nL=length(obj); if (nL==0) return("empty list")
		v.lis=list()
		if (is.null(names(obj))) ii=1:nL else ii=names(obj)
		for (i in 1:nL) {
			iobj=obj[[i]]
			v.lis = c(v.lis,list(view(iobj,n,last,random,print.console=FALSE,...)))
		}
		names(v.lis)=ii; return(v.lis) }
	showAll=function(obj){
		return(obj) }
	# End Subfunction------------------------------
	if (n==0) return("nada")
	n=abs(n) # coerce to positive
	if (is.data.frame(obj) || is.matrix(obj) || is.array(obj)) 
		viewed=showTab(obj,n,last,random,...)
	else if (is.list(obj)) 
		viewed=showLis(obj,n,last,random,print.console,...)
	else if (is.vector(obj) || is.integer(obj) || is.numeric(obj) || is.character(obj)) 
		viewed=showVec(obj,n,last,random,...)
	else viewed=showAll(obj)
	if (print.console) print(viewed)
	invisible(viewed)
}
#---------------------------------------------view


#viewCode-------------------------------2015-01-20
# View package R code on the fly.
#-----------------------------------------------RH
viewCode=function(pkg="PBSmodelling", funs, output=4, ...){
	eval(parse(text=paste("if(!require(",pkg,",quietly=TRUE)) stop(\"",pkg," package is required\")",sep="")))
	tdir <- tempdir(); tdir <- gsub("\\\\","/",tdir)                    # temporary directory for R
	if (is.element(pkg,loadedNamespaces())){
		penv=asNamespace(pkg); pkgOb=ls(penv,all.names=TRUE)             # objects in package including those in namespace
		bad=regexpr(".__",pkgOb); pkgOb=pkgOb[bad<0]                     # get rid of names starting with ".__"
		delim=":::" }
	else {
		pkgOb=ls(paste("package:",pkg,sep=""),all.names=TRUE)            # package objects
		delim="::" }
	bad=regexpr("[>=~@:/&%!\\|()[{^$*+?<-]",pkgOb); pkgOb=pkgOb[bad<0]  # get rid of weird names
	z=sapply(pkgOb,function(x){
		if (is.element(x,c("break","for","function","if","next","repeat","while"))) return(FALSE) # special words in pkg 'base'
		eval(parse(text=paste("is.function(",paste(pkg,x,sep=delim),")",sep=""))) } )
	pkgF=names(z)[z]                                                    # package functions
	if (length(pkgF)==0) {showAlert(paste(pkg,"has no recognizable functions")); return()}
	dots=regexpr("^\\.",pkgF); pkgF0=pkgF[dots==1]; pkgF1=pkgF[dots!=1]
	code=c(paste("#",pkg,"Functions"),paste("#",paste(rep("-",nchar(pkg)+10),collapse="")))
	pkgFuns=c(pkgF1,pkgF0)
	if (missing(funs)) funs=pkgFuns
	else {
		if (!is.null(list(...)$pat) && is.logical(list(...)$pat) && list(...)$pat)
			funs = findPat(funs,pkgFuns) }
	if (is.null(funs) || is.na(funs) || !is.character(funs) || all(funs=="")) {
		showAlert("Your choice for 'funs' is badly specified")
		return(invisible("Error: 'funs' badly specified")) }
	seeFuns=pkgFuns[is.element(pkgFuns,funs)]
	if (length(seeFuns)==0) {
		showAlert("Your choices yield no functions")
		return(invisible("Error: choices for 'funs' yield no functions")) }
	if (output==2) {
		unpackList(list(...))
		# `home' code from `help.start'
		home <- 
		if (!is.element("remote",ls(envir=sys.frame(sys.nframe()))) || is.null(remote)) {
			if (as.numeric(R.version$svn) >= 67548)
				httpdPort = tools::startDynamicHelp(NA)
			else {
				showAlert("Output 2 only available for R (>= 3.2.0)")
				stop("Output 2 only available for R (>= 3.2.0)", call. = FALSE)
			}
			httpdPort = tools::startDynamicHelp(NA)
			if (httpdPort > 0L) {
				if ("update" %in% ls(envir=sys.frame(sys.nframe())) && update)
					make.packages.html(temp = TRUE)
				paste0("http://127.0.0.1:", httpdPort) 
			}
			else stop("`viewCode' requires the HTTP server to be running", call. = FALSE)
		}
		else remote
		expr=paste("iloc=paste(\"",home,"/library/",pkg,"/html/00Index.html\"); ",sep="")
		expr=paste(expr,"index=readLines(iloc);",sep="")
		eval(parse(text=expr))
	}
	for (i in c(seeFuns)) {
		if (output %in% c(1,2)) {
			expr=paste("fun=\"",i,"\"; ",sep="")
			if (output==2) {
				expr=c(expr,paste("helptit=index[grep(\">",i,"<\",index)+1][1]; ",sep=""))
				expr=c(expr,"if (!is.na(helptit)) {")
				expr=c(expr,"fun=paste(fun,substring(helptit,5,nchar(helptit)-10),sep=\"\\t\") }; ") }
		}
		else if (output %in% c(3)) {
			expr=paste("fun=deparse(args(",pkg,delim,i,")); fun=fun[!fun%in%\"NULL\"]; ",sep="")
			expr=c(expr,"fun=paste(fun,collapse=\"\"); ")
			expr=c(expr,"fun=gsub(\" \",\"\",x=fun); ")
			expr=c(expr,paste("fun=gsub(\"function\",\"",i," \",x=fun); ",sep="")) }
		else if (output %in% c(4)) {
			expr=paste("fun=deparse(",pkg,delim,i,"); ",sep="")
			expr=c(expr,paste("fun[1]=paste(\"",i,"\",fun[1],sep=\" = \",collapse=\"\"); ",sep="")) }
		else expr="fun=\"\"; "
		expr=paste(c(expr,"code=c(code,fun)") ,collapse="")
		eval(parse(text=expr)) }
	fname=paste(tdir,"/",pkg,ifelse(output %in% 1:2,".txt",".r"),sep="")
	writeLines(code, fname); openFile(convSlashes(fname))
	invisible(code) }
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~viewCode


#writePBSoptions------------------------2012-12-04
# Save PBS options to a text file
#  fname - name of options file (or path to this file)
#-------------------------------------------ACB/RH
writePBSoptions=function(fname="PBSoptions.txt") {
	.initPBSoptions()
	tget(.PBSmod)
	if(fname!="PBSoptions.txt")
		#packList(".optionsFile",".PBSmod$.options",fname) #.PBSmod$.options$.optionsFile<<-fname
		.PBSmod$.options$.optionsFile <- fname
	#packList(".optionsChanged",".PBSmod$.options",NULL) #.PBSmod$.options$.optionsChanged<<-NULL
	.PBSmod$.options$.optionsChanged <- NULL
	tput(.PBSmod)
	saveOpt=.PBSmod$.options[-grep("^[.]", names(.PBSmod$.options))]
	writeList(saveOpt, fname)
}
#----------------------------------writePBSoptions


#============HIDDEN FUNCTIONS=====================


#.addslashes----------------------------2006-08-28
# Escapes special characters from a string, which can then be used in the "P" format
# if x has more than one element, then it will returned a nested characterVector
# eg: c("it's", "O K") becomes => "'it\'s' 'O K'"
# Arguments:
#  x - string to escape
#  asvector - if true, always force strings to be printed as charactervector - e.g. "h w" => "'h w'"
#----------------------------------------------ACB
.addslashes <- function(x, asvector = FALSE ) {
	#escape backslashes
	x <- gsub("\\\\", "\\\\\\\\", x)
	#escase doublequotes
	x <- gsub("\"", "\\\\\"", x)
	#escase singlequotes
	x <- gsub("'", "\\\\'", x)
	#convert into substrings if applicable
	if (length(x)>1 || asvector == TRUE ) {
		i<-append(grep("[ \t\\\\]+", x), grep("^$", x)) #indicies needing quotes
		x[i]<-paste("'", x[i], "'", sep="")
		x<-paste(x, collapse=" ")
	}
	else {
		#special case where it is a single word with no special chars
		if (!any(grep("[ \t\\\\]+", x)) && x!="")
			return(x)
	}
	return(paste("\"", x, "\"", sep="")) }
#--------------------------------------.addslashes


#.convertVecToArray---------------------2006-09-16
# Converts a vector to an Array
# Arguments:
#  x       - a vector of data to use to create array
#  d       - dimensions of the array
#  byright - if TRUE, vary indices by the most right number first
#                 ex) 1,1 - 1,2 - 1,3 - 2,1 - 2,2 - 2,3
#            if FALSE, varry by most left (R default)
#                 ex) 1,1 - 2,1 - 1,2 - 2,2 - 1,3 - 2,3  
#----------------------------------------------ACB
.convertVecToArray <- function(x,d, byright=FALSE, byrow=TRUE) {
	if (length(x)!=prod(d))
		stop("given vector x length does not match product of dimensions")
	#create an empty array with correct dimensionality
	y=vector(mode(x),prod(d));  dim(y)=d
	#iterate over every possible index
	pts<-.getArrayPts(d,byright=byright,byrow=byrow)
	for(i in 1:length(pts)) {
		arrIndex <- paste(pts[[i]],collapse=",")
		if (byrow) vecIndex=i
		else       vecIndex=.mapArrayToVec(pts[[i]], d, byright)
#print(paste(vecIndex,arrIndex,sep=" : "))
#browser();return()
		#map it to the appropriate place in the given X vector
		code = paste("y[", arrIndex, "] <- ", x[vecIndex], sep="")
		eval(parse(text=code))
	}
	return(y)
}
#-------------------------------.convertVecToArray


#.findSquare----------------------------2009-02-11
.findSquare=function(nc) {
	sqn=sqrt(nc); m=ceiling(sqn); n=ceiling(nc/m)
	return(c(m,n)) }
#--------------------------------------.findSquare


#.forceMode-----------------------------2009-02-11
# Forces a variable into a mode without showing any warnings
# Arguments:
#   x    - variable to convert
#   mode - mode to convert to
#-----------------------------------------------RH
.forceMode <- function(x, mode) {
	ex1=paste("xcon=as.",mode,"(x)",sep="")          # explicit conversion
	ex2=paste("xcon=",paste(x,collapse=""),sep="")   # character representation of object
	warn=options()$warn; options(warn=-1)
	try(eval(parse(text=ex1)),silent=TRUE)
	if (!exists("xcon",envir=environment()) || all(is.na(xcon)))
		try(eval(parse(text=ex2)),silent=TRUE)
	if (!exists("xcon",envir=environment()) || all(is.na(xcon))) xcon=x
	options(warn=warn)
	return(xcon) }
#---------------------------------------.forceMode


#.getArrayPts---------------------------2009-02-10
#  Returns all possible indices of an array
# Arguments:
#  d is a vector of integers specifing the dimensions
# output: a list of vectors of all possible indices
#-------------------------------------------ACB/RH
.getArrayPts <- function(d, byright=FALSE, byrow=TRUE) {
	x<-list(); m=length(d);
	if (m>2) mbyrow=c(2,1,3:m) else mbyrow=c(2,1)
	for(i in 1:length(d)) {
		if (byright) x[[i]] = d[i]:1
		else x[[i]] = 1:d[i] }
	if (byrow) {
		xx=expand.grid(c(x[2],x[-2])); xx=xx[,mbyrow] }
	else xx=expand.grid(x)
#browser();return()
	y<-list()
	for(i in 1:length(xx[[1]])) {
		z<-unlist(xx[i,])
		attributes(z)<-NULL
		y[[i]]<-z }
	return(y) }
#-------------------------------------.getArrayPts


#.initPBSoptions------------------------2012-12-04
# Called from zzz.R's .First.lib() intialization function
# Attach will place only a copy of the environment on the search path;
# so not useful for manipulating .PBSmod
#--------------------------------------------AE/RH
.initPBSoptions <- function() {
	if (exists(".PBSmod",envir=.GlobalEnv)){
		assign(".PBSmod",.PBSmod,envir=.PBSmodEnv)
		rm(.PBSmod,envir=.GlobalEnv) }                 # Can no longer modify user's global environment
	if (!exists(".PBSmod",envir=.PBSmodEnv))
		assign(".PBSmod",list(),envir=.PBSmodEnv)      #.GlobalEnv) #.PBSmod <<- list()
	tget(.PBSmod)
	if (is.null(.PBSmod$.options))
		#packList(".options",".PBSmod",list())          #.PBSmod$.options <<- list()
		.PBSmod$.options <- list()
	if (is.null(.PBSmod$.options$openfile))
		#packList("openfile",".PBSmod$.options",list()) #.PBSmod$.options$openfile <<- list()
		.PBSmod$.options$openfile <- list()
	tput(.PBSmod)
}
#----------------------------------.initPBSoptions


#.mapArrayToVec-------------------------2006-09-16
# Determines which index to use for a vector, when given an 
# N-dim index of an array.
# Arguments:
#  x       - array index (numeric vector)
#  d       - dimensions of the array
#  byright - if true, vary most right indices first, 
#            if false, vary by left (R default)
#----------------------------------------------ACB
.mapArrayToVec <- function(x,d, byright=FALSE) {
	x <- x - 1 #start counting at 0 instead of 1
	m <- length(x)
	if (m!=length(d))
		stop("given points (x), does not match lenght of given dimensions (d)")
	if (byright) {
		ind <- x[m]
		for(i in (m-1):1) {
			ind <- ind+x[i]*prod(d[(i+1):m]) }
	}
	else {
		ind <- x[1]
		for(i in 2:length(d)) {
			ind <- ind+x[i]*prod(d[(i-1):1]) 
#print(ind+1)
#browser();return()
			}
	}
	return(ind+1)

	#return(x[1] + d[1]*x[2])
}
#-----------------------------------.mapArrayToVec


#.removeFromList------------------------2008-10-06
# Remove items from a list.
#--------------------------------------------AE/RH
.removeFromList = function (l, items) {
	if (!length(l) || !length(items))  return(l)
	keep = l[!is.element(names(l),items)]
	return(keep) }
#----------------------------------.removeFromList


#.tclArrayToVector----------------------2008-10-06
.tclArrayToVector <- function( str )
{
	#strings (when multiple) are encoded as "{c:/program files/somefile.txt} c:/nospaces/isok.txt"
	file_vector <- c()
	quoted <- FALSE
	start <- -1 #negative means no current start
	chars <- strsplit( str, "" )[[ 1 ]]
	for( i in 1:nchar( str ) ) {
		if( quoted == TRUE ) {
			if( chars[ i ] == "}" ) {
				end <- i - 1
				s <- substr( str, start, end )
				file_vector <- c( file_vector, s )
				quoted <- FALSE
				start <- -1
			}
			#do nothing otherwise - just accecpt the input
		} else {
			#not quoted
			if( chars[ i ] == " " && start != -1 ) {
				end <- i - 1
				s <- substr( str, start, end )
				file_vector <- c( file_vector, s )
				start <- -1
			} else if( chars[ i ] == "{" ) {
				quoted <- TRUE
				start <- i + 1
			} else if( start < 0 ) {
				start <- i
			}
		}
	}
	if( start > 0 ) {
		end <- nchar( str )
		s <- substr( str, start, end )
		file_vector <- c( file_vector, s )
	}
	return( file_vector )
}
#--------------------------------.tclArrayToVector

#===== THE END ===================================
