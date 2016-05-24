#admb-----------------------------------2014-02-18
# Starts the primary GUI interface
# Authors: Jon T. Schnute, Rowan Haigh, Alex Couture-Beil
#-----------------------------------------------RH
admb <- function(prefix="", wdf="admbWin.txt", optfile="ADopts.txt", pathfile="ADpaths.txt"){
	.initOptions()
	#readADopts()
	if (!is.null(pathfile) && file.exists(pathfile))
		readADpaths(pathfile)
	else pathfile="ADpaths.txt"
	pkg="PBSadmb"
#browser();return()

	#TODO rename to something else - too similar to .PBSadmb
	assign("PBSadmb",list(pkg=pkg,call=match.call(),args=args(admb),useCols=NULL),envir=.PBSadmbEnv)

	pdir <- system.file(package=pkg)          # package directory
	wdir <- paste(pdir,"/win",sep="")         # window description file directory
	edir <- paste(pdir,"/examples",sep="")    # examples directory
	tdir <- tempdir()                         # temporary working directory
	twdf <- paste(tdir,"/",wdf,sep="")        # temporary window description file
	twdf <- convSlashes(twdf,os="unix")

	stripExt = function(x) { return(sub("[.].{1,3}$", "", x)) }

	win.filename <- paste(wdir,wdf,sep="/")
	temp <- readLines(win.filename)
	
	# insert examples into window description file menuitems
	etpl <- basename(Sys.glob(file.path(edir,"*.tpl"))) # TPL files in examples directory
	eprf <- stripExt(etpl)                              # strip off extensions
	enew <- character(0)
	edir <- gsub( "\\\\", "/", edir )
	for (i in eprf) 
		enew=c(enew,paste("menuitem label=",i," function=doAction action=\"copyFiles(`",
			i,".`,srcdir=`",edir,"`); convOS(paste(`",i,"`,c(`.tpl`,`.dat`,`.pin`,`.r`),sep=``))\"",sep=""))
	temp <- gsub("@nitems",length(eprf),temp)
	temp <- gsub("@menuitems",paste(enew,collapse="\n\t"),temp)

	temp = gsub("@wdf",twdf,temp)
	temp = gsub("@pathfile",pathfile,temp)
	
	#create the window (from temp string)
	temp <- unlist( strsplit(temp, "\n" ) )
	#createWin(twdf, TRUE)
	writeLines(temp,con=twdf)
	createWin(twdf)

	#set some values
	.load.prefix.droplist()
	loadOptionsGUI( atcall(.PBSadmb) )
	isOK <- .win.checkADopts()

	if( isOK == FALSE && .Platform$OS.type == "windows" ) {
		cat( "\nADMB or MinGW are not installed in the default location.\n" )
		cat( "If needed, follow install instructions from GUI menu <Install>\n" )
		cat( "If ADMB is installed on your system, you can manually set the ADMB or MinGW path values\n" )
		cat( "in the GUI to point to your own installations.\n" )
	}

	#TODO need centralized window variable init (is it done anywhere?)
	setWinVal( list( currentdir.values = getwd() ) )
	setWinVal( list( currentdir = getwd() ) )
	#setWinVal( list( optfile = optfile ) )
	setWinVal( list( optfile = pathfile ) )

	invisible() }
#---------------------------------------------admb

#=================================================
#              ADMB FUNCTIONS                     
#=================================================

#convAD---------------------------------2014-02-25
# Conver TPL file to CPP code.
#-------------------------------------------JTS/RH
convAD <- function(prefix, raneff=FALSE, safe=TRUE, dll=FALSE, debug=FALSE, logfile=TRUE, add=TRUE, verbose=TRUE, pathfile=NULL)
{
	#get path and name of program
	ext <- ifelse( .Platform$OS.type == "windows", ".exe", "" )
	prog <- ifelse( raneff == TRUE, "tpl2rem", "tpl2cpp" )
	prog <- paste(prog,ext,sep="")

	#add cmd flags
	flags <- c()
	if( dll )
		flags[ length( flags ) + 1 ] <- "-dll"
	if( safe )
		flags[ length( flags ) + 1 ] <- "-bounds"

	#build command string
	flags <- paste( flags, collapse=" " )

	cmd <- paste( prog, flags, prefix, sep=" " )
	if (.Platform$OS.type=="windows")
		cmd=shQuote(cmd)
		#cmd=.addQuotes(convSlashes(cmd))

	#add ADMB path to path env variable
	old_path <- Sys.getenv( "PATH" )
	.setPath(pathfile)

	#pre cmd run
	if (logfile & !add)
		startLog(prefix)
	if (verbose)
		cat(cmd,"\n")

	#run cmd
	tplout <- .callSys(cmd)
	tplout2 <- c(cmd,tplout)

	#post cmd run
	if (logfile) {
		#tplout2 <- c(cmd,tplout)
		appendLog(prefix, tplout2)
	}
	if (verbose)
		cat(tplout, sep="\n")

	#restore path
	Sys.setenv( PATH = old_path )
	invisible(tplout2)
}
.win.convAD=function(winName="PBSadmb")
{
	isOK=.win.checkADopts()
	if (!isOK)
		return()
	setWinVal(list(Mtime=matrix(NA,nrow=3,ncol=3)),winName=winName)
	time0=proc.time()[1:3]
	getWinVal(scope="L",winName=winName)
	convAD(prefix, raneff, safe, dll, debugsymbols, logfile, add, verbose)
	Ttime=round(proc.time()[1:3]-time0,2)
	setWinVal(list("Mtime[1,1]"=Ttime[1],"Mtime[1,2]"=Ttime[2],"Mtime[1,3]"=Ttime[3]),winName=winName) 
	if( .Platform$OS.type == "unix" )
		cat("\n> ")
	invisible(Ttime)
}

#compAD---------------------------------2014-02-25
# Apparently "raneff" doesn't influence the compile stage,
# but the argument is preserved here for future development.
#-------------------------------------------JTS/RH
compAD <- function(prefix, raneff=FALSE, safe=TRUE, dll=FALSE, debug=FALSE, logfile=TRUE, add=TRUE, verbose=TRUE, pathfile=NULL)
{
	#get path and name of program
	admbpath = getOptions(atcall(.PBSadmb),"admbpath")
	ext <- ifelse( .Platform$OS.type == "windows", ifelse(file.exists(paste(admbpath,"/bin/adcomp.cmd",sep="")),".cmd",".bat"), "" )
	#ext <- ifelse( .Platform$OS.type == "windows", ".bat", "" )
	prog <- paste( "adcomp", ext, sep="" )

	#add cmd flags
	if (is.null(getOptions(atcall(.PBSadmb),"admbver")))
		setADMBVer(gccver = NULL) # only get ADMB version
	admbvernum = .version(getOptions(atcall(.PBSadmb),"admbver"))
	flags <- c()
	if( dll )
		flags[ length( flags ) + 1 ] <- "-d"
	if( debug )
		flags[ length( flags ) + 1 ] <- "-g"
	if( safe && admbvernum < 11 )
		flags[ length( flags ) + 1 ] <- "-s"
	if( !safe && admbvernum >= 11 )
		flags[ length( flags ) + 1 ] <- "-f"
	if( raneff )
		flags[ length( flags ) + 1 ] <- "-r"

	#build command string
	flags <- paste( flags, collapse=" " )
	cmd <- paste( prog, flags, prefix, sep=" " )
	if (.Platform$OS.type=="windows")
		cmd=shQuote(cmd)
		#cmd=.addQuotes(convSlashes(cmd))

	#add ADMB path to path env variable
	old_path <- Sys.getenv( "PATH" )
	.setPath(pathfile)

	#pre cmd run
	if (logfile & !add)
		startLog(prefix)
	if (verbose)
		cat(cmd,"\n")

	#run cmd
	out <- .callSys(cmd)
	out2 <- c(cmd,out)

	#post cmd run
	if (logfile) {
		#out2 <- c(cmd,out)
		appendLog(prefix, out2)
	}
	if (verbose)
		cat(out, sep="\n")

	Sys.setenv( PATH = old_path )
	invisible(out2)
}
.win.compAD=function(winName="PBSadmb")
{
	isOK=.win.checkADopts()
	if (!isOK) return()
	time0=proc.time()[1:3]
	getWinVal(scope="L",winName=winName)
	compAD(prefix, raneff, safe, dll, debugsymbols, logfile, add, verbose)
	Ttime=round(proc.time()[1:3]-time0,2)
	setWinVal(list("Mtime[2,1]"=Ttime[1],"Mtime[2,2]"=Ttime[2],"Mtime[2,3]"=Ttime[3]),winName=winName) 
	if( .Platform$OS.type == "unix" ) cat("\n> ")
	invisible(Ttime)
}

#linkAD---------------------------------2014-02-25
# Links binaries into executable
#-------------------------------------------JTS/RH
linkAD <- function(prefix, raneff=FALSE, safe=TRUE, dll=FALSE, debug=FALSE, logfile=TRUE, add=TRUE, verbose=TRUE, pathfile=NULL)
{
	#get path and name of program
	admbpath = getOptions(atcall(.PBSadmb),"admbpath")
	ext <- ifelse( .Platform$OS.type == "windows", ifelse(file.exists(paste(admbpath,"/bin/adlink.cmd",sep="")),".cmd",".bat"), "" )
	#ext <- ifelse( .Platform$OS.type == "windows", ".bat", "" )
	prog <- paste( "adlink", ext, sep="" )

	#add cmd flags
	if (is.null(getOptions(atcall(.PBSadmb),"admbver")))
		setADMBVer(gccver = NULL) # only get ADMB version
	admbvernum = .version(getOptions(atcall(.PBSadmb),"admbver"))
	flags <- c()
	if( dll )
		flags[ length( flags ) + 1 ] <- "-d"
	if( debug )
		flags[ length( flags ) + 1 ] <- "-g"
	if( safe && admbvernum < 11 )
		flags[ length( flags ) + 1 ] <- "-s"
	if( !safe && admbvernum >= 11 )
		flags[ length( flags ) + 1 ] <- "-f"
	if( raneff )
		flags[ length( flags ) + 1 ] <- "-r"

	#build command string
	flags <- paste( flags, collapse=" " )
	cmd <- paste( prog, flags, prefix, sep=" " )
	if (.Platform$OS.type=="windows")
		cmd=shQuote(cmd)
		#cmd=.addQuotes(convSlashes(cmd))

	#add ADMB path to path env variable
	old_path <- Sys.getenv( "PATH" )
	.setPath(pathfile)

	#pre cmd run
	if (logfile & !add)
		startLog(prefix)
	if (verbose)
		cat(cmd,"\n")

	#run cmd
	out <- .callSys(cmd)
	out2 <- c(cmd,out)

	#post cmd run
	if (logfile) {
		#out2 <- c(cmd,out)
		appendLog(prefix, out2)
	}
	if (verbose)
		cat(out, sep="\n")

	Sys.setenv( PATH = old_path )
	invisible(out2)
}
.win.linkAD=function(winName="PBSadmb")
{
	isOK=.win.checkADopts()
	if (!isOK) return()
	time0=proc.time()[1:3]
	getWinVal(scope="L",winName=winName)
	linkAD(prefix, raneff, safe, dll, debugsymbols, logfile, add, verbose)
	Ttime=round(proc.time()[1:3]-time0,2)
	setWinVal(list("Mtime[3,1]"=Ttime[1],"Mtime[3,2]"=Ttime[2],"Mtime[3,3]"=Ttime[3]),winName=winName) 
	if( .Platform$OS.type == "unix" ) cat("\n> ")
	invisible(Ttime)
}

#makeAD---------------------------------2014-02-25
# Convert TPL file to CPP code.
# Compile CPP to object files.
# Links binaries into executable.
#-------------------------------------------JTS/RH
makeAD <- function(prefix, raneff=FALSE, safe=TRUE, dll=FALSE, debug=FALSE, logfile=TRUE, add=TRUE, verbose=TRUE, pathfile=NULL)
{
	convAD(prefix, raneff, safe, dll, debug, logfile, add, verbose, pathfile) # first call using pathfile is sufficient
	compAD(prefix, raneff, safe, dll, debug, logfile, add, verbose)
	linkAD(prefix, raneff, safe, dll, debug, logfile, add, verbose)
}
.win.makeAD=function(winName="PBSadmb") {
	isOK=.win.checkADopts()
	if (!isOK) return()
	.win.convAD( winName )
	.win.compAD( winName )
	.win.linkAD( winName )
}

runAD <- function(prefix, argvec="", logfile=TRUE, add=TRUE, verbose=TRUE)
{
	if( .Platform$OS.type == "windows" )
		p.exe <- paste(prefix,".exe",sep="")
	else
		p.exe <- paste("./", prefix, sep="" )
	if (logfile) {
		p.log=paste(prefix,".log",sep="")
		p.log.log=paste(p.log,".log",sep="")
		if (file.exists(p.log)) file.copy(p.log,p.log.log,overwrite=TRUE) }
	if (all(argvec==""))
		p.cmd <- p.exe
	else
		p.cmd <- paste(p.exe, paste(argvec,collapse=" ",sep=""), sep=" ")

	p.err <- paste("File",p.exe,"does not exist.\n",sep=" ")

	if (.Platform$OS.type=="windows") {
		p.cmd=shQuote(p.cmd)
		#p.cmd=.addQuotes(convSlashes(p.cmd))
	}
	if (file.exists(p.exe)) p.out <- .callSys(p.cmd) else p.out <- p.err

	if (logfile) {
		if (!add) startLog(prefix)
		else if (file.exists(p.log.log)) file.copy(p.log.log,p.log,overwrite=TRUE)
		appendLog(prefix, p.out) }
	if (verbose) cat(p.out, sep="\n")
	invisible(p.out)
}
.win.runAD=function(winName="PBSadmb")
{
	setWinVal(list(Rtime=matrix(NA,nrow=1,ncol=3)),winName=winName)
	time0=proc.time()[1:3]
	getWinVal(scope="L",winName=winName)
	runAD(prefix=prefix,argvec=argvec,logfile=logfile,add=add,verbose=verbose)
	Ttime=round(proc.time()[1:3]-time0,2)
	setWinVal(list("Rtime[1,1]"=Ttime[1],"Rtime[1,2]"=Ttime[2],"Rtime[1,3]"=Ttime[3]),winName=winName) 
	invisible(Ttime)
}

runMC <- function(prefix, nsims=2000, nthin=20, outsuff=".mc.dat",
                  logfile=FALSE, add=TRUE, verbose=TRUE) {
	outf <- paste(prefix,outsuff,sep="")

	arg1 <- paste("-mcmc",.asIs(nsims),"-mcsave",.asIs(nthin),sep=" ")

	arg2 <- "-mceval"

	runAD(prefix, arg1, logfile=logfile, add=add, verbose=verbose)

	p.out <- runAD(prefix, arg2, logfile=logfile, add=TRUE, verbose=verbose)

	writeLines(p.out,outf)
	invisible(p.out)
}
.win.runMC=function(winName="PBSadmb")
{
	time0=proc.time()[1:3]
	getWinVal(scope="L",winName=winName)
	arg1 <- paste("-mcmc",.asIs(nsims),"-mcsave",.asIs(nthin),sep=" ")

	setWinVal(list(argvec=arg1))
	runMC(prefix=prefix,nsims=nsims,nthin=nthin,logfile=logfile,add=add,verbose=verbose)
	Ttime=round(proc.time()[1:3]-time0,2)
	setWinVal(list("Rtime[1,1]"=Ttime[1],"Rtime[1,2]"=Ttime[2],"Rtime[1,3]"=Ttime[3]),winName=winName)
	#eval(parse(text="PBSadmb$useCols <<- NULL"))
	atget(PBSadmb); PBSadmb$useCols <- NULL; atput(PBSadmb)
	invisible(Ttime)
}

.win.run=function(winName="PBSadmb")
{
	getWinVal(scope="L",winName=winName)
	if (runType=="mcmc") .win.runMC()
	else .win.runAD() 
	if( .Platform$OS.type == "unix" ) cat("\n> ")
	invisible()
}

#=================================================
#              OPTIONS FUNCTIONS                  
#=================================================

.initOptions <- function()
{
	#don't re-init
	if(exists(".PBSadmb", envir=.PBSadmbEnv) && is(atcall(.PBSadmb),"PBSoptions"))
		return()
	readADopts()
}

#checkADopts----------------------------2014-02-18
# Checks path locations specified by ADMB options.
#-----------------------------------------------RH
checkADopts=function(opts=getOptions( atcall(.PBSadmb) ),
   check=c("admbpath","gccpath","editor"), warn=TRUE, popup=FALSE)
{
	sAF = options()$stringsAsFactors
	on.exit(options(stringsAsFactors=sAF))
	options(stringsAsFactors=FALSE)
	isWin = .Platform$OS.type=="windows" #; isWin=FALSE
	slash = ifelse(isWin, "\\", "/" )
	
	### assume g++ is always available on Unix machines, so don't check it.
	#if (!isWin) opts = opts[setdiff(names(opts),"gccpath")]

	### Check that .PSadmb has all required paths and that links point to actual files on the hard drive.
	### check for admb, mingw, and editor programs ###
	mess=list()
	for (i in names(opts)) {
		if (!any(i==check))
			next
		ii=ipath=opts[[i]]
		if (i=="admbpath") {
			ipath=paste(ii,slash,"bin",sep="")
			progs=paste(c("tpl2cpp","tpl2rem"),ifelse(isWin,".exe",""),sep="")
		}
		else if (i=="gccpath") 
			progs=paste("bin",slash,"g++",ifelse(isWin,".exe",""),sep="")
		else if (i=="editor") {
			ipath=dirname(ii)
			progs=basename(ii)
		}
		target=paste(ipath,progs,sep=slash)
		if (!isWin && i=="gccpath") {
			ipath = "/usr"; istatus = TRUE }
		else if (all(target==slash)) {
			progs="editor"; istatus = FALSE }
		else
			istatus = file.exists(target)
		names(istatus)=progs
#if (i=="editor") {browser();return()}
		mess[[ipath]]=istatus
	}
	
	ADstatus=all(unlist(mess)==TRUE)
	attr(ADstatus,"status")=mess
	vmess=unlist(mess)
	names(vmess)=paste(rep(names(mess),sapply(mess,length,simplify=FALSE)),
		unlist(sapply(mess,names,simplify=FALSE)),sep=slash)
	attr(ADstatus,"message")=vmess
	if (warn|popup) {
		if (all(vmess==TRUE)) {
			if(warn) cat("All programs found\n\n") }
			#if(popup) showAlert("All programs found","Continue","info") }
		else {
			badmess=paste("Programs not found:\n----------------------\n",paste(names(vmess)[!vmess],collapse="\n"),
				"\n\nAlter the path file (default 'ADpaths.txt') in the working directory.\n",
				"     ~~~OR~~~\n",
				"If using the PBSadmb GUI, alter the path entries in the Setup tab.\n\n",sep="")
			if (isWin && popup) {
				badmess <- paste( badmess, "If you need to install ADMB, follow links from Install dropdown menu.\n\n", sep="" )
			}
			if (warn) cat(badmess)
			if (popup) showAlert(badmess,"User action required","warning") 
		}
	}
	### check for sed.exe when all programs above are found ###
	if (isWin && all(vmess)) {
		opts.sed = c(list(msyspath=paste(opts[["gccpath"]],slash,"msys",slash,"bin",sep="")),opts[c("admbpath","gccpath")])
		sedmess=list()
		for (i in names(opts.sed)) {
			if (!any(i==c(check,"msyspath")))
				next
			ii=ipath=opts.sed[[i]]
			if (i=="msyspath")
				progs = paste("sed",ifelse(isWin,".exe",""),sep="")
			else if (i=="admbpath")
				progs=paste("bin",slash,"sed",ifelse(isWin,".exe",""),sep="")
			else if (i=="gccpath")
				progs=paste("bin",slash,"sed",ifelse(isWin,".exe",""),sep="")
			target=paste(ipath,progs,sep="/")
			istatus=file.exists(target)
			names(istatus)=progs
			sedmess[[ipath]]=istatus
		}
		smess=unlist(sedmess)
		names(smess)=paste(rep(names(sedmess),sapply(sedmess,length,simplify=FALSE)),
			unlist(sapply(sedmess,names,simplify=FALSE)),sep=slash)
		if (warn|popup) {
			if (any(smess)) {
				# if(warn) cat(paste("'sed",ifelse(isWin,".exe'",",")," program found\n\n",sep="")) 
			}
			else {
			#if (!any(smess)) {
				badsedmess=paste("Exception: the program `sed",ifelse(isWin,".exe",""),"` was not found on any of these paths:\n",
					paste(gsub(paste("sed",ifelse(isWin,".exe",""),sep=""),"",names(smess))[!smess],collapse="\n"),
					"\n\nPlace a copy of `sed",ifelse(isWin,".exe",""),"` on any one of the paths indicated above.\n\n",sep="")
				if (warn) cat(badsedmess)
				if (popup) showAlert(badsedmess,"User action required","warning") 
			}
		}
	}
	invisible(ADstatus)
}
.win.checkADopts=function(winName="PBSadmb")
{
	getWinVal(scope="L",winName=winName)
	chkstat=checkADopts(opts=list(admbpath=admbpath,gccpath=gccpath,editor=editor),popup=TRUE)
	#set label to OK/FIX with coloured background
	setWinVal(list(chkstat=ifelse(chkstat," OK"," Fix")),winName=winName)
	setWidgetColor( "chkstat", winName=winName, bg=ifelse(chkstat,"lightgreen","pink") )
	setWidgetColor( "checkbutton", winName=winName, bg=ifelse(chkstat,"moccasin","pink") )
	invisible(chkstat)
}

#makeADopts-----------------------------2009-08-12
# Sets the ADMB path directories (deprecated).
#----------------------------------------------ACB
makeADopts <- function( admbpath, gccpath, editor )
{
	cat( "THIS FUNCTION IS DEPRECATED - use setADMBPath\n" )
	setADMBPath( admbpath = admbpath, gccpath = gccpath, editor = editor )
}
.win.makeADopts=function(winName="PBSadmb")
{
	getWinVal(scope="L",winName=winName)
	setADMBPath(admbpath,gccpath,editor) 
	.win.setADMBVer() 
	.win.checkADopts()
	invisible()
}

#readADopts-----------------------------2009-08-12
# Read ADMB options from a file.
#-----------------------------------------------RH
readADopts <- function(optfile="ADopts.txt")
{
	# Create instance of option manager - use this to get/set/save/load options
	# First attempt to load options from the package, then attempt to load options from the current dir (which will override pkg options)
	pkg_fname = paste( system.file(package="PBSadmb"), "/ADopts.txt", sep="" )
	#eval(parse(text=".PBSadmb.pkgOptions <<- new( \"PBSoptions\", filename = pkg_fname, initial.options = list(admbpath=\"\", gccpath=\"\",editor=\"\"), gui.prefix=\"\" )"))
	.PBSadmb.pkgOptions <- new( "PBSoptions", filename = pkg_fname, initial.options = list(admbpath="", gccpath="",editor=""), gui.prefix="" )
	atput(.PBSadmb.pkgOptions)

	# Load from current dir, using pkgOptions as default values
	#eval(parse(text=".PBSadmb <<- new( \"PBSoptions\", filename = optfile, initial.options = getOptions( .PBSadmb.pkgOptions ), gui.prefix=\"\" )"))
	.PBSadmb <- new( "PBSoptions", filename = optfile, initial.options = getOptions( .PBSadmb.pkgOptions ), gui.prefix="" )

	.guessPath <- function( programs, includefilename = FALSE, failed = NULL ) {
		for( p in programs ) {
			found <- Sys.which( p )[ 1 ]
			if( found != "" ) {
				if( includefilename )
					return( found )
				return( dirname( found ) )
			}
		}
		return( failed )
	}

	#if( getOptions( .PBSadmb, "editor" ) == "" )
	#	setOptions( .PBSadmb, editor = .guessPath( c( "kate", "notepad" ), TRUE ) )
	atput(.PBSadmb)
	invisible()
}
.win.readADopts=function(winName="PBSadmb")
{
	getWinVal(scope="L",winName=winName)
	if (file.exists(optfile)) {
		readADopts(optfile=optfile)
		loadOptionsGUI( atcall(.PBSadmb) )
	} else {
		mess=paste("Options file '",optfile,"' does not exist",sep="")
		showAlert(mess)
	stop(mess) }
	invisible()
}

#readADpaths----------------------------2014-02-25
# Read ADMB paths from a simple 2-column file.
# Allows user to pass in a file easily made by hand.
# Assumes .PBSadmb options object exists.
#-----------------------------------------------RH
readADpaths = function(pathfile) {
	sAF = options()$stringsAsFactors
	on.exit(options(stringsAsFactors=sAF))
	options(stringsAsFactors=FALSE)
	if (!missing(pathfile) && !is.null(pathfile) && file.exists(pathfile)) {
		ufile = read.table(file=pathfile,header=FALSE,col.names=c("target","path"))
		uopts = split(ufile$path,ufile$target)
		atget(.PBSadmb)
		setOptions(.PBSadmb,uopts)
		atput(.PBSadmb)
	}
}
.win.readADpaths = function(winName="PBSadmb"){
	pathfile = getWinVal()$optfile
	readADpaths(pathfile)
	loadOptionsGUI( atcall(.PBSadmb) )
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~readADpaths

#saveADpaths----------------------------2014-02-27
# Save ADMB paths to a simple 2-column file.
# Assumes .PBSadmb options object exists.
#-----------------------------------------------RH
saveADpaths = function(pathfile) {
	sAF = options()$stringsAsFactors
	on.exit(options(stringsAsFactors=sAF))
	options(stringsAsFactors=FALSE)
	if (missing(pathfile) || is.null(pathfile)) pathfile="ADpaths.txt"
	isWin = .Platform$OS=="windows"
	Uopts = getOptions(atcall(.PBSadmb))
	upath = if (isWin) c("admbpath","gccpath","editor") else c("admbpath","editor")
	uopts = Uopts[upath]
	ufile = t(sapply(upath,function(x,u){
		c(x,paste(rep(" ",10 - nchar(x)),collapse=""),convSlashes(u[[x]],addQuotes=TRUE))
		},u=uopts,USE.NAMES=FALSE))
	write.table(ufile,file=pathfile,row.names=FALSE,col.names=FALSE,quote=FALSE,sep="")
}
.win.saveADpaths = function(winName="PBSadmb"){
	pathfile = getWinVal()$optfile
	saveADpaths(pathfile)
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~saveADpaths

#setupAD--------------------------------2014-02-26
# Command line initialization to read in path
# information and check for the presence of executables:
# admbpath -- tpl2cpp, tpl2rem
# gccpath  -- g++
# editor   -- e.g. notepad, but can be any valid editor.
#-----------------------------------------------RH
setupAD = function(pathfile){
	if (missing(pathfile) || is.null(pathfile) || !file.exists(pathfile))
	pathfile = "ADpaths.txt"
	.initOptions()
	readADpaths(pathfile)
	checkADopts()
}

#writeADopts----------------------------2009-08-12
# Writes ADMB options to a file.
#-----------------------------------------------RH
writeADopts <- function(optfile="ADopts.txt")
{
	#save to current dir
	atget(.PBSadmb)
	saveOptions( .PBSadmb, optfile )

	#save to pkg dir (don't change fname)
	opts <- getOptions( .PBSadmb )
	atget(.PBSadmb.pkgOptions)
	setOptions(.PBSadmb.pkgOptions,opts)
	atput(.PBSadmb.pkgOptions)
	#hack no longer needed: http://code.google.com/p/pbs-modelling/issues/detail?id=81 solved
	#tmp <- list( atcall(.PBSadmb.pkgOptions) )
	#tmp[ names( opts ) ] <- opts
	#do.call( setOptions, tmp )
	saveOptions( .PBSadmb.pkgOptions )
	return(invisible(NULL))
}
.win.writeADopts=function(winName="PBSadmb")
{
	isOK=.win.checkADopts()
	if (!isOK) return()
	getWinVal(scope="L",winName=winName)
	writeADopts(optfile=optfile) 
	invisible()
}

#setADMBPath----------------------------2009-08-12
# Sets the ADMB path directories.
#----------------------------------------------ACB
setADMBPath <- function( admbpath, gccpath, editor )
{
	.initOptions()
	atget(.PBSadmb)
	if( missing( admbpath ) == FALSE )
		setOptions( .PBSadmb, admbpath = admbpath )
	if( missing( gccpath ) == FALSE )
		setOptions( .PBSadmb, gccpath = gccpath )
	if( missing( editor ) == FALSE )
		setOptions( .PBSadmb, editor = editor )
	atput(.PBSadmb)
}
.win.setADMBPath=function(winName="PBSadmb")
{
	getWinVal(scope="L",winName=winName)
	setADMBPath(admbpath,gccpath,editor) 
	.win.checkADopts()
	invisible()
}


#.setPath-------------------------------2014-02-25
# Set the temporary environment path to run ADMB.
#-------------------------------------------ACB/RH
.setPath <- function(pathfile)
{
	path_sep <- .Platform$path.sep
	dir_sep <- ifelse( .Platform$OS.type == "windows", "\\", "/" )
	
	# User can specify a 2-column pathfile (no headers)
	if (!missing(pathfile) && !is.null(pathfile) && file.exists(pathfile))
		readADpaths(pathfile)

	admb_home <- getOptions( atcall(.PBSadmb), "admbpath" )
	admb_path <- paste( admb_home, "bin", sep = dir_sep )

	if( .Platform$OS.type == "windows" ) {
		gcc_path <- paste( getOptions( atcall(.PBSadmb), "gccpath"), "bin", sep = dir_sep )
		msys_path <- paste( getOptions( atcall(.PBSadmb), "gccpath"), "msys", "bin", sep = dir_sep )
		path <- paste(.normPath(admb_path,dir_sep), .normPath(msys_path,dir_sep), .normPath(gcc_path,dir_sep), sep=path_sep)
	} else {
		#linux must include original path so programs like cat, sed are found
		sys_path <- Sys.getenv( "PATH" )
		path <- paste(.normPath(admb_path,dir_sep), .normPath(sys_path,dir_sep), sep=path_sep)
	}
	Sys.setenv( PATH = path )
	Sys.setenv( ADMB_HOME = gsub("/*$","",gsub("\\\\*$","",.normPath(admb_home,dir_sep)) ) ) #ensure no trailing slash (`/') or (`\\') exists
}

#add a dir to path variable (only if it hasn't already been set)
.appendToPath <- function( path_to_add )
{
	path <- Sys.getenv( "PATH" )
	path_sep <- .Platform$path.sep
	if( any( unlist( strsplit( path, path_sep ) ) == path_to_add ) == FALSE ) {
		#path_to_add doesn't exist in path - append it, and reset env variable
		path <- paste( path, path_to_add, sep = path_sep )
		Sys.setenv( PATH = path )
	}
}

#setADMBVer-----------------------------2014-02-18
# Sets the ADMB versions.
#-----------------------------------------------RH
setADMBVer <- function( admbver, gccver )
{
	.initOptions()
	isWin = .Platform$OS.type=="windows" #; isWin=FALSE
	sayWhat = attributes(checkADopts(warn=FALSE))$status
	atget(.PBSadmb)
	opts = getOptions(.PBSadmb)
	if( !missing( admbver ) && is.null( admbver ) ) 
		junk = "do nothing"
	else if( !missing( admbver ) && !is.element(admbver,c("")) ) 
		setOptions( .PBSadmb, admbver = admbver )
	else {
		if(all(sayWhat[[1]])) { # check ADMB version
			if(file.exists(paste(opts["admbpath"],"/VERSION",sep="")))
				setOptions(.PBSadmb, admbver = readLines(paste(opts["admbpath"],"/VERSION",sep=""))[1] )
			else if (is.null(getOptions(.PBSadmb,"admbver")))
				setOptions(.PBSadmb, admbver = "10")
	}	}
	if( !missing( gccver ) && is.null( gccver ) ) 
		junk = "do nothing"
	else if( !missing( gccver ) && !is.element(gccver,c("")) )
		setOptions( .PBSadmb, gccver = gccver )
	else {
		if(all(sayWhat[[2]])) { # check g++ version
			cmd = "g++ --version"
			if (isWin) cmd = shQuote(paste(opts["gccpath"],"/bin/",cmd,sep=""))
			gccVer = .callSys(cmd)[1]
			gccpcs = strsplit(gccVer," ")[[1]]
			gccver = gccpcs[grep("\\.",gccpcs)][1]
#browser()
			#gccver =  paste(setdiff(strsplit(gsub("[[:alpha:]+()]","",gccver)," ")[[1]],""),collapse=".")
			setOptions(.PBSadmb, gccver = gccver)
	}	}
	atput(.PBSadmb)
}
.win.setADMBVer=function(winName="PBSadmb")
{
	getWinVal(scope="L",winName=winName)
	setADMBVer(admbver,gccver) 
	for (i in c("admbver","gccver")) {
		ival = getOptions(atcall(.PBSadmb),i)
		if (is.null(ival)) next
		eval(parse(text=paste("setWinVal(list(",i,"=getOptions(atcall(.PBSadmb),\"",i,"\")),winName=winName)",sep="")))
	}
	#.win.checkADopts()
	invisible()
}

showADargs <- function(prefix,ed=TRUE) {
	if( .Platform$OS.type == "windows" )
		p.exe <- paste(prefix,".exe", sep="")
	else
		p.exe <- paste("./", prefix, sep="") #TODO verify
		# p.exe <- paste(prefix,".exe", sep=""); # JTS: Wrong for Unix; should have been deleted
	p.arg <- paste(prefix,".arg", sep="")

	p.err <- paste("File",p.exe,"does not exist.\n",sep=" ")

	p.cmd <- paste(p.exe,"-?",sep=" ")

	#p.cmd=.addQuotes(convSlashes(p.cmd))
	p.cmd=shQuote(p.cmd)
	if (file.exists(p.exe)) p.out <- .callSys(p.cmd) else p.out <- p.err

	if (ed) {writeLines(p.out,p.arg)
		editADfile(p.arg)
	} else {
		cat(paste(p.out,collapse="\n"))
		cat(paste(p.arg,collapse="\n"))
	}
	invisible(p.out) }

.win.showADargs=function(winName="PBSadmb") {
	getWinVal(scope="L",winName=winName)
	showADargs(prefix=prefix) 
	invisible() }

.win.showGUIargs=function(winName="PBSadmb"){ # display the argument vector in the GUI based on radio selection
	getWinVal(scope="L",winName=winName)
	if (runType=="normal") setWinVal(list(argvec=""),winName=winName)
	else if (runType=="mcmc") setWinVal(list(argvec=paste("-mcmc",.asIs(nsims),"-mcsave",.asIs(nthin),sep=" ")),winName=winName)
	else if (runType=="lprof") setWinVal(list(argvec="-lprof"),winName=winName)
	else setWinVal(list(argvec=argvec),winName=winName) 
	invisible() }

#=================================================
#              FILE FUNCTIONS                     
#=================================================

#copyFiles------------------------------2013-03-25
# Copy files with specified prefixes and suffixes 
# from one location to another.
#-----------------------------------------------RH
copyFiles=function(prefix,suffix=NULL,srcdir=getwd(),dstdir=getwd(),ask=TRUE){
	if (missing(prefix)) return()
	if (is.null(prefix) || prefix=="*") prefix=""
	if (is.null(suffix) || suffix=="*") suffix=""
	prefix=gsub("\\.","\\\\\\.",prefix)
	suffix=gsub("\\.","\\\\\\.",suffix)
	npref=length(prefix)
	nsuff=length(suffix)
	fpatt=paste(rep(prefix,each=nsuff),rep(suffix,npref),sep="")
	fname=sapply(fpatt,function(x,dir){list.files(dir,x,ignore.case=TRUE)},dir=srcdir,simplify=FALSE)
	fname=unique(unlist(fname,use.names=FALSE))
	nfile=length(fname)
	fname=list.files(srcdir,pattern=fpatt,ignore.case=TRUE)
	if (nfile==0) return()
	copy.out=rep(FALSE,nfile)
	names(copy.out)=fname
	for (i in 1:nfile){
		fnam0=paste(srcdir,"/",fname[i],sep="")
		fnam1=paste(dstdir,"/",fname[i],sep="")
		if (!file.exists(fnam1)) {ask=FALSE; ovr=TRUE}
		if (ask) ovr=getYes(paste("Overwrite",fname[i],"?"))
		copy.out[i]=file.copy(fnam0,dstdir,overwrite=ovr) 
	}
	if (exists(".PBSmod",envir=.PBSmodEnv) && tcall(.PBSmod)$.activeWin=="PBSadmb") 
		setWinVal(list(prefix=substring(prefix,1,nchar(prefix)-2)),winName="PBSadmb")
	invisible(copy.out) }
#----------------------------------------copyFiles

#editADfile-----------------------------2013-03-25
# Use the specified editor to edit a file.
#-------------------------------------------JTS/RH
editADfile <- function(fname)
{
	if (!attributes(checkADopts(warn=FALSE))$status[[3]]) {cat("Valid editor not verified in Setup tab\n"); stop()}
	#if (!checkADopts(warn=FALSE)) {cat("Invalid options for PBSadmb\n"); stop()}
	#f.edit <- paste("start \"\"",.addQuotes(convSlashes(.ADopts$editor)),.addQuotes(convSlashes(fname)),sep=" ");
	if (.Platform$OS.type=="windows") {
		f.edit <- paste(shQuote(getOptions(atcall(.PBSadmb),"editor")),shQuote(fname),sep=" ")
		#f.edit <- paste(.addQuotes(convSlashes(getOptions(.PBSadmb,"editor"))),.addQuotes(convSlashes(fname)),sep=" ")
	} else {
		f.edit <- paste(getOptions(atcall(.PBSadmb),"editor"),fname,sep=" ")

	}
	f.err  <- paste("File",fname,"does not exist.\n",sep=" ")

	if (file.exists(fname)) {
		.callSys(edit=f.edit, wait=FALSE)
		cat(f.edit,"\n")
		f.out <- TRUE
	} else {
		cat(f.err)
		f.out <- FALSE
	}
	return(f.out)
}

editAD <- function(prefix, suffix=c(".tpl",".cpp",".log"))
{
	npref=length(prefix)
	nsuff=length(suffix)
	ed.out=logical(npref*nsuff)
	k=0
	for (i in 1:npref) {
		for (j in 1:nsuff) {
			fname=paste(prefix[i],suffix[j],sep="")
			k=k+1
			ed.out[k]=editADfile(fname)
		}
	}
	return(ed.out) }

.win.editAD=function(winName="PBSadmb") {
	getWinVal(scope="L",winName=winName)
	editAD(prefix=prefix)
	invisible() }
.win.editPLT=function() {
	pref=findPrefix(".plt")
	editAD(prefix=pref,suffix=".plt")
	invisible() }

startLog <- function(prefix)
{
	p.log <- paste(prefix, ".log", sep="");
	if (file.exists(p.log)) file.remove(p.log)

	dstamp <- date()

	line1 <- paste("Log file for ", prefix, " (", dstamp, ")\n", sep="")
	writeLines(line1, con=p.log)
}
.win.startLog=function(winName="PBSadmb")
{
	getWinVal(scope="L",winName=winName)
	startLog(prefix=prefix) 
	invisible()
}

appendLog <- function(prefix, lines)
{
	p.log <- paste(prefix, ".log", sep="")

	if (!file.exists(p.log)) startLog(prefix)

	cat(lines, file=p.log, sep="\n", append=TRUE)
}
.win.appendLog=function(winName="PBSadmb")
{
	getWinVal(scope="L",winName=winName)
	appendLog(prefix=prefix,lines=lines) 
	invisible()
}

#readRep--------------------------------2009-07-09
# Imports the generated reports in various ways
#-----------------------------------------------RH
readRep=function(prefix, suffix=c(".cor",".rep",".std",".mc.dat"), global=FALSE) {
	findFormat=function(dat){ # extracted from Alex's readList()
		for (i in 1:length(dat)) {
			if (!any(grep("^[ \t]*[#`]", dat[i]))) {
				if (any(grep("^[ \t]*structure", dat[i]))) fileformat="D"
				else if (any(grep("^[ \t]*list", dat[i]))) fileformat="R"
				else if (any(grep("^[ \t]*\\$", dat[i])))  fileformat="P"
				else fileformat="U" # unknown fileformat detected
				break
		}	}
	return(fileformat) }
	sAF = options()$stringsAsFactors
	on.exit(options(stringsAsFactors=sAF))
	options(stringsAsFactors=FALSE)
	if (missing(prefix) || any(is.null(c(prefix,suffix))) || any(c(prefix,suffix)=="")) return()
	flist=list()
	fname=paste(prefix,suffix,sep="")
	if (global) cat("R objects created:\n")
	for (i in fname) {
		if(!file.exists(i)) next
		mtime=file.info(i)$mtime  # date & time file was modified/made
		if (global) cat(paste("     ",i,"\n",sep=""))
		ii=substring(i,nchar(prefix)+2)
		contents=dat=readLines(i)
		ncont=length(contents)
		ff=findFormat(contents)
		if (ff=="P") dat=PBSmodelling::.readList.P(i)
		else if (ff=="D" || ff=="R") dat=eval(parse(i))
		else if (any(ii==c("cor","std","mc.dat"))) { # treat these as a data frames
			tcont=gsub("std dev","std",contents)
			if (ii=="cor") tcont[2]=paste("index name value std ",paste(1:(ncont-2),collapse=" "))
			tfile=paste(prefix,".tmp",sep="")
			writeLines(tcont,tfile)
			skip=ifelse(any(ii==c("cor")),1,0)
			header=ifelse(any(ii==c("cor","std")),TRUE,FALSE)
			fill=ifelse(any(ii==c("cor","std")),TRUE,FALSE)
			dat=read.table(tfile,skip=skip,header=header,fill=fill)
			#---look for MCMC names---
			if (ii=="mc.dat") { 
				report=NULL
				report=try(readList(paste(prefix,".rep",sep="")),silent=TRUE)
				if (!is.null(report) && any(names(report)=="mcest")) {
					if (length(report$mcest)==ncol(dat))
						dat=rbind(report$mcest,dat) }
				if (!is.null(report) && any(names(report)=="mcnames")) {
					if (length(report$mcnames)==ncol(dat))
						names(dat)=report$mcnames }
			}
			#---correlation matrix manipulation---
			if (ii=="cor") { 
				ldh=as.numeric(substring(tcont[1],regexpr("=",tcont[1])+1)) # log determinant of hessian
				name=NULL
				NAME=dat[,2]
				value=dat[,3]
				std=dat[,4]
				names(NAME)=1:length(NAME) # to keep the original order after split and sapply
				ldupe=split(NAME,NAME)
				ndupe=sapply(ldupe,function(x){
					if(length(x)>1) {
						xnam=names(x)
						len=length(x)
						npad=ifelse(any(len==10^(0:10)),1,0)
						npad=ceiling(log10(len))+npad
						x=paste(x,pad0(1:len,npad),sep="")
						names(x)=xnam 
					}
					return(x) },simplify=FALSE)
				for (j in names(ndupe)) name=c(name,ndupe[[j]])
				name=name[order(as.numeric(names(name)))]
				dat=dat[,setdiff(names(dat),c("index","name","value","std"))]
				names(dat)=name; row.names(dat)=name
				dat[is.na(dat)]=t(dat)[is.na(dat)]
				attr(dat,"determinant")=ldh
				attr(dat,"value")=value
				attr(dat,"std")=std
			} #---end correlation manipulation---
		}
		if (ii!="mc.dat") attr(dat,"contents")=contents
		attr(dat,"mtime")=mtime
		expr=paste("if (",global,") assign(\"",i,"\",dat,envir=.GlobalEnv); flist$",i,"=dat",sep="") 
		eval(parse(text=expr))
	}
	if (length(flist)==1) flist=flist[[1]] # return a single object, not a list of objects
	invisible(flist)
}
#------------------------------------------readRep

.win.readRep=function(winName="PBSadmb") {
	getWinVal(scope="L",winName=winName)
	act=getWinAct()[1]
	if (sum(toView)>0) {
		suffix=paste(".",names(toView)[toView],sep="")
		readRep(prefix=prefix,suffix=suffix,global=TRUE) 
	}
	if (pltView) {
		pref=findPrefix(".plt")
		if (length(pref)>0){
			for (i in pref) readRep(prefix=i,suffix=".plt",global=TRUE) } 
	}
	invisible()
}

.win.viewRep=function(winName="PBSadmb") {
	getWinVal(scope="L",winName=winName)
	act=getWinAct()[1]
	if (!is.null(act) && act=="allnone") {
		if (allnone==1) {
			toView=rep(TRUE,length(toView))
			pltView=TRUE
		}
		if (allnone==0) {
			toView=rep(FALSE,length(toView))
			pltView=FALSE
		} 
	}
	if (!is.null(act) && act=="check") allnone=3
	setWinVal(list(allnone=allnone,toView=toView,pltView=pltView),winName=winName)
	if (!is.null(act) && act=="open") {
		if (sum(toView)>0) {
			suffix=paste(".",names(toView)[toView],sep="")
			editAD(prefix,suffix) 
		}
		if (pltView) .win.editPLT() 
	}
	invisible()
}

#convOS---------------------------------2009-11-20
# Convert text files to the default format
# of the operating system.
#-----------------------------------------------RH
convOS = function(inam, onam=inam, path=getwd())
{
	if (missing(inam)) stop("Supply names(s) of text file(s) to convert")
	else if(length(inam)!=length(onam)) stop("Number of 'inam's does not match number of 'onam's")
	N=length(inam)
	inam=paste(path,inam,sep="/")
	onam=paste(path,onam,sep="/")
	for (i in 1:N) {
		if(file.exists(inam[i])) {
			idat=readLines(inam[i])
			writeLines(idat,con=onam[i]) 
	} }
	invisible()
}

#=================================================
#              PLOT FUNCTIONS                     
#=================================================

#plotMC---------------------------------2009-02-06
# Plots the MCMC output in various ways
#-----------------------------------------------RH
plotMC=function(prefix,act="pairs",pthin=1,useCols=NULL){
	if (is.null(prefix) || prefix=="") return()
	inFile=paste(prefix, ".mc.dat", sep="")
	if(!file.exists(inFile) || file.info(inFile)$size==0){
		showAlert(paste("Cannot find file", inFile, "in working directory."))
		return() 
	}
	if (!exists(inFile,envir=.GlobalEnv) || is.null(attributes(get(inFile))$mtime) ||
			file.info(inFile)$mtime!=attributes(get(inFile))$mtime)
		inData=readRep(prefix,".mc.dat")
	else inData=get(inFile)
	if (nrow(inData)==1) {
		showAlert(paste("No MCMC data for '",prefix,
			"'\n\nSpecify an 'mceval_phase()' in your 'tpl' file.",sep=""))
		return() 
	}
	x=1:nrow(inData)
	pthin=max(pthin,1)
	xthin=seq(1,nrow(inData),pthin)
	if (pthin>1) inData=inData[xthin,]
	if (!is.null(useCols)) inData=inData[,useCols]
	nc=ncol(inData)
	puce="#cc8899"
	clrs=c("blue", "red", "green", "magenta","navy",puce)
	clrs=rep(clrs,nc)[1:nc]
	resetGraph()
	
	panel=function(x,y) {
		len=length(x)
		points(x[2:len],y[2:len],pch=21,col="grey",bg="gainsboro",cex=0.8)
		points(x[1],y[1],pch=21,bg=puce,cex=1.2)
	}
	if (act=="pairs") pairs(inData, panel=panel,gap=0)
	if (act=="eggs") plotFriedEggs(inData)
	if (act=="acf") {
		expandGraph(mfrow=c(1,1),mar=c(3,3.5,.5,.5))
		plotACF(inData, clrs=rep(c("blue","red"),nc)[1:nc],lwd=ifelse(nc>6,1,2))
		mtext("Correlation",side=2,line=2.25,cex=1.2)
		mtext("Lags",side=1,line=1.5,cex=1) 
	}
	if (act=="trace") {
		inData=cbind(x=xthin,inData)
		expandGraph(mfrow=c(nc,1),mar=c(0,0,0,0),oma=c(4,4.5,.5,.5))
		for (i in 1:nc) {
			plotTrace(inData[,c(1,i+1)],clrs=clrs[i],xaxt="n")
			axis(1,labels=par()$mfg[1]==par()$mfg[3])
			mtext(names(inData)[i+1],side=2,line=3,cex=1) 
		}
		mtext("Sequential chain values",side=1,outer=TRUE,line=2,cex=1) 
	}
	if (act=="dens") {
		rc=PBSmodelling::.findSquare(nc)
		#sqn=sqrt(nc); m=ceiling(sqn); n=ceiling(nc/m)
		expandGraph(mfrow=c(rc[1],rc[2]),mar=c(2,2,0,0),oma=c(1,1.75,.5,.5),mgp=c(1.5,.2,0))
		for (i in 1:nc) {
			plotDens(inData[,i],clrs=clrs[i]) 
			addLabel(0.95,0.95,names(inData)[i],adj=c(1,1),cex=1.2) 
		}
		#normData=sweep(inData,2,apply(inData,2,calcGM),"/")
		mtext("Kernel Density",outer=TRUE,side=2,line=0.2,cex=1.2) 
		mtext("Parameter estimates",outer=TRUE,side=1,line=-0.5,cex=1)
	}
	invisible() }
#-------------------------------------------plotMC

.win.plotMC=function(winName="PBSadmb") {
	getWinVal(scope="L",winName=winName)
	if (is.null(prefix) || prefix=="") return()
	act=getWinAct()[1]
	if (is.null(act)) act="pairs"
	plotMC(prefix=prefix,act=act,pthin=pthin,useCols=atcall(PBSadmb)$useCols) 
	invisible() }

.win.viewCode=function(winName="PBSadmb",pkg="PBSadmb"){
	getWinVal(scope="L",winName=winName)
	tdir <- tempdir()
	tdir <- gsub("\\\\","/",tdir)  # temporary directory for R
	pkgN=match(paste("package:",pkg,sep=""),search()) # position number of package
	pkgO=ls(pos=pkgN) # package objects
	z=sapply(pkgO,function(x){f=get(x)
	is.function(f)}); pkgF=names(z)[z] # package functions
	code=c("# PBSadmb Functions","#------------------")
	for (i in pkgF) {
		expr=paste("fun=deparse(",pkg,"::",i,"); fun[1]=paste(\"",i,"\",fun[1],sep=\" = \",collapse=\"\"); code=c(code,fun)",sep="")
		eval(parse(text=expr)) 
	}
	fname=paste(tdir,"/",pkg,".r",sep="")
	writeLines(code, fname)
	editADfile(fname)
	invisible() }

#=================================================
#              CLEAN-UP FUNCTIONS                 
#=================================================

#cleanAD--------------------------------2009-07-22
# Clean files with given prefix.
#-------------------------------------------JTS/RH
cleanAD <- function(prefix=NULL) {
	ofile=findPat(c("^admodel\\.","^eigv\\.","^sims$","^variance$","^tmp_","^fmin\\.log$",
		"^tfile$","\\.tmp$","\\.bak$","^hessian\\.bin$","^hesscheck$","^diags$","^dgs2$"),
		list.files(all.files=TRUE,ignore.case=TRUE))
	if (is.null(prefix)) {
		tpl=findPrefix(".tpl")
		apat=paste("^",tpl,"\\.",sep="") # all tpls
		bpat=paste("\\.",c("tpl","dat","pin","r","pdf","mc.dat"),"$",sep="") # keep files
		afile=sapply(apat,function(x){list.files(pattern=x,ignore.case=TRUE)},simplify=FALSE)
		afile=as.vector(unlist(afile))
		bfile=sapply(bpat,function(x){list.files(pattern=x,ignore.case=TRUE)},simplify=FALSE)
		bfile=as.vector(unlist(bfile))
		pfile=setdiff(afile,bfile)
		ofile=union(ofile,pfile)  # remove duplicated names
		if (length(ofile)==0) ofile=NULL
		isOK=.cleanWD(files=ofile)
		if (!isOK) closeWin("cleanWindow") 
	}
	else {
		pfile=list.files(pattern=paste("^",prefix,"\\.",sep=""),ignore.case=TRUE)
		if (length(pfile)==0) psuff=NULL
		else {
			psuff=substring(pfile,nchar(prefix)+1)
			psuff=setdiff(psuff,c(".tpl",".dat",".pin",".r",".pdf",".mc.dat")) }
		ofile=setdiff(ofile,paste(prefix,psuff,sep="")) # remove duplicated names
		if (length(ofile)==0) ofile=NULL
		isOK=.cleanUp(prefix, suffix=psuff, files=ofile)
		if (!isOK) closeWin("cleanWindow") 
	}
}
#------------------------------------------cleanAD

.win.cleanAD=function(winName="PBSadmb") {
	getWinVal(scope="L",winName=winName)
	cleanAD(prefix=prefix) 
	invisible()
}

#.cleanUp-------------------------------2009-02-11
# Anisa's cleanProj function modified for flexibility.
#-----------------------------------------------RH
.cleanUp=function(prefix, suffix, files) {
	if (missing(suffix)) suffix = character(0)
	if (missing(files))  files  = character(0)
	rowLen = ceiling(sqrt(max(length(suffix), length(files))))
	if (rowLen == 0) return(invisible(FALSE))
	winDesc = c("window name=cleanWindow title=Clean",
		"grid 1 3",
		paste("entry name=cleanPrefix value=\"", prefix, "\" label=Prefix ",
			"mode=character width=12 font=\"bold 9\"", sep = ""),
		"button text=\">\" function=.win.findClean",
		"button text=refresh function=.cleanUpAgain",
		"label text=\"\n\nSuffixes to Clean\" font=\"bold 9\"", 
		PBSmodelling::.makeCleanVec("suff", suffix, rowLen), 
		"label text=\"\n\nFiles to Clean\" font=\"bold 9\"", 
		PBSmodelling::.makeCleanVec("file", files, rowLen), 
		"grid 1 3 relief=groove padx=4 pady=4", 
		"button function=.selectCleanBoxes action=1 text=\"Select All\" padx=4 pady=4", 
		"button function=.selectCleanBoxes action=0 text=\"Deselect All\" padx=4 pady=4", 
		"button function=.doClean text=Clean bg=aliceblue padx=4 pady=4")
	createWin(winDesc, astext = TRUE, env=PBSmodelling::.getHiddenEnv() ) 
	invisible(TRUE)
}

.cleanUpAgain=function(winName="cleanWindow"){
	cleanAD(getWinVal(winName=winName)$cleanPrefix)
	invisible()
}

.win.findClean=function(winName="cleanWindow"){
	choice=findPrefix(".tpl") 
	chooseWinVal(choice,"cleanPrefix",winname=winName) 
	invisible()
}

#.cleanWD-------------------------------2013-03-25
# Clean all potential grabage files.
#-----------------------------------------------RH
.cleanWD=function(files){ # Clean all nuisance files
	rowLen = ceiling(sqrt(length(files)))
	if (rowLen == 0) {
		try(closeWin("cleanWD"),silent=TRUE)
		return(invisible(FALSE)) 
	}
	winDesc = c("window name=cleanWD title=Clean",
		"label text=\"\n\nFiles to Clean\" font=\"bold 9\"",
		PBSmodelling::.makeCleanVec("file", files, rowLen),
		"grid 1 3 relief=groove padx=4 pady=4", 
		"button function=.selectCleanBoxes action=1 text=\"Select All\" padx=4 pady=4", 
		"button function=.selectCleanBoxes action=0 text=\"Deselect All\" padx=4 pady=4", 
		"button function=doAction text=Clean bg=aliceblue padx=4 pady=4 action=\".doCleanADMB(); closeWin(`cleanWD`)\"")
	createWin(winDesc, astext = TRUE, env=PBSmodelling::.getHiddenEnv() )
	invisible(TRUE)
}

#.doCleanADMB---------------------------2013-03-25
# Anisa's .doClean function modified for file names only
#-----------------------------------------------RH
.doCleanADMB=function () { 
	vec=getWinVal(scope="L")
	vecList=logical()
	for (i in names(vec)) vecList=c(vecList,vec[[i]])
	filenames = names(vecList)[vecList]
	if (!length(filenames)) 
		showAlert("No files to delete.")
	else if (getYes(paste("Delete ", paste(filenames, collapse = ", "), "?", sep = ""))) 
		file.remove(filenames)
	remaining = file.exists(filenames)
	if (sum(remaining)) 
		showAlert(paste("Failed to delete", paste(filenames[remaining], collapse = ", ")))
}


#=================================================
#              UTILITY FUNCTIONS                  
#=================================================

#atget----------------------------------2013-03-25
# Provide PBSadmb wrappers for PBSmodelling 
#  functions tget/tcall/tprint/tput/lisp
#-----------------------------------------------RH 
atget   = function(...) {tget  (..., penv=parent.frame(), tenv=.PBSadmbEnv)}
atcall  = function(...) {tcall (..., penv=parent.frame(), tenv=.PBSadmbEnv)}
atprint = function(...) {tprint(..., penv=parent.frame(), tenv=.PBSadmbEnv)}
atput   = function(...) {tput  (..., penv=parent.frame(), tenv=.PBSadmbEnv)}
alisp   = function(...) {lisp  (..., pos =.PBSadmbEnv)}


#suggestPath----------------------------2014-02-19
# Suggesta path for a specified program from PATH
# Arguments:
#  progs = vector of program names (without extension)
#  ipath = initial path by user to try before PATh directories
#  file_ext = alternative program extension if other
#    than `.exe` (primarily for Windows users)
#-----------------------------------------------RH
suggestPath <- function(progs, ipath=NULL, file_ext=NULL)
{
	isWin = .Platform$OS.type=="windows" #; isWin=FALSE
	path_sep = .Platform$path.sep
	file_sep = ifelse(isWin, "\\", "/" )
	if (is.null(file_ext))
		file_ext = ifelse(isWin, ".exe", "" )
	sys_path=Sys.getenv( "PATH" )
	paths = strsplit(sys_path,path_sep)[[1]]

	status = logical()
	pathos = list()
	for (aprog in progs) {
		aprog = paste(aprog,file_ext,sep="")
		inprog = sapply(c(ipath,paths),function(x){file.exists(paste(x,aprog,sep=file_sep))})
		pathos[[aprog]] = inprog
		if (!any(inprog)) {
			target  = "" 
			astatus = FALSE
		} else {
			target = names(inprog)[inprog][1]
			astatus = TRUE
		}
		names(astatus) = target
		status = c(status,astatus)
	}
	attr(status,"pathos") = pathos
	return(status)
}
#.win.suggestPath-----------------------2014-02-20
# Function called by GUI to suggest paths for setup.
#-----------------------------------------------RH
.win.suggestPath=function(winName="PBSadmb")
{
	getWinVal(scope="L",winName=winName)
	isWin = .Platform$OS.type=="windows" #; isWin=FALSE
	file_sep = ifelse(isWin, "\\", "/" )
	file_ext = ifelse(isWin, ".exe", "" )
	suggestions = list()
	nPath =function(x){
		if (is.null(x) || x=="") return(NULL)
		else gsub("(\\\\|/)$", "", .normPath(x) )
	}
	if (isWin) {
		admbpath.default = "C:\\admb"
		gccpath.default  = "C:\\mingw"
	} else {
		admbpath.default = "/usr/local/admb"
		gccpath.default  = "/usr"
	}

	admb_sugg = suggestPath("tpl2cpp",paste(c(nPath(admbpath),admbpath.default),"bin",sep=file_sep))
	if (admb_sugg)
		suggestions[["admbpath"]] = gsub(paste(file_sep,"bin",sep=ifelse(isWin,file_sep,"")),"",names(admb_sugg))
	else
		suggestions[["admbver"]] = ""
#browser();return()

	gcc_sugg = suggestPath("g++",paste(c(nPath(gccpath),gccpath.default),"bin",sep=file_sep))
	if (gcc_sugg)
		suggestions[["gccpath"]] = gsub(paste(file_sep,"bin",sep=ifelse(isWin,file_sep,"")),"",names(gcc_sugg))
	else
		suggestions[["gccver"]] = ""

	if (isWin) editors = c("Uedit32","notepad++","cedt","runemacs","notepad")
	else       editors = c("gedit","Vim","vi")
	edit_sugg = FALSE
	### what if user specifies an editor with some weird extension (e.g., non.sense.old, fake_emacs)
	if (!is.null(editor) && editor!="") {
		path_editor = dirname(editor)
		user_editor = basename(editor)
		edpcs = strsplit(user_editor,"\\.")[[1]]
		if (length(edpcs)==1) user_ext = ""
		else {
			user_ext = paste(".",rev(edpcs)[1],sep="")
			user_editor = paste(rev(rev(edpcs)[-1]),collapse=".")
		}
		edit_sugg = suggestPath(user_editor,nPath(path_editor),file_ext=user_ext)
	}
	if (edit_sugg)
		suggestions[["editor"]] = paste(names(edit_sugg)[edit_sugg][1],file_sep,user_editor,user_ext,sep="")
	else {
		edit_sugg = suggestPath(editors,nPath(editor))
		if (any(edit_sugg))
			suggestions[["editor"]] = paste(names(edit_sugg)[edit_sugg][1],file_sep,editors[edit_sugg][1],file_ext,sep="")
	}
#browser();return()
	if (length(suggestions)>0)
		setWinVal(suggestions,winName=winName)
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~suggestPath

.addQuotes=function(str){
	return(paste("\"",str,"\"",sep=""))
}

.asIs=function(x) {
	if (is.numeric(x)) x=format(x,scientific=FALSE)
	return(x)
}

.callSys <- function(..., wait=TRUE) { # note dots is not interpreted as expected, uses only first in list
	dots=list(...)
	if (.Platform$OS.type=="windows") {
		if("edit"%in%names(dots))
			dots[[1]]=paste("start \"\"",dots[[1]],sep=" ")
		out=shell(dots,intern=TRUE) }
	else {
		cmd <- unlist(list(...))
		#print( cmd )
		if( wait == FALSE )
			out=system(cmd,wait=FALSE)
		else
			out=system(cmd,intern=TRUE)
	}
	invisible(out)
}

.changeWD <- function( wd )
{
	if( missing( wd ) )
		wd <- selectDir() #tclvalue(tkchooseDirectory())
	if( wd != "" ) {
		currentdir.values <- sort( unique( c( wd, getWinVal()$currentdir.values ) ) )
		setWinVal( list( currentdir.values = currentdir.values ) )
		setWinVal( list( currentdir = wd ) )
		setwd( wd )
		.load.prefix.droplist()
	}
}

#called by droplist when user hits enter
.changeWDEnter <- function()
{
	wd <- getWinVal()$currentdir
	#remove trailing slash
	wd <- gsub("(\\\\|/)$", "", wd )
	if( file.exists( wd ) )
		.changeWD( wd )
	else {
		showAlert( paste( "unable to set working directory to \"", wd, "\" - does not exist", sep="" ) )
		.changeWD( getwd() )
	}
}

.chooseCols=function(winName="PBSadmb")
{
	getWinVal(scope="L",winName=winName)
	if (is.null(prefix) || prefix=="")
		return()
	inFile=paste(prefix, ".mc.dat", sep="")
	if(!file.exists(inFile)){
		showAlert(paste("Cannot find file", inFile, "in working directory."))
		return()
	}
	if (!exists(inFile,envir=.GlobalEnv) || is.null(attributes(get(inFile))$mtime) ||
			file.info(inFile)$mtime!=attributes(get(inFile))$mtime)
		inData=readRep(prefix,".mc.dat")
	else
		inData=get(inFile)
	flds=names(inData)
	nc=length(flds)
	#assign("PBSadmb",PBSadmb)
	useCols=atcall(PBSadmb)$useCols
	if (is.null(useCols) || length(useCols)!=nc)
		useCols=rep(TRUE,nc)
	#store feilds as a data.frame - for use with scrolling object widget
	choices <- as.data.frame( useCols )
	rownames( choices ) <- flds
	#converts data.frame scrolling object back into vector and saves to global var
	saveCols <- function()
	{
		choices <- getWinVal(winName="chooseCols")$choices
		atget(PBSadmb)
		PBSadmb$useCols=choices[[ 1 ]]
		atput(PBSadmb)
		#assign("PBSadmb",PBSadmb,envir=.GlobalEnv)
		closeWin("chooseCols")
	}
	toggleSelected <- function()
	{
		winName <- "chooseCols" 
		choices <- getWinVal(winName = winName )$choices
		if( any( choices[[1]] ) )
			choices[[1]] <- choices[[1]] & FALSE #force all false
		else
			choices[[1]] <- choices[[1]] | TRUE #force all true
		setWinVal( list( choices = choices ), winName = winName )
	}
	winDesc = c("window name=chooseCols title=Choose",
		"object choices rowshow=20",
		"grid 1 2",
		"button text=\"select all/none\" func=toggleSelected padx=\"0 10\"",
		"button text=OK bg=skyblue function=saveCols" )
		
	createWin(winDesc, astext = TRUE)
}

#given a vector of paths, return the last directory name of each path
.getDirName <- function( path )
{
	dirname <- gsub("\\\\", "/", path )
	return( unlist( lapply( strsplit( dirname, "/" ), function(x) x[length(x)] ) ) )
}

.load.prefix.droplist <- function()
# Repopulates droplist with all prefixes in current directory
# TODO we might want a "refresh" button on the GUI to call this
# ideally we could have the droplist call a function *BEFORE* 
# it drops down - thus refreshing the list
{
	choices <- findPrefix( ".tpl" )
	setWinVal( list( prefix.values = choices ) )
	if (any(choices=="vonb")) ch1 = grep("vonb",choices) else ch1 = 1
	setWinVal( list( prefix = choices[ ch1 ] ) )
}

.normPath = function(path, winslash="\\", mustWork=FALSE) {
	normalizePath( path, winslash, mustWork )
}

.version = function(x) {
	if (is.null(x) || is.numeric(x)) return(x)
	xpc = strsplit(x,split="\\.")[[1]]
	npc = length(xpc)
	xnu = as.numeric(paste(xpc[1:min(2,npc)],collapse="."))
	return(xnu)
}

.win.findTPL=function(suffix=".tpl",winName="PBSadmb"){ 
	choice=findPrefix(suffix) 
	chooseWinVal(choice,"prefix",winname=winName) 
	invisible()
}

#.win.viewCode--------------------------2009-02-13
# View the package R code on the fly.
#-----------------------------------------------RH
.win.viewCode=function(pkg="PBSadmb"){
	eval(parse(text=paste("if(!require(",pkg,",quietly=TRUE)) stop(\"",pkg," package is required\")",sep="")))
	tdir <- tempdir()
	tdir <- gsub("\\\\","/",tdir)                      # temporary directory for R
	pkgO=ls(paste("package:",pkg,sep=""),all.names=TRUE)                        # package objects
	z=sapply(pkgO,function(x){f=get(x);is.function(f)}); pkgF=names(z)[z] # package functions
	bad=regexpr("[\\|()[{^$*+?<-]",pkgF)
	pkgF=pkgF[bad<0]                # get rid of weird names
	if (length(pkgF)==0) {
		showAlert(paste(pkg,"has no recognizable functions"))
		return()
	}
	dots=regexpr("^\\.",pkgF)
	pkgF0=pkgF[dots==1]
	pkgF1=pkgF[dots!=1]
	code=c(paste("#",pkg,"Functions"),paste("#",paste(rep("-",nchar(pkg)+10),collapse="")))
	for (i in c(pkgF1,pkgF0)) {
		expr=paste("fun=deparse(",pkg,"::",i,"); fun[1]=paste(\"",i,"\",fun[1],sep=\" = \",collapse=\"\"); code=c(code,fun)",sep="")
		eval(parse(text=expr)) 
	}
	fname=paste(tdir,"/",pkg,".r",sep="")
	writeLines(code, fname)
	editADfile(fname)
	invisible()
}
#------------------------------------.win.viewCode

# functions called from window description files
#.win.onClose  = function(){ atcall(.onClose)() }
#.win.runModHelperQuit = function(){ atcall(.runModHelperQuit)() }

