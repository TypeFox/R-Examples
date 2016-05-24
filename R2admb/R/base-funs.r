#'Compile ADMB files, run, read output
#'
#'With various tests, calls the \code{admb} script to compile from a TPL file
#'to an executable, or runs the resulting executable
#'@usage compile_admb(fn,safe=FALSE,re=FALSE,
#' verbose=FALSE,
#' admb_errors=c("stop","warn","ignore"))
#' 
#'  run_admb(fn,verbose=FALSE,mcmc=FALSE,
#' mcmc.opts=mcmc.control(),profile=FALSE,
#' extra.args="",admb_errors=c("stop","warn","ignore"))
#' 
#' read_admb(fn,verbose=FALSE,profile=FALSE,
#' mcmc=FALSE,mcmc.opts=NULL,admbOut=NULL,checkterm=TRUE)
#'@aliases compile_admb run_admb read_admb
#'@export compile_admb run_admb read_admb
#'@param fn (character) name of TPL file, without extension
#'@param safe (logical) Compile in safe mode?
#'@param re (logical) Compile in random effects (ADMB-RE) mode?
#'@param profile (logical) Run likelihood profiles?
#'@param extra.args (character) extra arguments for ADMB run
#'@param mcmc (logical) run post-hoc MCMC?
#'@param mcmc.opts options for MCMC run (see \code{\link{mcmc.control}})
#'@param verbose (logical) Verbose output?
#'@param admb_errors (character) how to handle compilation/linking errors?
#'@param admbOut (character) ADMB run output for inclusion in \code{admb}
#'object (for internal use)
#'@param checkterm (logical) compute termination criteria (ratio of min/max
#'eigenvalue) and include it in the saved object?
#'@return \itemize{ \item \code{compile_admb} returns nothing (it has the side
#'effect of creating an executable) \item \code{run_admb} invisibly returns the
#'output produced by the ADMB run; it also produces output files on disk as a
#'side effect \item \code{read_admb} returns an object of class \code{admb},
#'containing as much information as possible gleaned from the output files
#'(parameter estimates, standard errors, variance-covariance matrix, profiles,
#'MCMC output) }
#'@note Compiling also sets executable mode.
#'@author Ben Bolker
#'@keywords misc
compile_admb <- function(fn,safe=FALSE,re=FALSE,verbose=FALSE,
		admb_errors=c("stop","warn","ignore")) {
	admb_errors <- match.arg(admb_errors)
	if (!file.exists(fn2 <- paste(fn,"tpl",sep=".")))
            stop("can't find TPL file ",fn2)
	test <- try(system("admb",intern=TRUE),silent=TRUE)
	if (inherits(test,"try-error")) stop("base admb command failed: run setup_admb(), or check ADMB installation")
	args <- ""
	if (re) args <- "-r"
	if (safe) args <- paste(args,"-s")
	if (verbose) cat("compiling with args: '",args,"' ...\n")
	res0 <- system(paste("admb",args,fn," 2>",paste(fn,".cout",sep="")),
			intern=TRUE)
	coutfile <- readLines(paste(fn,".cout",sep=""))
	if (verbose) {
		cat("compile output:\n",res0,"\n")
		cat("compile log:\n")
		cat(coutfile,sep="\n")
	}
	## sorting out the lines that come BEFORE the warnings
	admb_warnerr_index <- grepl("warning|error",coutfile)
	csplit <- split(coutfile,head(c(0,cumsum(admb_warnerr_index)),-1))
	wchunks <- which(sapply(lapply(csplit,grep,pattern="warning"),length)>0)
	echunks <- which(sapply(lapply(csplit,grep,pattern="error"),length)>0)
	if (length(wchunks)>0) {
		if (!verbose) {
			## figure we don't need these warnings
			## if we are spitting them out above anyway
			admb_warnings <- paste("from ADMB:",unlist(csplit[wchunks]))
			sapply(admb_warnings,warning)
		}
		csplit <- csplit[-wchunks]
	}
	Sys.chmod(fn,mode="0755")
	if (length(echunks)>0) {
		comperrmsg <- "errors detected in compilation: run with verbose=TRUE to view"
		if (admb_errors=="stop") stop(comperrmsg) else if (admb_errors=="warn") warning(comperrmsg)
	}
}
run_admb <- function(fn,verbose=FALSE,mcmc=FALSE,mcmc.opts=mcmc.control(),
		profile=FALSE,extra.args="",
		admb_errors=c("stop","warn","ignore")) {
	admb_errors <- match.arg(admb_errors)
	args <- ""
	if (mcmc) {
		if (is.null(mcmc.opts$mcmcpars)) stop("you must specify at least one parameter in 'mcmc.opts$mcmcpars' (see ?mcmc.control)")
		args <- paste(args,mcmc.args(mcmc.opts))
	}
	if (profile) args <- paste(args,"-lprof")
	if (!missing(extra.args)) {
		args <- paste(args,extra.args)
	}
	if (verbose) cat("running compiled executable with args: '",args,"'...\n")
	
	outfn <- paste(fn,"out",sep=".")
	
	if (.Platform$OS.type=="windows") {
		cmdname <- paste(fn,".exe",sep="")
		shellcmd <- shell
	} else {
		cmdname <- paste("./",fn,sep="")
		shellcmd <- system
	}
	if (!file.exists(cmdname)) stop("executable ",cmdname," not found: did you forget to compile it?")
	res <- shellcmd(paste(cmdname,args,">",outfn),intern=TRUE)
	
	outfile <- readLines(paste(fn,".out",sep=""))
	## replace empty res with <empty> ?
	if (mcmc) {
		## write MC info to disk so it will be retrievable ...
		mcinfofile <- file(paste(fn,"mcinfo",sep="."),"w")
		mctab <- unlist(mapply(function(x,y) {
							c(paste("# ",x),if (is.null(y))  "" else paste(y,collapse=" "))
						},names(mcmc.opts),mcmc.opts))
		writeLines(mctab,paste(fn,"mcinfo",sep="."))
	}
	if (verbose) {
		cat("Run output:\n",res,"\n",sep="\n")
		cat(outfile,"\n",sep="\n")
	}
	if (length(grep("^Error",outfile)>0)) {
		runerrmsg <- "errors detected in ADMB run: run with verbose=TRUE to view"
		if (admb_errors=="stop") stop(runerrmsg) else if (admb_errors=="warn") warning(runerrmsg)
	}
	invisible(res)
}
read_admb <- function(fn,verbose=FALSE,
		profile=FALSE,
		mcmc=FALSE,
		mcmc.opts=NULL,
		admbOut=NULL,
		checkterm=TRUE) {
	tpldat <- read_tpl(fn)  ## extract info from TPL file
	if (verbose) cat("reading output ...\n")
	parfn <- paste(fn,"par",sep=".")
	if (!file.exists(parfn)) stop("couldn't find parameter file ",parfn)
	L <- c(list(fn=fn,txt=admbOut),read_pars(fn))
	if (mcmc) {
		## if (checkparam!="write") {
		## warning("MCMC naming is probably wrong")
		## }
		## FIXME: get MCMC names -- how?
		if (file.exists(paste(fn,"hst",sep="."))) {
			L <- c(L,list(hist=read_hst(fn)))
		}
		if (is.null(mcmc.opts)) {
			## try to retrieve mc info from file
			mcinfofile <- paste(fn,"mcinfo",sep=".")
			if (file.exists(mcinfofile)) {
				w <- readLines(mcinfofile)
				wnames <- gsub("^# +","",w[seq(1,length(w),by=2)])
				wvals <- as.list(w[seq(2,length(w),by=2)])
				wvals[c(1,2,5)] <- as.numeric(wvals[c(1,2,5)])
				wvals[3:4] <- as.logical(wvals[3:4])
				wvals[[6]] <- strsplit(wvals[[6]]," ")[[1]]
				names(wvals) <- wnames
				mcmc.opts <- wvals
			} else warning("having difficulty retrieving MCMC info, will try to continue anyway")
		}
		## if (is.null(mcmc.opts) || is.null(mcmc.opts$mcmcpars) || nchar(mcmc.opts$mcmcpars)==0) {
		## pnames <- gsub("r_","",tpldat$info$sdnums$vname)
		## } else pnames <- mcmc.opts$mcmcpars
		sdinfo <- read.table(paste(fn,"std",sep="."),skip=1)
		pnames <- rep_pars(sdinfo[,2])
		## pnames <- grep("^r_.*",pnames,value=TRUE,invert=TRUE)
		sdreportvars <- as.character(tpldat$info$sdnums$vname)
		pnames <- setdiff(pnames,sdreportvars)
		if (is.null(mcmc.opts) || mcmc.opts[["mcsave"]]>0) {
                     ## exclude autogenerated sdreport variables
                        mcpnames <- pnames[!grepl("^r_",pnames)]
			L$mcmc <- read_psv(fn,names=mcpnames)
			## FIXME: account for mcmc2 if appropriate
			attr(L$mcmc,"mcpar") <- c(1,mcmc.opts[["mcmc"]],mcmc.opts[["mcsave"]])
		}
	}
	if (profile) {
		if (!is.null(tpldat$info$raneff)) {
			stop("something's wrong -- profiling is not implemented for random effects models")
		}
		profpars <- tpldat$info$profparms$vname
		L$prof <- lapply(profpars,read_plt)
		names(L$prof) <- gsub("p_","",profpars)  ## FIXME: maybe dangerous?
	}
	## compute termination criteria
	##  can we retrieve hessian directly???
	if (checkterm) {
		v <- with(L,vcov[seq(npar),seq(npar)])
		ev <- try(eigen(solve(v))$value,silent=TRUE)
		L$eratio <- if (inherits(ev,"try-error")) NA else min(ev)/max(ev)
	}
	class(L) <- "admb"
	L
}


