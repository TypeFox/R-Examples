##'Compile and/or run an ADMB model, collect output
##'
##'Compile an ADMB model, run it, collect output
##'
##'\code{do_admb} will attempt to do everything required to start from the model
##'definition (TPL file) specified by \code{fn}, the data list, and the list of
##'input parameters, compile and run (i.e. minimize the objective function of)
##'the model in AD Model Builder, and read the results back into an object of
##'class \code{admb} in R. If \code{checkparam} or \code{checkdata} are set to
##'"write", it will attempt to construct a DATA section, and construct or
##'(augment an existing) PARAMETER section (which may contain definitions of
##'non-input parameters to be used in the model). It copies the input TPL file
##'to a backup (.bak); on finishing, it restores the original TPL file and
##'leaves the auto-generated TPL file in a file called [fn]_gen.tpl.
##'
##'@param fn (character) base name of a TPL function, located in the working
##'directory
##'@param data a named list of input data variables (order must match TPL file): each element of the list can either be a single value, or a list containing elements
##'\itemize{
##' \item{value}{the value of the data}
##' \item{data_type}{character: possible values as in \code{\link{storage.mode}}, typically "integer" or "numeric": this overrides R2admb's attempts to guess whether variables are supposed to be integers or floats (default NA)}
##' }
##'@param params a named list of starting parameter values (order must match TPL file): each element of the list can either be a single value, or a list containing elements
##'\describe{
##' \item{value}{starting value of the parameter (default 0)}
##' \item{bounds}{two-element vector of lower and upper bounds}
##' \item{phase}{integer, specifying phase: not implemented yet}
##' }
##'@param bounds named list of 2-element vectors of lower and upper bounds for
##'specified parameters
##'@param phase named numeric vector of phases (not implemented yet)
##'@param re a named list of the identities and dimensions of any random effects
##'vectors or matrices used in the TPL file
##'@param data_type a named vector specifying (optional) data types for
##'parameters, in parname="storage mode" format (e.g.
##'\code{c(x="integer",y="numeric")})
##'@param safe (logical) compile in safe mode?
##'@param profile (logical) generate likelihood profiles? (untested!)
##'@param profile.opts (list) list of options, including
##' \itemize{
##' \item{pars}{vector of names of parameters to profile}
##' }
##'@param mcmc (logical) run MCMC around best fit?
##'@param mcmc.opts options for MCMC (see \code{\link{mcmc.control}} for details)
##'@param impsamp (logical) run importance sampling?
##'@param verbose (logical) print details
##'@param run.opts options for ADMB run (see \code{\link{run.control}} for details)
##'@param objfunname (character) name for objective function in TPL file (only relevant if \code{checkparam} is set to "write")
##'@param workdir temporary working directory (dat/pin/tpl files will be copied)
##'@param admb_errors how to treat ADMB errors (in either compilation or run): use "ignore" option at your own risk!
##'@param extra.args (character) extra argument string to pass to admb
##'@return An object of class \code{admb}.
##'@note 1. Mixed-case file names are ignored by ADMB; this function makes a
##'temporary copy with the file name translated to lower case. 2. Parameter
##'names containing periods/full stops will not work, because this violates C
##'syntax (currently not checked). 3. There are many, many, implicit
##'restrictions and assumptions: for example, all vectors and matrices are
##'assumed to be indexed starting from 1.
##' @export
##'@author Ben Bolker
##'@keywords misc
##'@examples
##'
##'\dontrun{
##'setup_admb()
##'file.copy(system.file("tplfiles","ReedfrogSizepred0.tpl",package="R2admb"),"tadpole.tpl")
##' tadpoledat <-
##'  data.frame(TBL = rep(c(9,12,21,25,37),each=3),
##'             Kill = c(0,2,1,3,4,5,0,0,0,0,1,0,0,0,0L),
##'             nexposed=rep(10,15))
##'m1 <- do_admb("tadpole",
##'              data=c(list(nobs=15),tadpoledat),
##'              params=list(c=0.45,d=13,g=1),
##'              bounds=list(c=c(0,1),d=c(0,50),g=c(-1,25)),
##'              run.opts=run.control(checkparam="write",
##'                checkdata="write",clean="all"))
##'m2 <- do_admb("tadpole",
##'              data=c(list(nobs=15),tadpoledat),
##'              params=list(c=list(0.45,bounds=c(0,1)),
##'                          d=list(13,bounds=c(0,50)),
##'                          g=list(1,bounds=c(-1,25))),
##'              run.opts=run.control(checkparam="write",
##'                checkdata="write",clean="all"))
##'unlink("tadpole.tpl")
##'}
##'
do_admb <- function(fn,
                    data=NULL,
                    params=NULL,
                    bounds=NULL,
                    phase=NULL,
                    re=NULL,
                    data_type=NULL,
                    safe=TRUE,
                    profile=NULL,
                    profile.opts=NULL,
                    mcmc=NULL,
                    mcmc.opts=mcmc.control(),
                    impsamp=FALSE,
                    verbose=FALSE,
                    run.opts=run.control(),
                    objfunname="f",
                    workdir=getwd(),
                    admb_errors=c("stop","warn","ignore"),
                    extra.args) {
    ## TO DO: check to see if executables are found
    ## MODULARIZE (separate sub-function):
    ##  1. check or construct input&data files, optionally write TPL file
    ##  2. compile (tpl -> rem/cpp -> binary)  [DONE]
    ##  3. run
    ##  4. retrieve & package output
    admb_errors <- match.arg(admb_errors)
    if (!missing(workdir)) {
        ## copy files to working directory
        ## FIXME: check on case-sensitivity??
        file.copy(list.files(pattern=paste(fn,"\\.(dat|pin|tpl)",sep="")),
                  workdir)
        cwd <- setwd(workdir)
        on.exit(setwd(cwd))
    }
    if (is.null(data) && is.null(params) && missing(run.opts)) {
        ## as far as we can tell there is a fully specified ADMB model,
        ## the user just wants to run it and retrieve the results -- not
        ## do any DAT/PAR/TPL construction
        checkparam <- checkdata <- "ignore"
    } else {
        checkparam <- run.opts["checkparam"]
        checkdata <- run.opts["checkdata"]
    }
    if (is.null(mcmc)) mcmc <- !missing(mcmc.opts)
    if (mcmc) {
        if (checkparam=="write" && !"mcmcpars" %in% names(mcmc.opts))
            stop("must specify mcmcpars when checkparam=='write' and mcmc is TRUE")
    }
    if (is.null(profile)) profile <- !missing(profile.opts)
    if (profile & !is.null(re)) {
        stop("profiling is not implemented for models with random effects")
    }
    if (profile && is.null(profile.opts$pars) && checkparam=="write")
        stop("must specify profpars when checkparam=='write' and profile is TRUE")
    if (is.null(re)) {
        if (mcmc && "mcmc2" %in% names(mcmc.opts)) stop("mcmc2 can only be used with random-effects models")
    }
    tplfile <- paste(fn,"tpl",sep=".")
    if (!file.exists(tplfile)) stop("can't find TPL file ",tplfile)
    if (run.opts["check_tpl"]) {
        tpldat <- read_tpl(fn)  ## extract info from TPL file
        tplinfo <- tpldat$info
        orig_tplfile <- tplfile
        ofn <- fn
        ## FIXME: need to make this work on MacOS/Windows, and safe
        if (test_OScase() && !tolower(fn)==fn) {
            warning("Base name converted to lower case for ADMB compatibility: copying TPL file")
            tplfile <- tolower(tplfile)
            fn <- tolower(fn)
            ## if (file.exists(tplfile)) stop("refusing to write over existing (lowercase) TPL file")
            file.copy(orig_tplfile,tplfile,overwrite=TRUE)
        }
        ## check PARAMETER section
        getvals <- function(L,el="bounds",getFirst=FALSE,valsOK) {
            res <- lapply(L,function(x) {
                if (is.list(x)) {
                    if (is.null(x[[el]]) && getFirst) {
                        x[[1]]  ## unnamed first arg
                    } else x[[el]] ## return arg or NULL
                } else {
                    if (valsOK) x else NULL  ## non-list
                }
            })
            res <- res[!sapply(res,is.null)]
            if (length(res)==0) res <- NULL
            res
        }
        if (!checkparam %in% c("write","ignore") && is.null(tplinfo$inits))
            stop("must specify PARAMETER section (or set 'checkparam' to 'write' or 'ignore')")
        if (checkparam!="ignore") {
            if (missing(bounds)) bounds <- getvals(params,"bounds",valsOK=FALSE)
            if (missing(phase)) phase <- getvals(params,"phase",valsOK=FALSE)
            if (any(unlist(lapply(params,names))=="re")) {
                ## replace TRUE with vector length
                params <- lapply(params,
                                 function(x) {
                                     if ("re" %in% names(x)) {
                                         len <- if ("value" %in% names(x)) {
                                             length(x[["value"]])
                                         } else if (names(x)[1]=="") {
                                             length(x[[1]])
                                         } else stop("first element must be unnamed or 'value' must be specified")
                                         x[["re"]] <- len
                                     }
                                     x
                                 })
            }
            if (missing(re)) re <- getvals(params,"re",valsOK=FALSE)
            params <- getvals(params,"value",getFirst=TRUE,valsOK=TRUE)
            dmsg <- check_section(ofn,tpldat,"inits",params,
                                  check=checkparam,
                                  bounds=bounds,
                                  data_type=data_type,
                                  phase=phase,
                                  secname="PARAMETER",
                                  objfunname=objfunname,
                                  re=re,
                                  mcmcpars=mcmc.opts$mcmcpars,
                                  profpars=profile.opts$pars)
        }
        if (!checkparam %in% c("write","ignore") && nchar(dmsg)>0) {
            if (checkparam=="stop") stop(dmsg)
            if (checkparam=="warn") warning(dmsg)
        } else if (checkparam=="write") {
            ## FIXME: break this out to a separate write_tpl function
            if (!is.null(tpldat$secs$PARAMETER)) {
                tpldat$secs$PARAMETER <- dmsg
            } else {
                ## insert immediately before PROCEDURE
                tpldat$secs <- append(tpldat$secs,list(PARAMETER=c(dmsg,"")),
                                      after=which(names(tpldat$secs)=="PROCEDURE")-1)
            }
            ## modifications to PROCEDURE section:
            ## need to assign MCMC reporting variables
            if (mcmc) {
                mcmcparnames <- gsub("^ +sdreport_(number|vector) r_","",
                                     gsub("\\(.*$","",
                                          dmsg[grep("^ +sdreport",dmsg)]))
                if (length(mcmcparnames)>0) {
                    tpldat$secs$PROCEDURE <- append(tpldat$secs$PROCEDURE,
                                                    indent(paste("r_",mcmcparnames,"=",mcmcparnames,";",sep="")))
                }
            }
            if (profile) {
                profparnames <- gsub("^ +likeprof_number p_","",
                                     gsub("\\(.*$","",
                                          dmsg[grep("^ +likeprof_",dmsg)]))
                tpldat$secs$PROCEDURE <- append(tpldat$secs$PROCEDURE,
						indent(paste("p_",profparnames,"=",profparnames,";",sep="")))
            }
        }
        ## check DATA section
        if (checkdata!="ignore") {
            dframes <- sapply(data,data.class)=="data.frame"
            if (any(dframes)) warning("attempted to convert data frame to matrix")
            data[dframes] <- lapply(data[dframes],as.matrix)
        }                              
        if (!checkdata %in% c("write","ignore") && is.null(tplinfo$data))
            stop("must specify DATA section (or set 'checkdata' to 'write' or 'ignore')")
        if (checkdata != "ignore") {
            dmsg <- check_section(ofn,tpldat,"data",data,
                                  check=checkdata,
                                  data_type=data_type,
                                  secname="DATA")
            if (!checkdata %in% c("write","ignore") && nchar(dmsg)>0) {
                if (checkdata=="stop") stop(dmsg)
                if (checkdata=="warn") warning(dmsg)
                
            } else if (checkdata=="write") {
                if (!is.null(tpldat$secs$DATA)) {
                    tpldat$secs$DATA <- dmsg
                } else {
                    ## insert immediately before PARAMETER
                    tpldat$secs <- append(tpldat$secs,list(DATA=c(dmsg,"")),
                                          after=which(names(tpldat$secs)=="PARAMETER")-1)
                }
            }
        }
        ##
        if (checkdata=="write" || checkparam=="write") {
            parnames <- c(names(data),names(params))
            badnames <- grep("\\.",parnames)
            if (length(badnames)>0) {
                for (i in badnames) {
                    old <- parnames[i]
                    new <- gsub("\\.","_",parnames[i])
                    tpldat$secs$PROCEDURE <- gsub(old,new,tpldat$secs$PROCEDURE)
                }
            }
            ## fn2 <- paste(fn,".tpl",sep="")
            fn <- paste(fn,"_gen",sep="")
            tplfile <- paste(fn,".tpl",sep="")
            ## fn2bak<- paste(fn,".tpl.bak",sep="")
            ## file.copy(fn2,fn2bak,overwrite=TRUE)
            ## FIXME: work with generated file, don't touch original!
            ## on exit, copy auto-generated file and restore original ...
            ## on.exit(file.copy(fn2,fn2gen,overwrite=TRUE),add=TRUE)
            ## on.exit(file.copy(fn2bak,fn2,overwrite=TRUE),add=TRUE)
            writeLines(do.call("c",tpldat$secs),con=tplfile)
            tpldat <- read_tpl(fn) ## get auto-generated info
        }
    }
    if (run.opts["write_files"]) {
        if (verbose) cat("writing data and parameter files ...\n")
        ## check order of data; length of vectors???
        dat_write(fn,data)
        ## check order of parameters ??
        ## add random effects to list of initial parameters
        if (!is.null(re)) {
            rv <- re[!names(re) %in% names(params)]
            params <- c(params,lapply(as.list(rv),rep,x=0))
        }
        pin_write(fn,params)
    }
    ## insert check(s) for failure at this point
    ## PART 2A: compile
    if (run.opts["compile"]) {
        compile_admb(fn,safe,re=!is.null(re),verbose)
    }
    ## PART 2B: run executable file
    if (run.opts["run"]) {
        res <- run_admb(fn,verbose,mcmc,mcmc.opts,profile,extra.args,admb_errors)
    } else res <- NULL
    if (run.opts["read_files"]) {
        ## pnames <- get_names(c(mcmc.opts[["mcmcpars"]],names(re)),
        ## tpldat$info)
        ## PART 3
        ## FIXME: we should be able to recover info about prof and MCMC par names from
        ##  somewhere ... re-read TPL files etc.?
        L <- read_admb(fn,verbose,
                       profile,
                       mcmc,
                       admbOut=res)
    }

    ## need to do this **AFTER** auto-generated versions get swapped around
    ## cover both cases
    ## FIXME: check -- why is this here? it's "on.exit" -- can't it go earlier?
    on.exit(clean_admb(fn,run.opts["clean_files"]),add=TRUE)
    ## check for NA/NaN in logLik, errors in text?
    if (run.opts["read_files"]) L else NULL
}
test_OScase <- function(dir=getwd()) {
    fn1 <- tempfile(tmpdir=dir)
    fn2 <- gsub("/([^/]+)$","/\\U\\1",fn1,perl=TRUE)
    if (!file.create(fn1)) stop("can't create temporary file")
    ## attempt to copy: will return FALSE if 'file already exists', i.e.
    ##  file system is case-insensitive.  Could also use file.rename(), but
    ##  would have to suppress warning msg
    res <- file.copy(fn1,fn2)
    unlink(fn1)
    unlink(fn2)
    res
}
