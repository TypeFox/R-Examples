##'Read in parameters from an AD Model Builder run
##'
##'Reads coefficients, standard errors, log-likelihoods, maximum gradients,
##'correlation and variance-covariance matrices from AD Model Builder output
##'files
##'
##'Given the output from an ADMB run on FOO.tpl, \code{read_pars} reads the
##'files FOO.par (parameters, log-likelihood, max gradient); FOO.std (standard
##'deviations); FOO.cor (correlations); FOO.rep (report variables);
##'admodel.hes for hessian; and admodel.cov for
##'covariance matrix.  \code{read_psv} reads the output of MCMC runs
##'
##'@rdname read_pars
##'@param fn (character) Base name of AD Model Builder
##'@param drop_phase (logical) drop negative-phase (fixed) parameters from results?
##'@return List containing the following elements
##' \itemize{
##' \item{coefficients}{parameter estimates}
##' \item{coeflist}{parameter estimates in list format, with proper shape (vectors, matrices, etc.)}
##' \item{se}{estimated standard errors of coefficients}
##' \item{loglik}{log-likelihood}
##' \item{maxgrad}{maximum gradient of log-likelihood surface}
##' \item{cor}{correlation matrix}
##' \item{vcov}{variance-covariance matrix}
##' \item{npar}{number of parameters}
##' \item{hes}{hessian matrix (only if no vcov matrix)}
##' }
##'@section Warning: The \code{coeflist} component is untested for data
##'structures more complicated than scalars, vectors or matrices (i.e. higher-dimensional or ragged arrays)
##'@author Ben Bolker
##'@seealso \code{\link{write_pin}}, \code{\link{write_dat}}
##'@keywords misc
##'@export read_pars

read_pars <- function (fn,drop_phase=TRUE) {
    ## see
    ##  http://admb-project.org/community/admb-meeting-march-29-31/InterfacingADMBwithR.pdf
    ## for an alternate file reader -- does this have equivalent functionality?
    ## FIXME: get hessian.bin ?
    ## get parameter estimates
    par_dat <- rs(fn,"par", skip = 1)
    tmp <- rs(fn, "par", what = "", comment.char="")
    ## COULD get parnames out of par file, but a big nuisance
    ##  for vectors etc.
    ## parnames <- gsub(":$","",tmp[seq(18,by=3,length=npar)])
    loglik <- as.numeric(tmp[11])
    maxgrad <- as.numeric(tmp[16])
    ## second pass to extract names from par file (ugh)
    pars <- read_pars0(paste(fn,"par",sep="."),par_dat)
    ## FIXME: refactor
    est <- pars$est
    parlist <- pars$parlist
    parnames0 <- pars$parnames0
    parnameslist <- pars$parnameslist
    if (!is.finite(loglik))
        warning("bad log-likelihood: fitting problem in ADMB?")
    nopar <- NULL ## in case nopar wasn't read from vcov
    ## if non-pos-def hessian, cor and std files will be missing ... but
    ##   we should still be able to retrieve some info
    sd_dat <- rt(fn,"std", skip = 1,as.is=TRUE)
    if (length(sd_dat)==1 && is.na(sd_dat)) {
        warning("std file missing: some problem with fit, but retrieving parameter estimates anyway")
        npar <- length(est)
        cormat <- vcov <- matrix(NA,nrow=npar,ncol=npar)
        std <- rep(NA,npar)
        sdrptvals <- numeric(0)
        sdparnames <- NULL
    } else {
        nsdpar <- nrow(sd_dat)
        ## need col.names hack so read.table knows how many
        ##  columns to read: ?read.table, "Details"
        ## FIXME: go gracefully if .cor missing?
        if (file.exists(paste(fn,"cor",sep="."))) {
            ncorpar <- length(readLines(paste(fn,"cor",sep=".")))-2
            cor_dat <- rt(fn,"cor", skip = 2, fill=TRUE, 
                          as.is=TRUE,col.names=paste("X",1:(4+ncorpar),sep=""))
            ## drop cors that are not parameters
            ## (have dropped mc parameters)
            cormat <- as.matrix(cor_dat[1:nsdpar,4+(1:nsdpar)])
            cormat[upper.tri(cormat)] <- t(cormat)[upper.tri(cormat)]
            ## be careful here -- need to adjust for phase<0 parameters,
            ##  which will be in parameter vector but not in
        } else {
            warning(".cor file not found")
            cormat <- NA
        }
        ##  sd
        sdparnames <- sd_dat[, 2]
        misspars <- setdiff(parnames0,sdparnames)
        ## only names of positive-phase parameters
        parnames2 <- unlist(parnameslist[!parnames0 %in% misspars])
        sdparnames <- c(parnames2,sdparnames[-seq_along(parnames2)])
        ## parnames <- c(parnames,sd_dat[-seq_along(parnames),2])
        if (any(duplicated(sdparnames))) {
            sdparnames <- rep_pars(sdparnames)
        }
        npar3 <- length(parnames2) ## positive-phase only
        if (drop_phase) {
            parlist <- parlist[!parnames0 %in% misspars]
            est <- unlist(parlist)
            names(est) <- parnames2
            npar <- npar3
        }
        std <- sd_dat[, 4]
        sdrptvals <- sd_dat[-(1:npar3),3]
        sdrptnames <- sd_dat[-(1:npar3),2]
        names(sdrptvals) <- sdrptnames
        ## less accurate, but available for all parameters (RE, extra)
        vcov <- cormat*outer(std,std)
        if (!is.matrix(cormat)) cormat <- vcov ## (missing cor file)
        if(file.exists("admodel.cov"))
        {
            ## see also:
            ## http://www.admb-project.org/examples/admb-tricks/covariance-calculations/r-code-to-interact-with-admodel.cov/view
            ## more accurate, binary info
            filen <- file("admodel.cov", "rb")
            nopar <- readBin(filen, what = "integer", n = 1)
            vcov2 <- readBin(filen, what = "double", n = nopar * nopar)
            hybrid_bounded_flag <- readBin(filen, "integer", 1)
            scale <- readBin(filen, "numeric", nopar)
            vcov2 <- matrix(vcov2, byrow = TRUE, ncol = nopar) *
                (scale %o% scale)
            close(filen)
            ## replace vcov[] with better values
            vcov[1:nopar,1:nopar] <- vcov2
        }
    }
    ##  check if sdparnames is null; jll modified to use ncol(vcov) for vcov
    if(!is.null(sdparnames)) {
        names(std) <- rownames(cormat) <- colnames(cormat) <- sdparnames
        colnames(vcov) <- rownames(vcov) <- sdparnames[seq(ncol(vcov))]
    }
    ##  JLL - file is padded at the end which is why read_admbbin doesn't
    ##   work; need to know the number of records.
    ##   For the hessian, nopar is number of rows and columns
    hes <- NULL
    if(all(is.na(vcov)) & file.exists("admodel.hes"))
    {
        filen <- file("admodel.hes", "rb")
        nopar <- readBin(filen, what = "integer", n = 1)
        hes <- readBin(filen, what = "double", n = nopar * nopar)
        hes <- matrix(hes, byrow = TRUE, ncol = nopar)
        colnames(hes) <- rownames(hes) <- sdparnames[seq(ncol(vcov))]
        close(filen)
    }
    reportfile <- paste(fn,"rep",sep=".")
    repvals <- list()
    if (file.exists(reportfile)) {
        if (file.info(reportfile)$size > 0){
            repvals <- read_rep(fn)
        }
    }
    npar_rep <- length(repvals$est)
    ## return npar as number of columns in vcov and npar_total as npar+re
    if (is.null(nopar)) nopar <- npar
    list(coefficients=c(est,sdrptvals,repvals$est),
         coeflist=c(parlist,repvals$parlist),
         se=std, loglik=-loglik, maxgrad=-maxgrad, cor=cormat, vcov=vcov,
         npar=nopar,           ## fixed-effect
         npar_re=npar-nopar,   ## random-effect
         npar_sdrpt=length(sdrptvals), ## sdreport
         npar_rep=npar_rep,     ## report
         npar_total=length(est)+length(sdrptvals)+npar_rep, ## all
         hes=hes)
}

rt <- function(f,ext,...) {
    fn <- paste(f,ext,sep=".")
    if (file.exists(fn)) read.table(fn,...) else NA
}

rs <- function(f,ext,comment.char="#",...) {
    fn <- paste(f,ext,sep=".")
    if (file.exists(fn)) scan(fn,
                              comment.char=comment.char,
                              quiet=TRUE,...) else NA
}


## FIXME: can I refactor slightly so I don't need to pass par_dat, or read twice?
read_pars0 <- function(fnext,par_dat,skip=1) {
    tmp2 <- readLines(fnext)
    if (skip>0) tmp2 <- tmp2[-(1:skip)]
    parlines <- grep("^#",tmp2)
    npar <- length(par_dat)      ## TOTAL number of numeric values
    npar2 <- length(parlines)  ## number of distinct parameters
    parlen <- count.fields(fnext,skip=skip)
    parlen2 <- count.fields(fnext,comment.char="",skip=skip)
    parnames0 <- parnames <- gsub("^# +","",gsub(":$","",tmp2[parlines]))
    parlist <- vector("list",npar2)
    parnameslist <- vector("list",npar2)
    names(parlist) <- parnames
    cumpar <- 1
    cumline <- 1
    pp <- c(parlines,length(tmp2)+1)
    ## reshape parameters properly
    parid <- numeric(npar2)
	## If report file didn't contain #value to identify report segments, 
	## it would cause:  Error in cumline:(cumline + nrows - 1) : NA/NaN argument; adding the if below 
	## solved that problem. J.Laake 7/1/2013
	if(npar2>0)
    for (i in seq(npar2)) {
        nrows <- diff(pp)[i]-1
        curlines <- cumline:(cumline+nrows-1)
        curlen <- sum(parlen[curlines])
        parvals <- par_dat[cumpar:(cumpar+curlen-1)]
                                        #    if (nrows==1) {
        if (curlen==1) {
            parnameslist[[i]] <- parnames[i]
        } else {
            parnameslist[[i]] <- numfmt(parnames[i],curlen)
        }
        parlist[[i]] <- parvals
                                        #    } else {
                                        #      parlist[[i]] <- matrix(parvals,nrow=nrows,byrow=TRUE)
                                        #      parnameslist[[i]] <- numfmt2(parnames[i],dim(parlist[[i]]))
                                        #    }
        ## cat(parnames[i],cumline,cumline+nrows-1,cumpar,cumpar+curlen-1,parnameslist[[i]],"\n")
        ## print(parlist[[i]])
        cumline <- cumline + nrows
        cumpar <- cumpar + curlen
    }
    ## FIXME: watch out, 'short param names' has now been overwritten by 'long param names'
    parnames <- unlist(parnameslist)
    
    ## parnames <- unname(unlist(mapply(function(x,len) {
    ## if (len==1) x else numfmt(x,len)  ## paste(x,1:len,sep=".")
    ##},
    ## parnames,parlen)))
    est <- unlist(par_dat)
    names(est) <- parnames
    return(list(est=est,parlist=parlist,parnames0=parnames0,parnameslist=parnameslist))
}

##'@rdname read_pars
##'@export read_psv
##'@param names (character) Names of variables
read_psv <- function(fn,names=NULL) {
    fn <- tolower(fn) ## arghv
    fn <- paste(fn,"psv",sep=".")
    if (!file.exists(fn)) stop("no PSV file found")
    ans <- read_admbbin(fn)
    if (is.null(names)) names <- paste("V",seq(ncol(ans)),sep="")
    if (length(names)<ncol(ans)) {
        warning("more columns than names: generating dummy names")
        names <- c(names,paste("V",seq(length(names)+1,ncol(ans)),sep=""))
    } else if (length(names)<ncol(ans)) {
        warning("more names than columns: truncating names")
        names <- names[seq(ncol(ans))]
    }
    colnames(ans) <- names
    ans <- as.data.frame(ans)
    ans
}


##'@rdname read_pars
##'@export read_rep
read_rep <- function(fn,names=NULL) {
    fn <- tolower(fn)
    fnx <- paste(fn,"rep",sep=".")
    ## check format!
    rr <- readLines(fnx)
    ## FIXME: could test for numerics more carefully ...
    numLines <- grepl("^[0-9. e]",rr)
    commentLines <- grepl("^#",rr)
    if (!all(numLines | commentLines) || !commentLines[1]) {
        warning("report file in non-standard format")
        if (all(numLines)) {
            return(list(est=scan(text=rr,quiet=TRUE)))
        }
    }
    par_dat <- rs(fn,"rep")
    if (!file.exists(fnx)) stop("no REP file found")
    read_pars0(fnx,par_dat,skip=0)
}

read_tpl <- function(f) {
    r <- readLines(paste(f,"tpl",sep="."))
    ## drop blank lines at head
    if (!nzchar(r[1]))
        r <- r[-(1:rle(!nzchar(r))$lengths[1])]
    ## secStart <- which(substr(r,1,1) %in% LETTERS)
    secStart <- grep("^ *[A-Z]+(_[A-Z]+)+$",r)
    calcLines <- grep("_CALCS *$",r)
    secStart <- setdiff(secStart,calcLines)
    ## which(substr(r,1,1) %in% LETTERS)
    if (length(secStart)==0) stop("tpl file must contain at least one section (capitalized header)")
    if (secStart[1]!=1) { ## add first (comments etc.) section
        secStart <- c(1,secStart)
    }
    nsec <- length(secStart)
    ## length (in lines) of each chunk
    L <- c(secStart[-1],length(r)+1)-secStart
    sec <- rep(1:nsec,L)
    splsec <- split(r,sec)
    ## not QUITE right: we get some stuff in here that is not
    ##  a "SECTION" but is SEPARABLE_FUNCTION or TOP_OF_MAIN_CALCS
    ##  or something ...
    splnames <- sapply(splsec,"[",1)
                                        #	names(splsec) <- ifelse(grepl("SECTION$",splnames),
                                        #			gsub("_.+","",splnames),
                                        #			splnames)
    names(splsec) <- gsub("_.+","",splnames)
    splsec_proc <- lapply(splsec,drop_calcs)
    L1 <- L2 <- NULL
    pp <- splsec_proc$PARAMETER
    ## EXPERIMENTAL:
    pp <- pp[!grepl("^ *!!",pp)]
    if (!is.null(pp)) {
        pp <- proc_var(pp,maxlen=7)
        type <- 1 ## kluge for R CMD check warnings; will be masked
        if (!is.null(pp)) {
            L1 <- with(pp,
                       list(inits=pp[grep("^init",type),],
                            raneff=pp[grep("^random",type),],
                            sdnums=pp[grep("^sdreport_number",type),],
                            sdvecs=pp[grep("^sdreport_vector",type),],
                            other=pp[grep("^init|random|sdreport",type,invert=TRUE),],
                            ## FIXME: don't know what I needed this for
                            ## sdvecdims <- gsub("^ +sdreport_vector[ a-zA-Z]+","",
                            ## gsub("[()]","",
                            ## grep( "^ +sdreport_vector",splsec$PARAMETER,
                            ## value=TRUE)))
                            profparms=pp[grep("^likeprof",type),]))
        }
    }
    pp <- splsec_proc$DATA
    if (!is.null(pp)) {
        pp <- proc_var(pp,maxlen=7)
        L2 <- with(pp,
                   list(data=pp[grep("^init",type),]))
    }
    L <- c(L1,L2)
    L <- L[!sapply(L,is.null)]
    list(secs=splsec,info=L[sapply(L,nrow)>0])
}

##'Read in ADMB profile file
##'
##'Read in the output from ADMB likelihood profiling stored in a \code{.plt}
##'file
##'
##'@usage read_plt(varname)
##'@export
##'@param varname (character) Name of profiled variable (base name of \code{.plot} file)
##'@return List containing the following elements:
##' \item{prof}{likelihood profile: a two-column matrix containing the parameter value and the corresponding likelihood (\emph{not} the log-likelihood or negative log-likelihood), scaled to integrate to 1.0}
##' \item{ci}{matrix of upper and lower confidence intervals at the 0.9, 0.95, and 0.975 levels}
##' \item{prof_norm}{likelihood profile, based on a normal approximation}
##' \item{cinorm}{confidence interval matrix, based on normal approximation}

read_plt <- function(varname) {
    fn <- paste(varname,"plt",sep=".")
    r <- readLines(fn)
    cisecline <- grep("Minimum width confidence limits",r)
    normline <- grep("Normal approximation$",r)
    t1 <- textConnection(r[3:(cisecline[1]-1)])
    prof1 <- matrix(scan(t1,quiet=TRUE),ncol=2,
                    byrow=TRUE,
                    dimnames=list(NULL,c("value","likelihood")))
    close(t1)
    t2 <- textConnection(r[cisecline[1]+(2:4)])
    ci1 <- matrix(scan(t2,quiet=TRUE),ncol=3,
                  byrow=TRUE,
                  dimnames=list(NULL,c("sig","lower","upper")))
    close(t2)
    t3 <- textConnection(r[(normline+1):(cisecline[2]-1)])
    profnorm <- matrix(scan(t3,quiet=TRUE),ncol=2,
                       byrow=TRUE,
                       dimnames=list(NULL,c("value","logliK")))
    close(t3)
    t4 <- textConnection(r[cisecline[2]+(2:4)])
    cinorm <- matrix(scan(t4,quiet=TRUE),ncol=3,
                     byrow=TRUE,
                     dimnames=list(NULL,c("sig","lower","upper")))
    close(t4)
    list(prof=prof1,ci=ci1,prof_norm=profnorm,ci_norm=cinorm)
}
## read a "standard" ADMB format binary file into R:
##  standard format is: 1 integer describing number
##  of (double) values per vector
## FIXME: check bin sizes 
read_admbbin <- function(fn) {
    f <- file(fn,open="rb")
    nv <- readBin(f,"int")
    fs <- file.info(fn)$size
    isize <- 4; dsize <- 8
    m <- matrix(readBin(f,"double",n=(fs-isize)/dsize),byrow=TRUE,
                ncol=nv)
    close(f)
    m
}
