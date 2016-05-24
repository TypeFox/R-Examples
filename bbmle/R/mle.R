## require(methods,quietly=TRUE)  ## for independence from stats4
## require(numDeriv,quietly=TRUE) ## for hessian()

call.to.char <- function(x) {
    ## utility function
    x <- as.list(x)
    if (length(x)>1) x <- x[c(2,1,3)]
    paste(sapply(x,as.character),collapse="")
}

## FIXME: problem with bounds and formulae!
calc_mle2_function <- function(formula,
                               parameters,
                               links,
                               start,
                               parnames,
                               use.deriv=FALSE,
                               data=NULL,
                               trace=FALSE) {
  ## resid=FALSE
  ##  stub: what was I going to use this for ???
  ##  returning residuals rather than mle (e.g. for minpack.nls??)
  RHS <- formula[[3]]
  ddistn <- as.character(RHS[[1]])
  if (ddistn=="dnorm" && !("sd" %in% names(RHS))) {
      warning("using dnorm() with sd implicitly set to 1 is rarely sensible")
  }
  if (ddistn=="dnbinom" && !("mu" %in% names(RHS))) {
  }
  ## need to check on variable order:
  ## should it go according to function/formula,
  ##   not start?
  if (!is.list(data)) stop("must specify data argument",
                           " (as a list or data frame)",
                           " when using formula argument")
  vecstart <- (is.numeric(start))
  if (vecstart) start <- as.list(start) ## expand to a list
  if (missing(parnames) || is.null(parnames)) {
    parnames <- as.list(names(start))
    names(parnames) <- names(start)
  }
  ## hack
  if (!missing(parameters)) {
    ## linear model specified for some parameters
    vars <- as.character(sapply(parameters,"[[",2))
    if (length(parameters)>1) {
      models <-  sapply(parameters,function(z) call.to.char(z[[3]]))
    } else {
      models <- as.character(parameters)
    }
    models <- gsub(" ","",models)
    parameters <- parameters[models!="1"]
    npars <- length(parameters)
    if (npars==0) { ## no non-constant parameters
      parameters <- mmats <- vpos <- NULL
    } else {
      ## BUG IN HERE SOMEWHERE, FIXME: SENSITIVE TO ORDER OF 'start'
      mmats <- list()
      vpos <- list()
      pnames0 <- parnames
      names(parnames) <- parnames
      for (i in seq(along=parameters)) {
        vname <- vars[i]      ## name of variable
        p <- parameters[[i]]  ## formula for variable
        p[[2]] <- NULL
        mmat <- model.matrix(p,data=data)     
        pnames <- paste(vname,colnames(mmat),sep=".")
        parnames[[vname]] <- pnames ## insert into parameter names
        vpos0 <- which(pnames0==vname)
        vposvals <- cumsum(sapply(parnames,length))
        ## fill out start vectors with zeros or replicates as appropriate
        if (length(start[[vname]])==1) {
            if (length(grep("-1",models[i])>0)) {
                start[[vname]] <- rep(start[[vname]],length(pnames))
            } else {
                start[[vname]] <- c(start[[vname]],rep(0,length(pnames)-1))
            }
        }
        ## fix: what if parameters are already correctly specified?
        startpos <- if (vpos0==1) 1 else vposvals[vpos0-1]+1
        vpos[[vname]] <- startpos:vposvals[vpos0]
        mmats[[vname]] <- mmat
      }
    }
  } else parameters <- vars <- mmats <- vpos <- NULL
  if (!missing(links)) {
    stop("parameter link functions not yet implemented")
    for (i in length(links)) {
    }
  }
  parnames <- unlist(parnames)
  start <- as.list(unlist(start)) ## collapse/re-expand (WHY?)
  names(start) <- parnames
  arglist <- as.list(RHS[-1]) ## delete function name
  arglist$parameters <- NULL
  arglist1 <- c(list(x=formula[[2]]),arglist,list(log=TRUE))
  arglist1  ## codetools check kluge
  fn <- function() {
      ## is there a better way to do this?
      ## need to look for parameters etc.
      pars <- unlist(as.list(match.call())[-1])
      if (!is.null(parameters)) {
          for (.i in seq(along=parameters)) {
              assign(vars[.i],mmats[[.i]] %*% pars[vpos[[.i]]])
          }
      }
    ## if (is.null(data) || !is.list(data))
    ## stop("data argument must be specified when using formula interface")
    ## BUG/FIXME: data evaluates to 'FALSE' at this point -- regardless of whether
    ## it has been specified
    ## FIXME: how to make this eval() less fragile???
    ## sys.frame(sys.nframe()) specifies the number of the *current* frame
    ## ... envir=data,enclos=parent.frame()
      ## this actually works OK: fails enigmatically if we
      ## 
    arglist2 <- lapply(arglist1,eval,envir=data,
                       enclos=sys.frame(sys.nframe()))
    if (use.deriv) {
      stop("use.deriv is not yet implemented")
      ## browser()
      ## minor hack -- should store information otherwise -- could have
      ##  different numbers of arguments for different distributions?
      LLform <- get(gsub("^d","s",as.character(RHS[[1]])))(NA,NA)$formula
      avals <- as.list(formula[[3]][-1])
      for (i in seq_along(avals))
        LLform <- gsub(names(avals)[i],avals[[i]],LLform)
      r <- eval(deriv(parse(text=LLform),parnames),envir=c(arglist2,data))
    } else {
      r <- -sum(do.call(ddistn,arglist2))
    }
    ## doesn't work yet -- need to eval arglist in the right env ...
    ## if (debugfn) cat(unlist(arglist),r,"\n")
    if (trace) cat(pars,r,"\n")
    r
  }
  npars <- length(parnames)
  flist <-  vector("list",npars)
  names(flist) <- parnames
  ## add additional parnames?
  ## browser()
  ## flist <- c(flist,setdiff(names(arglist),c("x","log",... ?))
  formals(fn) <- flist
  if (vecstart) start <- unlist(start)
  list(fn=fn,start=start,parameters=parameters,
       fdata=list(vars=vars,mmats=mmats,vpos=vpos,
       arglist1=arglist1,ddistn=ddistn,parameters=parameters),
       parnames=parnames)
}

## need logic that will identify correctly when
## we need to pass parameters as a vector
mle2 <- function(minuslogl,
                 start,  ## =formals(minuslogl),
                 method,
                 optimizer,
                 fixed=NULL,
                 data=NULL,
                 subset=NULL,
                 default.start=TRUE, 
                 eval.only = FALSE,
                 vecpar = FALSE,
                 parameters=NULL,
                 parnames=NULL,
                 skip.hessian=FALSE,
                 hessian.opts=NULL,
                 use.ginv=TRUE,
                 trace=FALSE,
                 browse_obj=FALSE,
                 gr,
                 optimfun,
                 ...) {

  if (missing(method)) method <- mle2.options("optim.method")
  if (missing(optimizer)) optimizer <- mle2.options("optimizer")
  L <- list(...)
  if (optimizer=="optimize" && (is.null(L$lower) || is.null(L$upper)))
      stop("lower and upper bounds must be specified when using
'optimize'")
  if (inherits(minuslogl,"formula")) {
    pf <- function(f) {if (is.null(f))
                         {  ""
                          } else {
                            paste(f[2],"~",
                                  gsub(" ","",as.character(f[3])),sep="")
                          }
                     }
    if (missing(parameters)) {
      formula <- pf(minuslogl)
    } else {
      formula <- paste(pf(minuslogl),
                       paste(sapply(parameters,pf),collapse=", "),sep=": ")
    }
    tmp <- calc_mle2_function(minuslogl,parameters,
                              start=start,
                              parnames=parnames,
                              data=data,trace=trace)
    minuslogl <- tmp$fn
    start <- tmp$start
    fdata <- tmp$fdata
    parameters <- tmp$parameters
  } else {
    formula <- ""
    fdata <- NULL
  }
  call <- match.call()
  call.orig <- call
  ## ?? still not sure this is the best thing to do, but:
  ##   evaluate all elements of call
  ##    to make sure it will still function in new environments ...
  ## call[-1] <- lapply(call[-1],eval.parent)
  ## call[-1] <- lapply(call[-1],eval,envir=parent.frame(),enclos=parent.frame(2))
  ## FAILS if embedded in a funny environment (e.g. called from lapply)
  ##  why do we need this in the first place?
  ## FIXME: change update(), profile() to re-fit model properly
  ##  rather than evaluating call(), or generally find a less-fragile
  ##  way to do this.  Reverting to original form for now.
  call$data <- eval.parent(call$data)
  call$upper <- eval.parent(call$upper)
  call$lower <- eval.parent(call$lower)
  call$gr <- eval.parent(call$gr)
  ## FIX based on request from Mark Clements
  ## call$control$parscale <- eval.parent(call$control$parscale)
  ## call$control$ndeps <- eval.parent(call$control$ndeps)
  ## call$control$maxit <- eval.parent(call$control$maxit)
  call$control <- eval.parent(call$control)
  if(!missing(start))
    if (!is.list(start)) {
      if (is.null(names(start)) || !is.vector(start))
        stop("'start' must be a named vector or named list")
      ## do we want this or not???
      vecpar <- call$vecpar <- TRUE  ## given a vector start: set vecpar=TRUE
      start <- as.list(start)
    }
  ## also check parnames(minuslogl)?
  if (missing(start) && default.start) start <- formals(minuslogl)
  if (!is.null(fixed) && !is.list(fixed)) {
    if (is.null(names(fixed)) || !is.vector(fixed))
      stop("'fixed' must be a named vector or named list")
    fixed <- as.list(fixed)
  }
  if (!is.null(data) && !is.list(data)) ##  && !is.environment(data)) 
    stop("'data' must be a list")
  nfix <- names(unlist(namedrop(fixed)))
  if (!is.null(parnames(minuslogl))) {
    nfull <- parnames(minuslogl)
    fullcoef <- vector("list",length(nfull))
    names(fullcoef) <- nfull
  } else {
    fullcoef <- formals(minuslogl)
    nfull <- names(fullcoef)
  }
  if(any(! nfix %in% nfull))
    stop("some named arguments in 'fixed' are not arguments to the specified log-likelihood function")
  if (length(nfix)>0) start[nfix] <- NULL
  fullcoef[nfix] <- fixed
  ## switched namedrop() from outside to inside sapply ?
  nstart <- names(unlist(sapply(namedrop(start),eval.parent)))
  fullcoef[! nfull %in% nfix & ! nfull %in% nstart ] <- NULL  ## delete unnecessary names
  nfull <- names(fullcoef)
  lc <- length(call$lower)
  lu <- length(call$upper)
  npnfix <- sum(!nfull %in% nfix)
  if (!npnfix==0 && (lu>npnfix || lc>npnfix )) {
    warning("length mismatch between lower/upper ",
            "and number of non-fixed parameters: ",
            "# lower=",lc,", # upper=",lu,", # non-fixed=",npnfix)
  }
  template <- lapply(start, eval.parent)  ## preserve list structure!
  if (vecpar) template <- unlist(template)
  start <- sapply(namedrop(start), eval.parent) # expressions are allowed; added namedrop
  nstart <- names(unlist(namedrop(start)))
  ## named <- length(names(fullcoef))
  oo <- match(nstart, names(fullcoef))
  if (any(is.na(oo)))
    stop("some named arguments in 'start' are not arguments to the specified log-likelihood function")
  ## if (named)
  start <- start[order(oo)]
  ## rearrange lower/upper to same order as "start"
  ## FIXME: use names to rearrange if present
  fix_order <- function(c1,name,default=NULL) {
      if (!is.null(c1)) {
          if (length(unique(c1))>1) {  ## not all the same
            if (is.null(names(c1)) && length(unique(c1))>1) {
              warning(name," not named: rearranging to match 'start'")
              oo2 <- oo
            } else oo2 <- match(names(unlist(namedrop(c1))),names(fullcoef))
            c1 <- c1[order(oo2)]
          }
        } else c1 <- default
      c1
  }
  call$lower <- fix_order(call$lower,"lower bounds",-Inf)
  call$upper <- fix_order(call$upper,"upper bounds",Inf)
  call$control$parscale <- fix_order(call$control$parscale,"parscale")
  call$control$ndeps <- fix_order(call$control$ndeps,"ndeps")
  if (is.null(call$control)) call$control <- list()
  ## attach(data,warn.conflicts=FALSE)
  ## on.exit(detach(data))
  denv <- local(environment(),c(as.list(data),fdata,list(mleenvset=TRUE)))
  ## denv <- local(new.env(),c(as.list(data),fdata,list(mleenvset=TRUE)))
  argnames.in.data <- names(data)[names(data) %in% names(formals(minuslogl))]
  args.in.data <- lapply(argnames.in.data,get,env=denv)
  names(args.in.data) <- argnames.in.data
  args.in.data  ## codetools kluge
  objectivefunction <- function(p){
      if (browse_obj) browser()
      l <- relist2(p,template) ## redo list structure
      ## if (named)
      names(p) <- nstart[order(oo)] ## make sure to reorder
      ## ??? useless, comes after l is constructed ???
      l[nfix] <- fixed
      ##    cat("p\n"); print(p)
      ## cat("l\n"); print(l)
      ##    cat("data\n"); print(data)
      if (vecpar) {
          ## if (named)
          l <- namedrop(l[nfull])
          l <- unlist(l)
          args <- list(l)
          args <- c(list(l),args.in.data)
      } else { args <- c(l,args.in.data)
           }
      ## eval in environment of minuslogl???
      ## doesn't help, environment(minuslogl) is empty by this time
      ## cat("e3:",length(ls(envir=environment(minuslogl))),"\n")
      ## hack to remove unwanted names ...
      do.call("minuslogl",namedrop(args))
  } ## end of objective function
  objectivefunctiongr <-
    if (missing(gr)) NULL else
        function(p) {
          if (browse_obj) browser()
          l <- relist2(p,template) ## redo list structure
          names(p) <- nstart[order(oo)] ## make sure to reorder
          l[nfix] <- fixed
          if (vecpar) {
            l <- namedrop(l[nfull])
            l <- unlist(l)
            args <- list(l)
            args <- c(list(l),args.in.data)
          } else { args <- c(l,args.in.data)
                 }
          v <- do.call("gr",args)
          if (is.null(names(v))) {
              if (length(v)==length(l) && !is.null(tt <- names(l))) {
                  ## try to set names from template
                  vnames <- tt
              } else if (length(v)==length(p) && !is.null(tt <- names(p))) {
                  ## try to set names from params
                  vnames <- tt
              } else if (!is.null(tt <- parnames(minuslogl))) {
                  ## names were set as an attribute of the function
                  vnames <- tt
              } else vnames <- names(formals(minuslogl))
              if (length(vnames)!=length(v))
                  stop("name/length mismatch in gradient function")
              names(v) <- vnames
          }
          v[!names(v) %in% nfix] ## from Eric Weese
        } ## end of gradient function
  ## FIXME: try to do this by assignment into appropriate
  ##    environments rather than replacing them ...
  ## only set env if environment has not been previously set!
  if (!("mleenvset" %in% ls(envir=environment(minuslogl)))) {
      newenv <- new.env(hash=TRUE,parent=environment(minuslogl))
      d <- as.list(denv)
      mapply(assign,names(d),d,
             MoreArgs=list(envir=newenv))
      environment(minuslogl) <- newenv
      if (!missing(gr)) {
          newenvgr <- new.env(hash=TRUE,parent=environment(minuslogl))
          mapply(assign,names(d),d,
                 MoreArgs=list(envir=newenvgr))
      }
  }
  if (length(start)==0 || eval.only) {
    if (length(start)==0) start <- numeric(0)
    optimizer <- "none"
    skip.hessian <- TRUE
    oout <- list(par=start, value=objectivefunction(start),
                 hessian = matrix(NA,nrow=length(start),ncol=length(start)))
  } else {
    oout <- switch(optimizer,
                   optim = {
                       arglist <- list(...)
                       arglist$lower <- arglist$upper <-
                         arglist$control <- NULL
                       do.call("optim",
                               c(list(par=start,
                                      fn=objectivefunction,
                                      method=method,
                                      hessian=FALSE,
                                      gr=objectivefunctiongr,
                                      control=call$control,
                                      lower=call$lower,
                                      upper=call$upper),
                                 arglist))
                   },
                   optimx = {
                     ## don't ask, will get us into
                     ##   dependency hell
                     ## require("optimx")
                     arglist <- list(...)
                     arglist$lower <- arglist$upper <-
                       arglist$control <- NULL
                     do.call("optimx",
                             c(list(par=start,
                                    fn=objectivefunction,
                                    method=method,
                                    hessian=FALSE,
                                    gr=objectivefunctiongr,
                                    control=call$control,
                                    lower=call$lower,
                                    upper=call$upper),
                               arglist))
                   },
                   nlm = nlm(f=objectivefunction, p=start, hessian=FALSE, ...),
                     ##!skip.hessian,
                   ## <http://netlib.bell-labs.com/cm/cs/cstr/153.pdf>
                   nlminb = nlminb(start=start,
                     objective=objectivefunction, hessian=NULL, ...),
                   constrOptim = constrOptim(theta=start,
                     f=objectivefunction, method=method, ...),
                   optimize=,
                   optimise= optimize(f=objectivefunction,
                     interval=c(call$lower,call$upper), ...),
                   user = {
                     arglist <- list(...)
                     arglist$lower <- arglist$upper <-
                       arglist$control <- NULL
                     do.call(optimfun,
                             c(list(par=start,
                                    fn=objectivefunction,
                                    method=method,
                                    hessian=FALSE,
                                    gr=objectivefunctiongr,
                                    control=call$control,
                                    lower=call$lower,
                                    upper=call$upper),
                               arglist))
                   },
                   stop("unknown optimizer (choices are 'optim', 'nlm', 'nlminb', 'constrOptim', 'user', and 'optimi[sz]e')")
                 )
  }
  optimval <- switch(optimizer,
                     optim= , constrOptim=, optimx=, user=, none="value",
                     nlm="minimum",
                     optimize=, optimise=, nlminb="objective")
  if (optimizer=="optimx") {
    fvals <- oout[["value"]]
    conv <- oout[["convcode"]]
    ## best <- if (!any(conv==0)) {
    best <- which.min(fvals)
    ##    } else {
    ## fvals <- fvals[conv==0]
    ## which.min(fvals)
    ## }
    oout <- list(par=as.numeric(unlist(oout[best,1:attr(oout,"npar")])),
                 value=fvals[best],
                 convergence=conv[best],
                 method.used=attr(oout,"details")[,"method"][[best]])
    ## FIXME: should do profiles only with best method for MLE?
  }
  if (optimizer=="nlm") {
    oout$par <- oout$estimate
    oout$convergence <- oout$code
  }
  if (optimizer %in% c("optimise","optimize")) {
    oout$par <- oout$minimum
    oout$convergence <- 0 ## can't detect non-convergence
  }
  if (optimizer %in% c("nlminb","optimise","optimize") ||
      ## optimizer (bobyqa?) may have stripped names -- try to restore them!
      is.null(names(oout$par))) {
    names(oout$par) <- names(start)
  }

  ## FIXME: worry about boundary violations?
  ## (if we're on the boundary then the Hessian may not be useful anyway)
  ##
  if (length(oout$par)==0) skip.hessian <- TRUE
  if (!skip.hessian) {
    if ((!is.null(call$upper) || !is.null(call$lower)) &&
        any(oout$par==call$upper) || any(oout$par==call$lower))
      warning("some parameters are on the boundary: variance-covariance calculations based on Hessian may be unreliable")
  }
  namatrix <- matrix(NA,nrow=length(start),ncol=length(start))
  if (!skip.hessian) {
    psc <- call$control$parscale
    if (is.null(psc)) {
      oout$hessian <- try(hessian(objectivefunction,oout$par,method.args=hessian.opts))
    } else {
      tmpf <- function(x) {
        objectivefunction(x*psc)
      }
      oout$hessian <- try(hessian(tmpf,oout$par/psc,method.args=hessian.opts))/outer(psc,psc)
    }
  }
  if (skip.hessian || inherits(oout$hessian,"try-error"))
      oout$hessian <- namatrix
  coef <- oout$par
  nc <- names(coef)
  if (skip.hessian) {
    tvcov <- matrix(NA,length(coef),length(coef))
  } else {
    if (length(coef)) {
      if (use.ginv) {
        tmphess <- try(MASS::ginv(oout$hessian),silent=TRUE)
      } else {
        tmphess <- try(solve(oout$hessian,silent=TRUE))
      }
      if (class(tmphess)=="try-error") {
        tvcov <- matrix(NA,length(coef),length(coef))
        warning("couldn't invert Hessian")
      } else tvcov <- tmphess
    } else {
      tvcov <- matrix(numeric(0),0,0)
    }
  }
  dimnames(tvcov) <- list(nc,nc)
  min <-  oout[[optimval]]
  ##  if (named)
  fullcoef[nstart[order(oo)]] <- coef
  ## else fullcoef <- coef
  ## compute termination info
  ## FIXME: should we worry about parscale here??
  if (length(coef)) {
    gradvec <- if (!missing(gr)) {
        objectivefunctiongr(coef)
    } else {
        if (inherits(tt <- try(grad(objectivefunction,coef),silent=TRUE),
                     "try-error")) NA else tt
    }
    oout$maxgrad <-  max(abs(gradvec))
    if (!skip.hessian) {
        if (inherits(ev <- try(eigen(oout$hessian)$value,silent=TRUE),
                     "try-error")) ev <- NA
        oout$eratio <- min(ev)/max(ev)
    }
  }
  if (!is.null(conv <- oout$conv) &&
      ((optimizer=="nlm" && conv>2) ||
       (optimizer!="nlm" && conv!=0))) {
      ## warn of convergence failure
      if (is.null(oout$message)) {
          cmsg <- "unknown convergence failure: refer to optimizer documentation"
          if (optimizer=="optim") {
              if (conv==1) cmsg <- "iteration limit 'maxit' reached"
              if (conv==10) cmsg <- "degenerate Nelder-Mead simplex"
          } else if (optimizer=="nlm") {
              if (conv==3) cmsg <- "last global step failed to locate a point lower than 'estimate': see ?nlm"
              if (conv==4) cmsg <- "iteration limit exceeded"
              if (conv==5) cmsg <- "maximum step size 'stepmax' exceeded five consecutive times: see ?nlm"
          }
      } else cmsg <- oout$message
      warning(paste0("convergence failure: code=",conv," (",cmsg,")"))
  }
  m <- new("mle2", call=call, call.orig=call.orig, coef=coef, fullcoef=unlist(fullcoef), vcov=tvcov,
           min=min, details=oout, minuslogl=minuslogl, method=method,
           optimizer=optimizer,data=as.list(data),formula=formula)
  attr(m,"df") = length(m@coef)
  if (!missing(data)) attr(m,"nobs") = length(data[[1]])
  ## to work with BIC as well
  m
}


get.mnames <- function(Call) {
    xargs <- which(names(Call) %in% names(formals(ICtab))[-1])
    mnames <- as.character(Call)[c(-1,-xargs)]
    if (length(mnames)==1) {
        g <- get(mnames)
        if (is.list(g) && length(g)>1) {
            if (is.null(names(g))) mnames <- paste("model",1:length(g),sep="")
            else mnames <- names(g)
            if (any(duplicated(mnames))) stop("model names must be distinct")
        }
    }
    mnames
}
  

mle2.options <- function(...) {
single <- FALSE
args <- list(...)
  setvals <- !is.null(names(args))
  if (!length(args)) args <- names(.Mle2.options)
if (all(unlist(lapply(args, is.character)))) 
     args <- as.list(unlist(args))
  if (length(args) == 1) {
    if (is.list(args[[1]]) | is.null(args[[1]])) 
      args <- args[[1]]
    else if (!setvals)
      single <- TRUE
  }
  if (setvals) {
    .Mle2.options[names(args)] <<- args
    value <- .Mle2.options[names(args)]
  } else value <- .Mle2.options[unlist(args)]
   if (single) value <- value[[1]]
if (setvals) invisible(value) else value
}


.Mle2.options = list(optim.method="BFGS",confint = "spline",optimizer="optim")


## .onLoad <- function(lib, pkg) require(methods)


## (not yet) replaced by relist?
## reconstruct list structure:
##   v is a vector, l is the original list
##      to use as a template
relist2 <- function(v,l) {
 if (is.list(v)) v <- unlist(v)
 if (!all(sapply(l,mode)=="numeric")) {
     stop("can't relist non-numeric values")
   }
   lens = sapply(l,length)
   if (all(lens==1))
     return(as.list(v))
   l2 <- split(v,rep(1:length(l),lens))
   names(l2) <- names(l)
   l3 <- mapply(function(x,y) {
     if (!is.null(dim(y))) {
       z=array(x,dim(y)); dimnames(z)=dimnames(y); z
     } else {
       z=x; names(z)=names(y); z
     }
   },l2,l,SIMPLIFY=FALSE)
 names(l3) <- names(l)
   l3
 }

namedrop <- function(x) {
if (!is.list(x)) x
  for (i in seq(along=x)) {
    ## cat(i,length(x),"\n")
    n = names(x[[i]])
    lx = length(x[[i]])
    if (!is.null(n)) {
      if (lx==1) {
        names(x[[i]]) <- NULL
      } else if (length(unique(n))<lx) {
        names(x[[i]]) <- 1:lx
      }
    } ## !is.null(names(x[[i]]))
  } ## loop over elements
  x
}

"parnames<-" <- function(obj,value) {
attr(obj,"parnames") <- value
  obj
}

parnames <- function(obj) {
attr(obj,"parnames")
}


### TEST OF NLMINB
if (FALSE) {
x <- 0:10
y <- c(26, 17, 13, 12, 20, 5, 9, 8, 5, 4, 8)
mle2(y~dpois(lambda=ymax/(1+x/xhalf)),start=list(ymax=25,xhalf=3.06),
     optimizer="nlminb",fixed=list(ymax=38.76),lower=c(0,0),trace=TRUE)

f = calc_mle2_function(y~dpois(lambda=ymax/(1+x/xhalf)),
  start=list(ymax=25,xhalf=3.06))
f2 = function(xhalf) {
  -sum(dpois(y,38.76/(1+x/xhalf),log=TRUE))
}
optim(f2,par=3.06,method="BFGS")
## optim(f2,par=3.06,method="L-BFGS-B",lower=0) ## error
nlminb(objective=f2,start=3.06)
nlminb(objective=f2,start=3.06,lower=0)
nlminb(objective=f2,start=3.06,lower=1e-8)
}

strwrapx <- function(x, width = 0.9 * getOption("width"),
indent = 0, exdent = 0, 
prefix = "", simplify = TRUE,
    parsplit= "\n[ \t\n]*\n", wordsplit = "[ \t\n]") 
{
    if (!is.character(x)) 
      x <- as.character(x)
    indentString <- paste(rep.int(" ", indent), collapse = "")
    exdentString <- paste(rep.int(" ", exdent), collapse = "")
    y <- list()
    ## split at "+" locations
    plussplit = function(w) {
      lapply(w,
             function(z) {
               plusloc = which(strsplit(z,"")[[1]]=="+")
               plussplit =
                 apply(cbind(c(1,plusloc+1),
                             c(plusloc,nchar(z,type="width"))),
                       1,
                       function(b) substr(z,b[1],b[2]))
               plussplit
             })}
    ## ugh!
    z <- lapply(strsplit(x, parsplit),
                function(z) { lapply(strsplit(z,wordsplit),
                                   function(x) unlist(plussplit(x)))
                })
    ## print(str(lapply(strsplit(x,parsplit),strsplit,wordsplit)))
    ## print(str(z))
    for (i in seq_along(z)) {
        yi <- character(0)
        for (j in seq_along(z[[i]])) {
            words <- z[[i]][[j]]
            nc <- nchar(words, type = "w")
            if (any(is.na(nc))) {
                nc0 <- nchar(words)
                nc[is.na(nc)] <- nc0[is.na(nc)]
            }
            if (any(nc == 0)) {
                zLenInd <- which(nc == 0)
                zLenInd <- zLenInd[!(zLenInd %in% (grep("\\.$", 
                  words) + 1))]
                if (length(zLenInd) > 0) {
                  words <- words[-zLenInd]
                  nc <- nc[-zLenInd]
                }
            }
            if (length(words) == 0) {
                yi <- c(yi, "", prefix)
                next
            }
            currentIndex <- 0
            lowerBlockIndex <- 1
            upperBlockIndex <- integer(0)
            lens <- cumsum(nc + 1)
            first <- TRUE
            maxLength <- width - nchar(prefix, type = "w") - 
                indent
            while (length(lens) > 0) {
                k <- max(sum(lens <= maxLength), 1)
                if (first) {
                  first <- FALSE
                  maxLength <- maxLength + indent - exdent
                }
                currentIndex <- currentIndex + k
                if (nc[currentIndex] == 0) 
                  upperBlockIndex <- c(upperBlockIndex, currentIndex - 
                    1)
                else upperBlockIndex <- c(upperBlockIndex, currentIndex)
                if (length(lens) > k) {
                  if (nc[currentIndex + 1] == 0) {
                    currentIndex <- currentIndex + 1
                    k <- k + 1
                  }
                  lowerBlockIndex <- c(lowerBlockIndex, currentIndex + 
                    1)
                }
                if (length(lens) > k) 
                  lens <- lens[-(1:k)] - lens[k]
                else lens <- NULL
            }
            nBlocks <- length(upperBlockIndex)
            s <- paste(prefix, c(indentString, rep.int(exdentString, 
                nBlocks - 1)), sep = "")
            for (k in (1:nBlocks)) {
              s[k] <- paste(s[k],
                            paste(words[lowerBlockIndex[k]:upperBlockIndex[k]], 
                                  collapse = " "), sep = "")
            }
            s = gsub("\\+ ","+",s)  ## kluge
            yi <- c(yi, s, prefix)
        }
        y <- if (length(yi)) 
            c(y, list(yi[-length(yi)]))
        else c(y, "")
    }
    if (simplify) 
        y <- unlist(y)
    y
  }



