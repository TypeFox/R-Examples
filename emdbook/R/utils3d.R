is_installed <- function(pkg) {
    pkg %in% installed.packages()[,"Package"]
}

apply2d <-   function(fun,x,y,...,use_plyr=TRUE,.progress="none") {
    if (is.character(fun)) fun <- get(fun)
    if (!use_plyr) {
        a <- apply(expand.grid(x,y),1,function(z) { fun(z[1],z[2],...)})
    } else {
        a <- aaply(as.matrix(expand.grid(x,y)),1,
                   function(z) {fun(z[1],z[2],...)},
                   .progress=.progress)
    }
    return(matrix(a,nrow=length(x)))
}

## FIXME: think about the meaning of ... more carefully
curve3d <- function (expr, from=c(0,0), to=c(1,1),
                     n = c(41,41),
                     xlim, ylim, add = FALSE, 
                     xlab=varnames[1],ylab=varnames[2],
                     zlab = NULL,
                     log = NULL, 
                     sys3d = c("persp","wireframe","rgl","contour","image",
                       "none"),
                     varnames = c("x","y"),
                     use_plyr=TRUE,
                     .progress="none",
                     ...) 
{
    vars <- lapply(as.list(varnames),parse,file="",n=NULL)
    sys3d <- match.arg(sys3d)
    if (add && !(sys3d %in% c("contour","rgl")))
      stop("can only add contour or rgl to 3D plots")
    n = rep(n,length.out=2)
    if (!missing(xlim) && !missing(ylim)) {
      if (!missing(from) | !missing(to)) {
        stop("must specify one of (xlim,ylim) or (from,to)")
      }
      from <- c(xlim[1],ylim[1])
      to <- c(xlim[2],ylim[2])
    }
    sexpr <- substitute(expr)
    if (is.name(sexpr)) {
        fcall <- paste(sexpr, "(",varnames[1],",",varnames[2],")",sep="")
        expr <- parse(text = fcall)
        if (is.null(zlab))
          zlab <- fcall
  } else {
    if (!(is.call(sexpr) && 
          (match(varnames[1], all.vars(sexpr), nomatch = 0) || 
           match(varnames[2], all.vars(sexpr), nomatch = 0))))
      stop(paste("'expr' must be a function or an expression containing '",
           varnames[1],"' and '",varnames[2],"'",sep=""))
    expr <- sexpr
    if (is.null(zlab)) 
      zlab <- deparse(sexpr)
  } 
  lg <- if (length(log)) 
    log
  if (length(lg) == 0) 
    lg <- ""
    x <- if (lg != "" && "x" %in% strsplit(lg, NULL)[[1]]) {
        if (any(c(from[1], to[1]) <= 0)) 
            stop("'from[1]' and 'to[1]' must be > 0 with log=\"x\"")
        exp(seq(log(from[1]), log(to[1]), length = n))
      }  else seq(from[1], to[1], length = n[1])
    y <- if (lg != "" && "y" %in% strsplit(lg, NULL)[[1]]) {
        if (any(c(from[2], to[2]) <= 0)) 
            stop("'from[2]' and 'to[2]' must be > 0 with log=\"y\"")
        exp(seq(log(from[2]), log(to[2]), length = n))
    } else seq(from[2], to[2], length = n[2])
    tmpfun <- function(x,y) {
        env <- list(x,y)
        names(env) <- varnames
        ## SKIP the inside of curve3d, go back one more level ...
        eval(expr, envir = env, enclos = parent.frame(2))
    }
    z <- apply2d(tmpfun,x,y,use_plyr=use_plyr,.progress=.progress)
    switch(sys3d,
           persp=persp(x,y,z,xlab=xlab,ylab=ylab,zlab=zlab,...),
           contour=contour(x,y,z,xlab=xlab,ylab=ylab,add=add,...),
           image=image(x,y,z,xlab=xlab,ylab=ylab,...),
           none=NA,
           wireframe={
               print(wireframe(z,row.values=x,col.values=y,
                               xlab=xlab,ylab=ylab,zlab=zlab,...))},
           rgl={
              if (!is_installed("rgl")) {
                  stop("please install the rgl package first ...")
              } 
              rgl::persp3d(x,y,z,
                        xlab=xlab,
                        ylab=ylab,zlab=zlab,add=add,...)})
  invisible(list(x=x,y=y,z=z))
}

gridsearch2d <- function(fun,v1min,v2min,v1max,v2max,n1=20,n2=20,
                         logz=FALSE,sys3d=c("both","contour","image"),...) {
  sys3d=match.arg(sys3d)
  redraw = function(v1vec,v2vec,m) {
    if (sys3d=="contour") {
      contour(v1vec,v2vec,m,xlab="",ylab="",...)
    } else image(v1vec,v2vec,m,xlab="",ylab="",...)
    if (sys3d=="both") contour(v1vec,v2vec,m,add=TRUE,...)
  }
  recalc = function(v1vec,v2vec,logz) {
    m = apply2d(fun,v1vec,v2vec)
    mindm = pmax(diff(sort(m)),1e-10)
    if (logz) m = log10(m-min(m)+mindm)
    m
  }
  v1vec = seq(v1min,v1max,length=n1)
  v2vec = seq(v2min,v2max,length=n2)
  m = recalc(v1vec,v2vec,logz)
  stop1 = FALSE
  first = TRUE
  while (!stop1) {
    redraw(v1vec,v2vec,m)
    if (!first) {
      resp = readline("Continue (y|n)? ")
      stop1 = (toupper(resp)=="N")
    }
    first = FALSE
    if (!stop1) {
      cat("click on box corners\n")
      z = lapply(locator(2),sort)
      rect(z$x[1],z$y[1],z$x[2],z$y[2])
      resp = readline("OK (y|n)? ")
      if (toupper(resp)=="Y") {
        v1min=z$x[1]
        v1max=z$x[2]
        v2min=z$y[1]
        v2max=z$y[2]
        v1vec = seq(v1min,v1max,length=n1)
        v2vec = seq(v2min,v2max,length=n2)
        m = recalc(v1vec,v2vec,logz)
      }
    }
  }
  resp = readline("get point (y|n)? ")
  if (toupper(resp)=="Y") {
    cat("click on point\n")
    z=locator(1)
    return(z)
  } else invisible(NULL)
}

metropSB <- function(fn,start,deltap=NULL,
                      scale=1,rptfreq=-1,
                      acceptscale=1.01,rejectscale=0.99,
                      nmax=10000,retvals=FALSE,retfreq=100,
                      verbose=FALSE,...) {
  ## initialization
  ndim <- length(start)  # number of parameters/dimensions
  p <- start             # current parameter vector = starting vector
  minp <- start          # parameter vector corresponding to minimum value so far
  val <- fn(p,...)     # current function value
  minval <- val          # minimum value so far
  if (retvals) {
    info <- matrix(nrow=round(nmax/retfreq),ncol=3*ndim+3)
    dimnames(info) <- list(NULL,c(paste("p",1:ndim,sep=""),
                                  paste("minp",1:ndim,sep=""),
                                  paste("deltap",1:ndim,sep=""),
                                        "val","minval","accept"))
    # save info on parameter values, function value, whether accepted or not, current best parameters,
    # current best value, current jump size
  }
  it <- 1                # iteration count
  if (is.null(deltap)) deltap <- p*0.05   # default starting jump size(s)
  while (it<=nmax) {
    oldp <- p                             # save current parameters and function value
    oldval <- val
    p <- p + runif(ndim,-1,1)*deltap      # perturb current values
    val <- fn(p,...)                    # new function value
    dval <- val-oldval                    # change
    saveinfo <- (retvals && it %% retfreq == 0)
    savecount <- it %/% retfreq
    if (saveinfo) {
      info[savecount,-ncol(info)] <- c(p,minp,deltap,val,minval)
    }
    if (verbose)
      if (it %% 100 == 0)
         cat(it,"\n")
    if (dval<0 || (exp(-dval*scale)>runif(1))) {
      ## accept new value
      if (saveinfo) info[savecount,ncol(info)] <- 1
      if (val<minval) {  # if better than best so far, record new minimum
	minval <- val
	minp <- p
      }
      deltap <- deltap*acceptscale   # increase jump size
    } else {
      ## replace old values
      if (saveinfo) info[savecount,ncol(info)] <- 0
      val <- oldval
      p <- oldp
      deltap <- deltap*rejectscale   # decrease jump size
    }
    it <- it+1
    if ((rptfreq>0) && (it %% rptfreq)==1)
      cat("it=",it," value=",val," min. val.=",minval,"\n")
  }
  if (retvals)
    return(list(minimum=minval,estimate=minp,funcalls=it,retvals=info))
  else
    return(list(minimum=minval,estimate=minp,funcalls=it))
}

contour3d <- function(x,y,z,contourArgs=NULL,...) {
    if (!is_installed("rgl")) {
        stop("must install the rgl package to use contour3d")
    }
  if (is.list(x)) {
    if (!all(sort(names(x))==c("x","y","z")))
      stop("list should contain components 'x', 'y', 'z'")
    z <- x$z
    y <- x$y
    x <- x$x
  }
  ccc = do.call(contourLines,c(list(x,y,z),contourArgs))
  ccc = lapply(ccc, function(x) {
    list(x=x$x,y=x$y,z=rep(x$level,length(x$x))) })
  invisible(mapply(rgl::lines3d,ccc,MoreArgs=list(...)))
  invisible(ccc)
}
