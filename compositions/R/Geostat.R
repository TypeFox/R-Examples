


logratioVariogram <- function(comp,
                          loc,
                          maxdist=max(dist(loc))/2,
                          nbins=20,
                          dists=seq(0,maxdist,length.out=nbins+1),
                          bins=cbind(dists[-length(dists)],dists[-1]),
                          azimuth=0,
                          azimuth.tol=180  
                          ) {
  lcomp <- log(ifelse(is.NMV(comp),comp,NaN))
  if( !is.matrix(loc))
    loc <- as.matrix(loc)
  stopifnot(nrow(loc)==nrow(comp))
  stopifnot(ncol(bins)==2)
  nbins = nrow(bins)
  D     = ncol(lcomp)
  labs   = colnames(comp)
  erg <- .C("gsiCGSvariogram",
            dimZ=gsiInt(dim(lcomp),2),
            Z   =lcomp,
            dimX=gsiInt(dim(loc),2),
            X=gsiDouble(loc,prod(dim(loc))),
            nbinsK=gsiInt(dim(bins),2),
            bins=bins,
            aziVec=gsiDouble(if( length(azimuth)== 1 ) c(sin(azimuth/180*pi),cos(azimuth/180*pi)) else azimuth,ncol(loc)),
            cosMin=gsiDouble(cos(azimuth.tol/180*pi),1),
            vg=numeric(nbins*D^2),
            h =numeric(nbins*D^2),
            n =integer(nbins*D^2),
            NAOK=TRUE
            )
  with(erg,
       structure(
                 list(
                      vg=structure(vg,
                        dim=c(nbins,D,D),
                        dimnames=list(NULL,labs,labs)
                        ),
                      h =structure(h,
                        dim=c(nbins,D,D),
                        dimnames=list(NULL,labs,labs)
                        ),
                      n =structure(h,
                        dim=c(nbins,D,D),
                        dimnames=list(NULL,labs,labs)
                        )
                      ),
                 class="logratioVariogram"
                 )
       )

}

cgram2vgram <- function(cgram) {
  function(h,...) {
    cgram(0*h,...)-cgram(h,...)
  }
}

vgram2lrvgram <- function(vgram) {
  function(h,...) {
    vg <- vgram(h,...)
    lr <- .C("gsiCGSvg2lrvg",
       dimVg=gsiInt(dim(vg),3),
       vg   =gsiDouble(vg),
       lr   =numeric(prod(dim(vg)))
       )$lr
    attributes(lr)<-attributes(vg)
    lr
  }
}


plot.logratioVariogram <- function(x,...,type="l",lrvg=NULL,fcols=2:length(lrvg),oma=c(4, 4, 4, 4),gap=0,ylim=NULL) {

  D <- dim(x$vg)[2]
  lab <- dimnames(x$vg)[[2]]
  hmax <- max(x$h)
  hh <- seq(0,hmax,length.out=101)[-1]
  opar <- par(mfrow = c(D,D), 
            mar = rep.int(gap/2, 4), oma = oma)
  on.exit({
    par(opar)
  }, add = TRUE)
  if( !is.list(lrvg) && !is.null(lrvg) )
    lrvg <- list(lrvg)
  lrvgH <- lapply(lrvg,function(lrvg) {
    lrvg(hh)
    })
  if( is.null(ylim) )
    ylim <- c(0,max(max(x$vg),unlist(sapply(lrvgH,max),recursive = FALSE)))
  for(i in 1:D)
    for(j in 1:D) {
      if( i != j ) {
        h <- x$h[,i,j]
        vg<- x$vg[,i,j]
        plot(h,vg,type="l",...,ylim=ylim,
             yaxt=if(j==1) "s" else "n",
             xaxt=if(i==D) "s" else "n")
        if( length(lrvgH) > 0 )
          for(k in 1:length(lrvgH)) 
            lines(hh,lrvgH[[k]][,i,j],col=fcols[k])
      } else {
        plot.new()
        gsi.textpanel(lab=lab[i])
      }
    }
}

vgmGetParameters <- function(vg,envir=environment(vg)) {
  form <- lapply(formals(vg)[-1],eval,envir=envir)
  fi   <- if(length(form)>0) 1:length(form) else c()
  nams <- lapply(fi,function(i) {
    foi <- eval(form[[i]])
    myName <- names(form)[i]
    nam <- names(foi)
    if( length(nam) != length(foi))
      nam <- if(length(foi)>0) 1:length(foi) else c()
    if( length(nam)== 1)
      myName
    else
      paste(myName,nam,sep=".")
  }
                 )
  structure(unlist(form,recursive = FALSE),names=unlist(nams,recursive = FALSE))
}

vgmSetParameters <- function(vg,p) {
  form <- formals(vg)[-1]
  fi   <- if(length(form)>0) 1:length(form) else c()
  nn   <- cumsum(c(1,sapply(form,function(x) length(eval(x)))))
  for(i in fi) {
    formals(vg)[[i+1]]<-p[nn[i]:(nn[i+1]-1)]
  }
  vg
}

vgmGof<- function (p = vgmGetParameters(vg), emp, vg, mode = "log") 
{
    h <- apply(emp$h, c(1), mean)
    vg2 <- vgmSetParameters(vg, p)
    vgv <- vgram2lrvgram(vg2)(h)
    D <- dim(emp$vg)[2]
    N <- dim(emp$vg)[1]
    if (mode == "log") {
      tk <- rep(diag(D),each=N)<0.5
      ratio <- (emp$vg[tk]/vgv[tk])
      sum(log(ratio)^2 * sqrt(emp$n[tk]))
    } else sum((emp$vg - vgv)^2 * sqrt(emp$n))
}

vgmFit2lrv <- function(emp,vg,...,mode="log",psgn=rep(-1,length(param)),print.level=1) {
  h <- apply(emp$h,c(1),mean)
  param <- vgmGetParameters(vg)
  pa <- ifelse(psgn>0,log(param),ifelse(psgn<0,param,sqrt(param)))
  D <- dim(emp$vg)[2]
  N <- dim(emp$vg)[1]
  gof <- function(p) {
    vg2 <- vgmSetParameters(vg,ifelse(psgn>0,exp(p),ifelse(psgn<0,p,p^2)))
    vgv <- vgram2lrvgram(vg2)(h)
    if( mode=="log" ) {
      tk <- rep(diag(D),each=N)<0.5
      ratio <- (emp$vg[tk]/vgv[tk])
      sum(log(ratio)^2 * sqrt(emp$n[tk]))
    } else
      sum((emp$vg-vgv)^2*sqrt(emp$n))
  }
  erg <- nlm(gof,param,...,print.level=print.level)
  p <- erg$estimate
  list(nlm=erg,vg=vgmSetParameters(vg,ifelse(psgn>0,exp(p),ifelse(psgn<0,p,p^2))))
}
vgmFit<-function(emp,vg,...,mode="log",psgn=rep(-1,length(param)),print.level=1) {
  .Deprecated("vgmFit2lrv")
  h <- apply(emp$h,c(1),mean)
  param <- vgmGetParameters(vg)
  pa <- ifelse(psgn>0,log(param),ifelse(psgn<0,param,sqrt(param)))
  D <- dim(emp$vg)[2]
  N <- dim(emp$vg)[1]
  gof <- function(p) {
    vg2 <- vgmSetParameters(vg,ifelse(psgn>0,exp(p),ifelse(psgn<0,p,p^2)))
    vgv <- vgram2lrvgram(vg2)(h)
    if( mode=="log" ) {
      tk <- rep(diag(D),each=N)<0.5
      ratio <- (emp$vg[tk]/vgv[tk])
      sum(log(ratio)^2 * sqrt(emp$n[tk]))
    } else
      sum((emp$vg-vgv)^2*sqrt(emp$n))
  }
  erg <- nlm(gof,param,...,print.level=print.level)
  p <- erg$estimate
  list(nlm=erg,vg=vgmSetParameters(vg,ifelse(psgn>0,exp(p),ifelse(psgn<0,p,p^2))))
}

gsih2Dist <- function(h,A=NULL) {
  if( is.data.frame(h) )
    h <- as.matrix(h)
  if( is.matrix(h) )
    h <- norm(rmult(h))
  h
}

# The isotropic spherical variogram
vgram.sph <- function( h , nugget = 0, sill = 1, range= 1,... ) {
  "Sphaerisches Variogramm"
  h <- gsih2Dist(h)
  s <- sill-nugget
  ifelse(h>range*1E-8,
         nugget+s*ifelse(h>range,1,1.5*h/range-0.5*(h/range)^3)
         ,0)
}

# The isotropic exponential variogram
vgram.exp <- function( h , nugget = 0, sill = 1, range= 1,... ) {
  "Exponentielles Variogramm"
  s <- sill-nugget
  r <- -range/log(0.1)
  h <- gsih2Dist(h)
  ifelse(h>range*1E-8,nugget+s*(1-exp(-h/r)),0)
}

# The isotropic gaussian variogram
vgram.gauss <- function( h , nugget = 0, sill = 1, range= 1,... ) {
  "Gausssches Variogramm"
  h <- gsih2Dist(h)
  s <- sill-nugget
  r <- range/sqrt(-log(0.1))
  ifelse(h>range*1E-8,nugget+s*(1-exp(-(h/r)^2)),0)
}

# The isotropic linear variogram
vgram.lin <- function( h , nugget = 0, sill = 1, range= 1,...  ) {
  "Lineares Variogramm"
  h <- gsih2Dist(h)
  s <- sill-nugget  
  ifelse(h>range*1E-8,nugget+h*s/range,0)
} # sill und range geben einen Punkt an der durchlaufen wird

# The isotropic power variogram
vgram.pow <- function( h , nugget = 0, sill = 1, range= 1,... ) {
  "Power Variogramm"
  h <- gsih2Dist(h)
  ifelse(h>range*1E-8,nugget+sill*h^range,0)
} # sill und range unsinnige parameter

vgram.nugget <- function( h, nugget = 1 ,...,tol=1E-8) {
  h <- gsih2Dist(h)
  ifelse(h>tol,nugget,0)
}


  
parametricPosdefMat <- function(p) {
  D = 0.5*(-1+sqrt(1+8*length(p))) # p=(D^2-D)/2+D 0=D^2+D-2p
  mat <- matrix(0,nrow=D,ncol=D)
  mat[row(mat)<=col(mat)]<-p
  t(mat) %*% mat
}

parametricPosdefClrMat <- function(p) {
  D = 0.5*(-1+sqrt(1+8*length(p))) # p=(D^2-D)/2+D 0=D^2+D-2p
  mat <- matrix(0,nrow=D,ncol=D)
  mat[row(mat)<=col(mat)]<-p
  ilrvar2clr(t(mat) %*% mat)
}


parametricRank1Mat <- function(p) {
  p %o% p
}

parametricRank1ClrMat <- function(p) {
  ilrvar2clr(p %o% p)
}

parameterRank1Mat <- function(A) {
  with(eigen(A,symmetric=TRUE),vectors[,1]*sqrt(max(values[1],0)))
}

parameterRank1ClrMat <- function(A) {
  with(eigen(clrvar2ilr(A),symmetric=TRUE),vectors[,1]*sqrt(max(values[1],0)))
}


parameterPosdefMat <- function(A) {
  R <- chol(A,pivot=FALSE)
  #pivot <- attr(R,"pivot")
  #oo    <- order(pivot)
  #R<-R[,oo]
  R[row(R)<=col(R)]
}

parameterPosdefClrMat <- function(A) {
  R <- chol(clrvar2ilr(A),pivot=FALSE)
  #pivot <- attr(R,"pivot")
  #oo    <- order(pivot)
  #R<-R[,oo]
  R[row(R)<=col(R)]
}


gsiLMCRassum <- function(form) {
  lf = length(form)
  if( lf == 1 )
    return(list(form))
  if( lf == 2 )
    return(list(form))
  if( lf == 3 )
    if(is.call(form) ) 
      if(as.character(form[[1]])=="+" )  
        return( c(Recall(form[[2]]),Recall(form[[3]]) ) )
      else if(as.character(form[[1]])=="-")
        return( c(Recall(form[[2]]),form[c(1,3)]))
  return(list(form))
}

gsiLMCRasprod <- function(form) {
  lf = length(form)
  if( lf != 3 )
    return(list(form))
  if(is.call(form) ) 
    if(as.character(form[[1]])=="*" )  
      return( c(Recall(form[[2]]),Recall(form[[3]]) ) )
  return(list(form))
}



gsiLMCRas1vg <-function(form,nr="",D,envir=globalenv()) {
  if( is.call(form) ) {
    tag <- as.character(form[[1]])
    name <- paste(tag,nr,sep="")
    if( tag %in% c("sph","exp","gauss","lin","pow")) {
      func <- paste("vgram",tag,sep=".")
      stopifnot(length(form)<=2)
      range = if( length(form)>=2 ) form[[2]] else 1
      rangeVar = paste("r",name,sep="")
      params<-structure(list(range),names=rangeVar)
      .call <- as.call(list(as.name(func),h=as.name("h"),range=as.name(rangeVar)))
      list(name,params=params,call=.call,full=FALSE)
    } else if(tag %in% c("nugget")) {
      func <- paste("vgram",tag,sep=".")
      stopifnot(length(form)<=1)
      params<-list()
      .call <- as.call(list(as.name(func),h=as.name("h")))
      list(name,params=params,call=.call,full=FALSE)
    } else {
      .call <- form
      params <- as.list(form)[-1]
      list(name=name,params=params,call=.call,full=NA)
    }
  } else {
    return(NULL)
  }
}

gsiLMCRasMatTerm <- function(form,nr,D,envir=globalenv()) {
  if( is.call(form) ) { #calls
    tag <- as.character(form[[1]])
    if( tag == "function" ) { # Call and use parameters
      params <- as.list(form[[2]])
      .call <- as.call(c(list(form),lapply(names(form[[2]]),as.name)))
    } else if(tag == "PSD") { # call -> a matrix parameter is given
      sillVar <- paste("sPSD",nr,sep="")
      sill <- parameterPosdefClrMat(eval(form[[2]],envir=envir))
      params  <- structure(list(sill),names=sillVar)
      .call <- as.call(list(as.name("parametricPosdefClrMat"),as.name(sillVar)))
    } else if(tag == "R1") { # rank 1 call-> a matrix parameter is given
      sillVar <- paste("sR",nr,sep="")
      sill <- parameterPosdefClrMat(eval(form[[2]],envir=envir))
      params  <- structure(list(sill),names=sillVar)
      .call <- as.call(list(as.name("parametricRank1ClrMat"),as.name(sillVar)))
    } else if(tag == "S"){ # scalar sill factor
      sillVar <- paste("sS",nr,sep="")
      mat <- deparse(eval(form,envir=envir))
      sill = 1
      params  <- structure(list(sill),names=sillVar)
      .call <- as.call(list(as.name("("),as.call(list(as.name("*"),as.name(sillVar),mat))))
    } else { # Anything else, is used verbatim
      params <- list()
      .call <- form
    }

  } else { # values or keywords
    if( is.name(form) || (is.character(form) && length(form)==1 )) {
      tag <- as.character(form)
      if( tag=="PSD") { # keyword: allow as positive semidefinit parameter
        sillVar <- paste("sPSD",nr,sep="")
        sill    <- parameterPosdefClrMat(diag(D))
        params  <- structure(list(sill),names=sillVar)
        .call <- as.call(list(as.name("parametricPosdefClrMat"),as.name(sillVar)))
      } else if( tag=="R1" ) { # keyword allow a rank 1 parameter
        sillVar <- paste("sR",nr,sep="")
        sill    <- rep(1,D-1)
        params  <- structure(list(sill),names=sillVar)
        .call <- as.call(list(as.name("parametricRank1ClrMat"),as.name(sillVar)))
      } else { # anything else verbatim
        params <- list()
        .call <- form
      }
    } else {
      params <- list()
      .call <- form
    }
  }
  list(name="",params=params,call=.call,full =TRUE)
}

gsiLMCRasterm <- function(form,nr,D,envir=globalenv()) {
  pr <- gsiLMCRasprod(form)
  lpr <- length(pr)
  if( lpr == 1 ) {
    erg <- gsiLMCRas1vg(pr[[1]],nr,D,envir)
    if( is.na(erg$full) || erg$full )
      return(erg)
    sillVar <- paste("sPSD",nr,sep="")
    sill    <- parameterPosdefClrMat(diag(D))
    params  <- c(erg$params,structure(list(sill),names=sillVar))
    .call   <- as.call(list(as.name("%o%"),erg$call,as.call(list(as.name("parametricPosdefClrMat"),as.name(sillVar)))))
    list(name=erg$name,params=params,call=.call,full=TRUE)
  } else if(lpr == 2 ) {
    erg <- gsiLMCRas1vg(pr[[2]],nr,D,envir)
    mt  <- gsiLMCRasMatTerm(pr[[1]],nr,D,envir)
    params  <- c(erg$params,mt$params)
    .call   <- as.call(list(as.name("%o%"),erg$call,mt$call))
    list(name=erg$name,params=params,call=.call,full =TRUE)
  } else stop("Currently only simple products supported in CompLinModCoReg")
}



gsiLMCRterms <- function(formula,D,envir=environment(formula)) {
  form  <- formula[[2]] # remove ~
  lterms <-gsiLMCRassum(form)
  terms  <-mapply(gsiLMCRasterm,lterms,1:length(lterms),MoreArgs=list(D=D,envir=envir),SIMPLIFY=FALSE)
  spc <- list(names="name",params="params",calls="call",full="full")                
  erg<-c(lapply(spc,function(n) lapply(terms,function(x,n) x[[n]],n=n)),
    terms=list(terms)
    )
  erg$param<-unlist(erg$params,recursive = FALSE)
  erg$envir=envir
  erg$D=D
  .call <- erg$calls[[1]]
  if( length(erg$calls)>1 )
    for(i in 2:length(erg$calls))
      .call <- as.call(list(as.name("+"),.call,erg$calls[[i]]))
  erg$call<-.call
  erg
}

CompLinModCoReg <- function(formula,comp,D=ncol(comp),envir=environment(formula)) {
  vg <- function(h) {}
  LMCRterms <- gsiLMCRterms(formula,D,envir=envir)
  formals(vg)<- c(formals(vg),LMCRterms$param)
  body(vg)<-LMCRterms$call
  environment(vg)<-envir
  vg
}




compOKriging <- function(comp,X,Xnew,vg,err=FALSE) {

  Y <- comp
  n <- nrow(Y)
  nnew <- nrow(Xnew)
  D <- ncol(Y)
  F <- matrix(1,nrow=n)
  f <- matrix(1,nrow=nnew)
  hxx <- X[rep(1:n,n),,drop=FALSE]-X[rep(1:n,each=n),,drop=FALSE]
  hxy <- X[rep(1:n,nnew),,drop=FALSE]-Xnew[rep(1:nnew,each=n),,drop=FALSE]
  Gamma <- vg(hxx)
  gamma <- vg(hxy)
  err <- as.integer(if(err) 1 else 0)
  erg <- .C(
     gsiCGSkriging,
     zDim=gsiInt(dim(Y),2),
     z   =gsiDouble(Y),
     FDim=gsiInt(dim(F)),
     F   =gsiDouble(F),
     GammaDim=gsiInt(dim(Gamma),3),
     Gamma=gsiDouble(Gamma),
     fdim=gsiInt(dim(f),2),
     f   =gsiDouble(f),
     gammaDim=gsiInt(dim(gamma),3),
     gamma=gsiDouble(gamma),
     predDim=gsiInt(c(nnew,D),2),       
     pred =numeric(nnew*D),
     err  =if( err==1 ) numeric(nnew*D*D) else numeric(1),
     computeErrors=err)
 pred <- matrix(erg$pred,ncol=D)
 colnames(pred)<-colnames(comp)
 list(X=Xnew,
      Z=acomp(pred),
      err=if(err) aperm(structure(erg$err,dim=c(D,D,nnew)),c(3,1,2)) else NULL
      )

}



tryDebugger <- function (dump = last.dump) 
{
    debugger.look <- function(.selection) {
        for (.obj in ls(envir = dump[[.selection]], all.names = TRUE)) assign(.obj, 
            try(get(.obj, envir = dump[[.selection]])))
        cat(gettext("Browsing in the environment with call:\n   "), 
            calls[.selection], "\n", sep = "")
        rm(.obj, .selection)
        browser()
    }
    if (class(dump) != "dump.frames") {
        cat(gettext("'dump' is not an object of class 'dump.frames'\n"))
        return(invisible())
    }
    err.action <- getOption("error")
    on.exit(options(error = err.action))
    if (length(msg <- attr(dump, "error.message"))) 
        cat(gettext("Message: "), msg)
    n <- length(dump)
    calls <- names(dump)
    repeat {
        cat(gettext("Available environments had calls:\n"))
        cat(paste(1:n, ": ", calls, sep = ""), sep = "\n")
        cat(gettext("\nEnter an environment number, or 0 to exit  "))
        repeat {
            ind <- menu(as.character(calls))
            if (ind <= n) 
                break
        }
        if (ind == 0) 
            return(invisible())
        debugger.look(ind)
    }
}


