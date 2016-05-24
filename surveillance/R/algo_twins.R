######################################################################
# Experimental version -- integrating the twins program into
# the surveillance package
######################################################################

algo.twins <- function(disProgObj,
                       control=
                       list(burnin=1000, filter=10,
                            sampleSize=2500, noOfHarmonics=1,
                            alpha_xi=10, beta_xi=10, psiRWSigma=0.25,
                            alpha_psi=1, beta_psi=0.1, nu_trend=FALSE,
                            logFile="twins.log"))
{
  if (inherits(disProgObj, "sts"))
      disProgObj <- sts2disProg(disProgObj)

  if (ncol(disProgObj$observed)>1) {
    stop("algo.twins() only handles univariate time series of counts")
  }
  
  ## Determine period from data
  T <- as.integer(disProgObj$freq)

  ## set default values (if not provided in control)
  if(is.null(control[["burnin",exact=TRUE]]))
    control$burnin <- 1000
  if(is.null(control[["filter",exact=TRUE]]))
    control$filter <- 10
  if(is.null(control[["sampleSize",exact=TRUE]]))
    control$sampleSize <- 2500
  if(is.null(control[["alpha_xi",exact=TRUE]]))
    control$alpha_xi <- 10
  if(is.null(control[["beta_xi",exact=TRUE]]))
    control$beta_xi <- 10
  if(is.null(control[["psiRWSigma",exact=TRUE]]))
    control$psiRWSigma <- 0.25
  if(is.null(control[["alpha_psi",exact=TRUE]]))
    control$alpha_psi <- 1
  if(is.null(control[["beta_psi",exact=TRUE]]))
    control$beta_psi <- 0.1
  if(is.null(control[["nu_trend",exact=TRUE]]))
    control$nu_trend <- FALSE
  if(is.null(control[["logFile",exact=TRUE]]))
    control$logFile <- "twins.log"
  if(is.null(control[["noOfHarmonics",exact=TRUE]]))
    control$noOfHarmonics <- 1
  

  nfreq <- control$noOfHarmonics
  control$logFile2 <- paste(control$logFile,"2",sep="")

  ## Call the C code
  x <- disProgObj$observed
  n <- as.integer(dim(x)[1])
  I <- as.integer(dim(x)[2])

  ## with(control, res <- .C(...)) is not valid R syntax!!
  res <- 
      with(control,
           .C("twins",
              x=as.integer(x),
              n=n,
              I=I,
              logFile=logFile,
              logFile2=logFile2,
              burnin=as.integer(burnin), 
              filter=as.integer(filter),
              sampleSize=as.integer(sampleSize),
              alpha_xi=as.double(alpha_xi),
              beta_xi=as.double(beta_xi),
              T=as.integer(T),
              nfreq=as.integer(nfreq),
              psiRWSigma=as.double(0.25),
              alpha_psi=as.double(alpha_psi),
              beta_psi=as.double(beta_psi),
              nu_trend=as.integer(nu_trend),
              PACKAGE="surveillance"))

  ## Log files
  results <- read.table(control$logFile,header=T,na.strings=c("NaN","-NaN"))
  results2 <- read.table(control$logFile2,header=T,na.strings=c("NaN","-NaN"))
  acc <- read.table(paste(control$logFile,".acc",sep=""),col.names=c("name","RWSigma","acc"))

  rownames(acc) <- acc[,1]
  acc <- acc[,-1]
  
  ## Nothing is returned by the function - result is not a
  ## standard survObj
  result <- structure(list(control=control,
                           disProgObj=disProgObj,
                           logFile=results,
                           logFile2=results2),
                      class="atwins")
  return(result)
}




######################################################################
# Adapted the functions form figures.R
######################################################################

## Helper functions to make list of Z and the means of X,Y and omega
make.pois <- function(obj) {
    n <- nrow(obj$disProgObj$observed)
    m<-list()
    m$n <- n
    m$Z <- obj$disProgObj$observed
    m$X <- numeric(n)
    m$Y <- numeric(n)
    m$omega <- numeric(n)
    ## Read means at each time instance
    Vars <- c("X","Y","omega") 
    for (t in 1:n) {
        for (v in Vars) {
            m[[v]][t] <- obj$logFile2[,paste(v,".",t,".",sep="")]
        }
    }
    return(m)
}

pois.plot <- function(m.results,...) {
    plotorder <- c(expression(Z),expression(Y),expression(X))
    plotcols <- c(1,"red","blue")
    lwd <- c(1,3,3)
    sts <- disProg2sts(m.results$disProgObj)

    ## Make default legend if nothing else is specified.
    if (!"legend.opts" %in% names(list(...))) {
        plot(sts,legend.opts=list(x="topleft",legend=paste(plotorder),lwd=lwd,col=plotcols,horiz=TRUE,y.intersp=0,lty=1,pch=NA),...)
    } else {
        plot(sts,...)
    }
    
    ## Add Y and X lines
    for (i in 2:length(plotorder)) {
        lines(1:(m.results$n)+0.5,m.results[[paste(plotorder[i])]][c(2:m.results$n,m.results$n)],type="s",col=plotcols[i],lwd=lwd[i])
    }
}

## makes list of gamma, zeta and nu
make.nu <- function(obj) {
    n <- nrow(obj$disProgObj$observed)
    samplesize <- obj$control$sampleSize
    frequencies <- obj$control$noOfHarmonics # instead of just always "1" !
    
    season <- obj$disProgObj$freq
    basefrequency <- 2 * pi / season
    
    ## optionally also get the linear time trend coefficient
    withTrend <- obj$control$nu_trend

    ## this list will be returned at the end
    m<-list()        
    
    ## first get all the gamma's from the logFile matrix into nicer elements of
    ## the list m
    for (j in 0:(2*frequencies + withTrend)) {
        m$gamma[[j+1]] <- numeric(samplesize)
        m[["gamma"]][[j+1]] <- obj$logFile[,paste("gamma",".",j,".",sep="")]
    }

    ## zeta is a list which has one element for each time point (vector of samples)
    m$zeta<-list()
    
    ## for all time points:
    for (t in 1:n) {

        ## start with the intercept
        m$zeta[[t]]<-m$gamma[[1]]

        ## add all harmonic terms
        for(j in 1:frequencies){
            m$zeta[[t]] <- m$zeta[[t]] + m$gamma[[2*j]]*sin(basefrequency*j*(t-1)) + m$gamma[[2*j+1]]*cos(basefrequency*j*(t-1)) 
        }
        
        ## and (optionally) finally add the linear trend
        if(withTrend)
        {
            m$zeta[[t]] <- m$zeta[[t]] + m$gamma[[2*frequencies + 2]] * (t - n/2)
        }
    }

    ## nu is the analogous list with the exponentiated zeta's
    m$nu<-list()
    for (t in 1:n) {
        m$nu[[t]]<-exp(m$zeta[[t]])
    }

    ## also copy the number of harmonics
    m$frequencies <- frequencies

    ## and return
    return(m)
}

## Function to plot median, and quantiles over time for m.par (m.par is list of n vectors, x is time)
tms.plot <-function(x,m.par,xlab="",ylab="",ylim=FALSE,...){
    m<-list()
    n<-length(m.par)
    m$median<-numeric(n)
    for (t in 1:n) {
        m$median[t]<- median(m.par[[t]])
        m$q025[t]<- quantile(m.par[[t]],0.025)
        m$q975[t]<- quantile(m.par[[t]],0.975)
    }
    if(!ylim){
        ymin<-min(m$q025)
        ymax<-max(m$q975)
        ylim=c(ymin,ymax)
    }

    plot(x-1,m$q975[x],type="l",col="red",main="",xlab=xlab,ylab=ylab,ylim=ylim,...) 
    lines(x-1,m$median[x],type="l")
    lines(x-1,m$q025[x],type="l",col="red")
}

######################################################################
# Function to plot an atwins object -- currently not
# properly documented
######################################################################

plot.atwins <- function(x, which=c(1,4,6,7), ask=TRUE,...) {
    
    ## Extract from the 3 dots
    if(is.null(which)) {
        which <- c(1,4,6,7)
    } 
    if(is.null(ask)) {
        ask <- TRUE
    }
    
    ## Make list of X,Y,Z,omega means of results2
    m.results <-make.pois(x)
    m.results$disProgObj <- x$disProgObj
    
    ## Make list of results of  gamma, zeta and nu
    nu<-make.nu(x)

    
    ## Plots
    show <- rep(FALSE,7)
    show[which] <- TRUE
    par(ask=ask)
    
    if (show[1]) {
        par(mfcol=c(1,1))
        pois.plot(m.results,...)
    }

    if (show[2]) {
        ## make room for 2 * (frequencies + 1) panels
        par(mfcol=c(2,nu$frequencies+1))

        ## and plot all gamma coefficients (possibly including the linear time
        ## trend coef)
        for(j in seq_along(nu$gamma)) {
            plot(nu$gamma[[j]],type="l",ylab=paste("gamma",j - 1,sep=""))
        }
    }

    if (show[3]) {
        par(mfcol=c(1,1))
        plot(x$logFile$K,type="l",ylab=expression(K))
        plot(x$logFile$xilambda,type="l",ylab=expression(xi))
        plot(x$logFile$psi,type="l",ylab=expression(psi))
    }

    if (show[4]) {
        par(mfcol=c(1,2))
        acf(x$logFile$K,lag.max = 500,main="",xlab=expression(K))
        acf(x$logFile$psi,lag.max = 500,main="",xlab=expression(psi))
    }

    if (show[5]) {
        par(mfcol=c(1,1))
        tms.plot(2:m.results$n,nu$nu,xlab="time")
    }

    if (show[6]) {
        par(mfcol=c(1,2))
        hist(x$logFile$K,main="",xlab=expression(K),prob=TRUE,breaks=seq(-0.5,max(x$logFile$K)+0.5,1))
        hist(x$logFile$psi,main="",xlab=expression(psi),prob=TRUE,nclass=50)
    }

    if (show[7]) {
        par(mfcol=c(1,1))
        hist(x$logFile$Znp1,main="",xlab=expression(Z[n+1]),prob=TRUE,breaks=seq(-0.5,max(x$logFile$Znp1)+0.5,1))
    }
}

