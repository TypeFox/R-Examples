#To illustrate the use of the package, the \code{SMI} dataset included in the \pkg{SDD} package and already analyzed in \citet{Bagn:DeCa:Punz:Dete:2013} and \citet{Bagn:Punz:Usin:2013} is considered.
#Data consist of $n=660$ daily returns of the Swiss Market Index spanning the period from August 12$^{\text{th}}$, 2009, to March 6$^{\text{th}}$, 2012 (the share prices used to compute the daily returns are downloadable from \url{http://finance.yahoo.com/}).
#########################
## Dependence Diagrams ##
#########################

ADF   <- function(x,                 # numeric vector related to the observed time series (or the series of residuals from a model)
                  dtype=c("ADF","CADF", "RPADF", "DeltaADF", "ACF"),         # options: ADF, CADF, RPADF, DeltaADF       
                  lag.max=floor(10 * log10(length(x))),            # the maximum lag considered
                  alpha=0.05,         # the significance level of the tests of lag independence (related to each bar)
                  num.clas,           # ADF:number of classes in each contingency table (as default, see equation (2) in the paper)
                  B=99,               # DeltaADF:number of permutation
                  bandwidth,          # DeltaADF:bandwidth
                  delta="Delta_1",    # DeltaADF:option to select the divergence measure when dtype="Delta_ADF" is used
                  fres=".Perm",       # DeltaADF:name of a function used to make permutations
                  fdenest=".denest",  # DeltaADF:name of a function used for density estimation 
                  fdiv,          # DeltaADF:name of a function used to compute divergence
                  argacf,        # List with pptional acf arguments
                  R=1:lag.max,
                  p.adjust.method =  p.adjust.methods,
                  plot=TRUE,          # Display plot
                 ...                  # plot method arguments
){
  if(missing(x)) {stop("error: no value specified for argument x")} 
  if(missing(lag.max)) 
    lag.max <- floor(10 * log10(length(x)))  # rule of thumb
  lag.max <- min(lag.max, length(x) - 1)
  if(lag.max < 0)  stop("'lag.max' must be at least 0")
  dtype <- match.arg(dtype)
  delta <- match.arg(delta,c("Delta_1", "Delta_0.5", "Delta_2", "Delta_3", "Delta_4", "Delta_SD", "Delta_L1", "Delta_ST","Delta_fdiv"), several.ok=TRUE)
  series <- deparse(substitute(x))
  result <- eval(parse(text=paste0("result <- .",dtype,".compute(x=x,lag.max=lag.max,num.clas=num.clas,alpha=alpha,B=B,dtype=dtype,bandwidth=bandwidth,fres=fres,fdenest=fdenest,delta=delta,fdiv=fdiv, argacf=argacf, plot=plot,series=series, R= R,...)")))
  class(result) <-"SDD"  
  if (plot) {
    if (is.null(list(...)$main)) plot(result, main="",...)
    else plot(result,...)
  }
  p.adjust.method <- match.arg(p.adjust.method)
  result[["R"]] <- R
  result[["p.adjust"]] <- p.adjust(result$res$pvalue[R], method = p.adjust.method)
  result[["p.adjust.method"]] <-  p.adjust.method
  invisible(result)
} 


.ACF.compute<- function (x,lag.max,num.clas,alpha,B,dtype,bandwidth,argacf,main,plot,R,series,...){
  if (missing(argacf)) {argacf <- list()}
  argacf$x <- x
  argacf$plot <- FALSE
  acf <-  do.call("acf", argacf)
  
  acf$acf = acf$acf[-1, , , drop = FALSE]
  acf$lag = acf$lag[-1, , , drop = FALSE]
  alpha <- ifelse (is.null(argacf$ci),0.05, 1 - argacf$ci)
  pvalue  <- 2 * (1 - pnorm( abs(acf$acf), 0, 1/sqrt( length(x) ) ) )
  pstar <- (pvalue<alpha)*(-1/(2*alpha)*pvalue+1)+(pvalue>=alpha)*(1-pvalue)/(2*(1-alpha))
  
  res <- data.frame(
    lag  = acf$lag,
    vbar = acf$acf, 
    pvalue   = pvalue,
    pstar    = pstar,
    crit.val = qnorm(1-alpha/2)/sqrt(length(x)),
    n = (length(x)-1):(length(x)-max(acf$lag))
  )
  acf$Portmanteau.test <- 1-pchisq(sum(res$vbar[R]^2)*length(x),df=length(R))
  acf$series <- series
  acf$res <- res
  acf$dtype    <- dtype
  acf$alpha   <- alpha
  return(acf)
 
}

################################################################################
#####  dtype ADF                                                        #########
################################################################################
.ADF.compute <- function (x,lag.max,num.clas,alpha,B,dtype,series,R,...){
  index <- length(x)
  if(missing(num.clas)){
    ks       <- floor(sqrt((index-lag.max)/5))
    kp       <- floor(sqrt(4*(2*(((index-1)^2)/(qnorm(alpha)^2)))^(1/5)))
    num.clas <- max(2,min(ks,kp))
  }
  
  if(num.clas>sqrt((index-lag.max)/5)){
    cat("warning:num.clas and/or lag.max maybe too large for the lenght of x\n")
  }
  
  tab.freq      <- NULL
  tab.freq.teo  <- NULL
  test          <- NULL
  comp.test     <- NULL
  pvalue     <- numeric(lag.max)
  ADF        <- numeric(lag.max)
  n          <- numeric(lag.max)
  
  for(i in 1:lag.max){
    
    y             <- cut2(x[(i+1):index],g=num.clas)
    z             <- cut2(x[1:(index-i)],g=num.clas)
    tab.freq[[i]] <- table(y,z)
    n[i]          <- sum(tab.freq[[i]])       
    cont.table    <- cbind(tab.freq[[i]],Total=rowSums(as.matrix(tab.freq[[i]])))
    cont.table    <- rbind(cont.table,Total=colSums(cont.table))    
    tab.freq.teo[[i]] <- outer(rowSums(as.matrix(tab.freq[[i]])),colSums(as.matrix(tab.freq[[i]])),"*")/n[i]
    comp.test[[i]]    <- (tab.freq[[i]]-tab.freq.teo[[i]])^2/tab.freq.teo[[i]]    
    test[[i]]     <- chisq.test(tab.freq[[i]])
    pvalue[i]     <- test[[i]]$p.value
    ADF[i]        <- test[[i]]$statistic
  }
  
  df         <- test[[i]]$parameter
  critvalue  <- qchisq(1-alpha,df=df)
  
  pstar <- (pvalue<alpha)*(-1/(2*alpha)*pvalue+1)+(pvalue>=alpha)*(1-pvalue)/(2*(1-alpha))
  result <- list()
  res <- data.frame(
    lag  = 1:lag.max,
    vbar = ADF, 
    pvalue   = pvalue,
    pstar    = pstar,
    crit.val = critvalue,
    n=n
    )
  Portmanteau.test <- 1-pchisq(q=sum(res$vbar[R]),df=(num.clas-1)^2*length(R))
  return(list(
      res      = res,
      dtype    = dtype,
      num.clas = num.clas,
      alpha    = alpha,
      df       = df,
      Portmanteau.test = Portmanteau.test,
      series=series
    )
    )
}                               

################################################################################
#####  dtype CADF                                                     #########
################################################################################
.CADF.compute <- function (x,lag.max,num.clas,alpha,B,dtype,series,R,...){
  ADF             <- .ADF.compute(x,lag.max,num.clas,alpha,B,dtype="ADF",series,R=R)
  CADF.ADF      <- sqrt(ADF$res$vbar/(ADF$res$n*(ADF$num.clas-1)))
  CADF.crit.val <- sqrt(ADF$res$crit.val/(ADF$res$n*(ADF$num.clas-1)))
  res <- data.frame(
    lag  = 1:lag.max,
    vbar    = CADF.ADF,
    pvalue   = ADF$res$pvalue,
    crit.val=CADF.crit.val
  )
  return(list(
    res      = res,
    dtype     = dtype,
    num.clas = ADF$num.clas,
    alpha    = alpha,
    df       = ADF$df,
    Portmanteau.test = ADF$Portmanteau.test,
    series=series
  ) )
}
################################################################################
#####  dtype RPADF                                                      #########
################################################################################
.RPADF.compute <- function (x,lag.max,num.clas,alpha,B,dtype,fres,series,...){ 

    ADF     <- .ADF.compute(x,lag.max,num.clas,alpha,B,dtype="ADF",series,R=1)
    
    lambda     <- numeric(lag.max)
    RP         <- numeric(lag.max)
    
    lambda     <-  sapply(1:lag.max,.fmin,ADF,k=1.5,g=0.5)
 
    RP         <- 1-pchisq(qchisq(1-alpha,df=ADF$df),df=ADF$df,ncp=lambda)
    
    critvalue  <- rep(qchisq(1-alpha,df=ADF$df[1]),lag.max)
       
    res <- data.frame(
      lag  = 1:lag.max,
      vbar    = RP,
      pvalue   = ADF$res$pvalue,      
      xmin = lambda
    )
    return(list(
      res      = res,
      dtype     = dtype,
      num.clas = ADF$num.clas,
      alpha    = alpha,
      df       = ADF$df,
      Portmanteau.test = ADF$Portmanteau.test,
      series=series
    ) )
}  
################################################
## Objective function to be minimized for the ##
## Estimation of the noncentrality parameter  ##
################################################
.fmin <- function(i,ADF,k,g) optimize(.fobj, c(0, k*ADF$res$vbar[i]+0.001), tol = 10^(-10),TestOss=ADF$res$vbar[i],df=ADF$df,gamma=g)$minimum
.fobj <- function(lambda,TestOss,df,gamma=0.5) (pchisq(q=TestOss,df=df,ncp=lambda)-gamma)^2
  
################################################################################
#####  dtype DeltaADF                                                   #########
################################################################################
.DeltaADF.compute <- function (x,lag.max,num.clas,alpha,B,dtype,bandwidth, fres, fdenest,delta,fdiv,series,R=1:lag.max,...)
{ 
  if (missing(bandwidth))  bandwidth <- .lcv(x)
  P         <- eval(parse(text=paste(fres,"(x,B=B)"))) #Perm(x,B=B) # (B x n)-matrix of permuted samples
  pvalue    <- numeric(lag.max)
  result    <- list()
      deltap <- sapply(1:lag.max,.p_value,x,P=P,delta=delta,fdenest=fdenest,bandwidth=bandwidth,fdiv=fdiv)
      deltapm   <- matrix(unlist(deltap),nrow=(B+2),ncol=lag.max) 
      pvalue <- deltapm[1,]
      pstar  <- (pvalue<alpha)*(-1/(2*alpha)*pvalue+1)+(pvalue>=alpha)*(1-pvalue)/(2*(1-alpha))
      Portmanteau.test <- ((rank(-rowSums(deltapm[2:(B+2),R]),ties.method="random")-1)/B)[1] #Portmanteau-test
      res    <- data.frame(
      lag  = 1:lag.max,
      vbar = pstar,
      pvalue   = pvalue
      )
    result <- list(
      res      = res,
      dtype     = dtype,
      delta     = delta,
      alpha    = alpha,
      bandwidth = bandwidth,
      Portmanteau.test = Portmanteau.test,
      series=series
    )
  return(result)
  
}

# Likelihood Cross-Validation
.lcv  <-function(x){
  
  if(sum(is.na(x)) > 0)  dati <- x[-which(is.na(x))]  # Missing data are removed 
  else  dati <- x
  
  cv <- function(x,dat=dati) {
    n  <- length(dat)
    cv <- sapply(1:n, function(d) .khleave(d,dati,exp(x)))
    cv <- 1/n*sum(log(cv))
    return(-cv)
  }
  
  iniz <- log(bw.nrd0(dati))
  st   <- nlm(cv,iniz)$estimate
  bw   <- exp(st)
  
  return(bw)
  
}

######################
## Permuted samples ##
######################

.Perm <- function(x,  # time series
                  B   # number of permutations
){
  n <- length(x)
  P <- array(0,c(B,n),dimnames=list(paste("Perm.",1:B,sep=""),1:n))
  for(b in 1:B){
    P[b,] <- sample(x,n,replace=FALSE)
  }
  return(P)
}

###################################################
## p-value from the selected divergence measure  ## 
## with the Gaussian-Kernel density estimator    ##
###################################################

.p_value <- function(lag,  # considered lag
                     x,      # time series
                     P,      # matrix of permutations arising from Perm()
                     delta,   # used divergence measure (Delta_1, Delta_0.5, Delta_2, Delta_3, Delta_4, Delta_SD, Delta_L1, Delta_ST) 
                     fdenest,  # name of a function used for density estimation 
                     bandwidth,
                     fdiv
){
 
  P           <- rbind(x,P)
  B           <- nrow(P)    # number of permutations
  Svector     <- numeric(B)
  Svector <- sapply(1:B,.stattestKERN, x=P,lag=lag,bandwidth=bandwidth,delta=delta,fdenest=fdenest, fdiv=fdiv)
  Svector1 <- (rank(-Svector,ties.method="random")-1)/B
  p.value <- Svector1[1]
  
  return(list(p.value=p.value,Svector=Svector))
  #return(p.value=p.value)
}


################################################################################
######### Divergence and Distace Measures ######################################
################################################################################
.stattestKERN<-function(j,
                        x,                   # time series
                        lag=1,               # considered lag
                        bandwidth,                # the bandwidth 
                        delta="Delta_1", # Delta_1, Delta_0.5, Delta_2, Delta_3, Delta_4, Delta_SD, Delta_L1, Delta_ST
                        fdenest=".denest",
                        fdiv
){
  x  <- x[j,]
  DE <- eval(parse(text=paste(fdenest,"(x,m=lag,bandwidth=bandwidth)"))) 
  
  fi <- DE$fi
  gi <- DE$gi
  ri <- fi/gi
  ngi <- (DE$ngi)^2  
  
  na <- function(x){return(x[which(abs(x)<Inf)])}
  
  if(delta=="Delta_1")   a<-sum(na(log(ri)*fi))/ngi
  if(delta=="Delta_0.5") a<-sum((fi^0.5-gi^0.5)^2)/ngi
  if(delta=="Delta_2")   a<-sum(na((ri-1)*fi))/ngi
  if(delta=="Delta_3")   a<-1/2*sum(na((ri^2-1)*fi))/ngi
  if(delta=="Delta_4")   a<-1/3*sum(na((ri^3-1)*fi))/ngi
  if(delta=="Delta_SD")  a<-sum((fi-gi)^2)/ngi
  if(delta=="Delta_L1")  a<-sum(abs(fi-gi))/ngi
  if(delta=="Delta_ST")  a<-sum((fi-gi)*fi)/ngi
  if(delta=="Delta_fdiv")    a<-eval(parse(text=paste(fdiv,"(fi=fi, gi=gi,ngi=ngi)"))) 
  return(a)
}
#########################################################################
####  univariate and bivariate density estimation    ####################
#########################################################################

.denest <- function(x,         # time series
                    m=1,       # considered lag
                    ngrid=100, # number of points in the grid
                    bandwidth       # bandwidth 
){
  
  n   <- length(x)
  xt  <- x[(m+1):n]
  y   <- x[1:(n-m)]
  z   <- cbind(xt,y)
  if(sum(is.na(z[,1])) > 0) z <- z[-which(is.na(z[,1])),]  # missing "row" data are removed
  if(sum(is.na(z[,2])) > 0) z <- z[-which(is.na(z[,2])),]
  ngm <- length(xt)  # n-lag
  ngi <- ngrid
  
  # Univariate Densities
  f1    <- sm.density(x, h = bandwidth, ngrid = ngi , display="none") # delete automatically missing data
  grid1 <- f1$eval.points
  
  # Joint Density
  f <- sm.density(z, h=c(bandwidth,bandwidth),eval.points=cbind(grid1,grid1),display="none")$estimate
  
  # Independence Density
  g <- f1$estimate%*%t(f1$estimate)
  
  return(list(
    fi  = f,
    gi  = g,
    ngi = ngi
  )
  )
  
}
# -- khleave --
# definition of the univariate kernel

.khleave <- function(i,x,h) {
  
  n <- length(x)
  b <- 1/((n-1)*h*sqrt(2*pi))*sum(exp(-(x[i]-x[-i])^2/(2*h^2)))
  return(b)
  
}

########################################################################################################
#####   Plot method for SDD class objects                                                      #########
########################################################################################################

plot.SDD <- function(x,
                     norm=FALSE,         # if TRUE the "normalized" p-values of the ADF are plotted 
                     stability=FALSE,    # 
                     step=5,             # the step between x-ticks (default=5) 
                     ...)
{
    eval(parse(text=paste0(".",x$dtype,".plot(x,norm=norm,stability=stability,step=step,...)")))
}  

.ADF.plot <- function(x,norm,stability,step,ylim,ylab,...){
  if(norm==TRUE){
    if (missing(ylab)) ylab <- "Transformed p-values"
    x$res[2] <- x$res$pstar
    .gen.plot(x, step=step, ylim=ylim, ylab=ylab,...)
    abline(a=1/2,b=0,col="blue",lty=2)
  }
  else {
    if (missing(ylab)) ylab <- "ADF"
    if(missing(ylim)) ylim <- c(0,max(1.2*x$res$crit.val,1.2*x$res[[2]]))
    .gen.plot(x, step=step, ylim=ylim, ylab=ylab,...)
    abline(a=x$res$crit.val,b=0,col="blue",lty=2)
  }
}
.ACF.plot <- function(x,norm,stability,step,ylim,ylab,...){
    if (missing(ylab)) ylab <- "ACF"
    if(missing(ylim)) ylim <- c(min(-1.2*x$res$crit.val,1.2*x$res[[2]]),max(1.2*x$res$crit.val,1.2*x$res[[2]]))
    .gen.plot(x, step=step, ylim=ylim, ylab=ylab,...)
    abline(a=x$res$crit.val,b=0,col="blue",lty=2)
    abline(a=-x$res$crit.val,b=0,col="blue",lty=2)
}

.CADF.plot <-  function(x,norm,stability,step,ylim,ylab,...){
  if (missing(ylab)) ylab <- "CADF"
  .gen.plot(x,ylab=ylab,step=step, ylim=ylim,...)
  points(1:max(x$res$lag), x$res$crit.val, col="blue",pch=16,cex=0.7)
}
.RPADF.plot <-  function(x,norm,stability,step,ylim,ylab,...){
  if (missing(ylab)) ylab <- "RP-ADF"
  .gen.plot(x,ylab=ylab, ylim=ylim, step,...)
  abline(a=1/2,b=0,col="blue",lty=2)
  abline(a=x$alpha,b=0,col="orange",lty=3)
  if (stability){
    lag.max <- max(x$res$lag)
    rect(1:lag.max-0.21, rep(0,lag.max), 1:lag.max+0.21, x$res$vbar, density = NULL, angle = 45, border = NULL, lty=1, lwd=2, col="gray")
  }
}
.DeltaADF.plot <-  function(x,norm,stability,step,ylim,ylab,...){
  if (missing(ylab)){
    if (x$delta=="Delta_0.5") x$delta <- "Delta_1/2"
    ylab <-parse(text=paste0("Delta[",strsplit(x$delta,"_")[[1]][2],"]-ADF")) 
  }
  .gen.plot(x,ylab=ylab,ylim=ylim, step=step,...)
  abline(a=1/2,b=0,col="blue",lty=2)
}

.gen.plot<- function(x, step,ylim,ylab,xlab,axes,type,lwd,main,...){
  if(missing(ylim)) ylim <- c(0,1)
  if(missing(xlab)) xlab <- "Lag"
  if(missing(axes)) axes <- FALSE
  if(missing(type)) type <- "h"
  if(missing(lwd)) lwd <- 1
  if(missing(main)) main <- x$series
  plot(x$res$lag,x$res[[2]],main=main,xlab=xlab,ylab=ylab,ylim=ylim, axes = axes,type = type,lwd=lwd,...)
  axis(1, at = seq(0,max(x$res$lag),by=step),labels = seq(0,max(x$res$lag),by=step)) # at = 1:lag.max, label = 1:lag.max
  axis(2)
  box(col = "black")
  abline(a=0,b=0)
}

print.SDD <- function(x, digits=3, ...)
{
  invisible(x)
  msg <- paste0(x$dtype, " bars ",ifelse(!is.null(x$delta),paste0("(with method ",x$delta,") "),""),"for series '",x$series,"'")
  cat("\n",msg,"\n")
  adf <-drop (x$res[[2]])
  names(adf) <- format(x$res[[1]])
  print(adf,digits = digits, ...)
  cat("\nSimultaneous Test \n")
  cat("  adjustment method:", x$p.adjust.method)
  cat("\n  adjusted p-value: ", round(min(x$p.adjust),3),"\n",sep="")
  cat("\nPortmanteau Test ")
  if (x$dtype %in% c("ACF")) cat("(Box-Pierce test)")
  cat("\n  p-value: ", round(x$Portmanteau.test,3),"\n",sep="")
  cat(paste0("\nTested lags \n  ",paste(x$R,collapse = ","), "\n"))
}


