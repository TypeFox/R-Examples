# periodogram class
new.periodogram <- methods::setClass("periodogram",representation(info="list"),contains="data.frame")

# periodogram wrapper
periodogram <- function(data,CTMM=NULL,T=NULL,dt=NULL,res=1,fast=NULL)
{
  # listify for general purpose code
  if(class(data)=="telemetry" || class(data)=="data.frame") { data <- list(data)  }
  if(!is.null(CTMM)) { if(class(CTMM)=="ctmm") { CTMM <- list(CTMM) } }
  
  # detrend the mean
  for(i in 1:length(data))
  {
    if(is.null(CTMM[[i]]))
    {
      data[[i]]$x <- data[[i]]$x - mean(data[[i]]$x)
      data[[i]]$y <- data[[i]]$y - mean(data[[i]]$y)
    }
    else # use the better result if provided
    {
      data[[i]]$x <- data[[i]]$x - CTMM[[i]]$mu[1] 
      data[[i]]$y <- data[[i]]$y - CTMM[[i]]$mu[2] 
    }
  }
  
  # intelligently select algorithm
  if(is.null(fast))
  {
    if(all(sapply(data,function(d){length(d$t)<10000}))) { fast <- FALSE }
    else { fast <- TRUE }
  }
  
  # default worst temporal resolution
  dts <- sapply(data,function(d) { stats::median(diff(d$t)) })
  if(is.null(dt)) { dt <- max(dts) }
  
  # default best sampling period
  Ts <- sapply(data,function(d) { last(d$t)-d$t[1] })
  if(is.null(T)) { T <- max(Ts) }
  n <- round(T/dt) + 1
  T <- (n-1)*dt
  
  if(fast)
  { lsp <- lapply(data,function(d) periodogram.fast(d,T=T,dt=dt,res=res)) }
  else
  { lsp <- lapply(data,function(d) periodogram.slow(d,T=T,dt=dt,res=res)) }

  f <- lsp[[1]]$f
  
  # DOF of one frequency in each LSP
  DOF <- (dts/dt)*2*(sapply(data,function(d){length(data$t)})-1)/length(f)

  # average the periodograms if there are multiple
  LSP <- rowSums( sapply(1:length(lsp),function(i){DOF[i]*lsp[[i]]$LSP}) )
  SSP <- rowSums( sapply(1:length(lsp),function(i){DOF[i]*lsp[[i]]$SSP}) )
    
  # DOF of one frequency in ave LSP
  DOF <- sum(DOF)
  
  LSP <- LSP/DOF
  SSP <- SSP/DOF

  result <- data.frame(LSP=LSP,SSP=SSP,DOF=DOF,f=f)
  result <- new.periodogram(result,info=mean.info(data))
  
  return(result)
}

# FFT periodogram code
periodogram.fast <- function(data,T=NULL,dt=NULL,res=1)
{
  t <- data$t
  x <- data$x
  y <- data$y
  
  # Nyquist frequency
  F <- 1/(2*dt)
  
  # construct simple time grid
  dt <- dt
  t <- t - grid.init(t,dt)
  t <- round(t/dt)
  t <- t - (t[1]-1) # observation indices
  n <- last(t)
  
  W <- rep(0,n)
  X <- rep(0,n)
  Y <- rep(0,n)

  W[t] <- 1  
  X[t] <- x
  Y[t] <- y

  # default sampling period
  if(!is.null(T)) { n <- round(T/dt) + 1 }
  # inflate padding to increase frequency resolution
  if(res>1) { n <- n*res }
  T <- (n-1)*dt # effective period
  
  # frequency resolution
  df <- 1/(2*(T+dt))
  
  # frequency grid
  f <- seq(df,F-df,df)
  
  # double padded fourier transforms
  W <- FFT(pad(W,2*n))
  X <- FFT(pad(X,2*n))
  Y <- FFT(pad(Y,2*n))

  # double frequency argument
  # this code will fail for large n
  # W2 <- rep(W,4*n)[seq(1,4*n,2)]
  # same thing but without rep
  W2 <- c(W,W)[seq(1,4*n,2)]
  W2 <- Conj(W2)
    
  # LSP and SSP denominator
  DEN <- W[1]^2 - abs(W2)^2
  
  LSP <- W[1]*(abs(X)^2+abs(Y)^2) - W2*(X^2+Y^2)
  LSP <- Re(LSP/DEN/2)
  
  SSP <- W[1]*abs(W)^2 - W2*W^2
  SSP <- Re(SSP/DEN)
  
  # Stop before Nyquist periodicity
  n <- length(f)+1
  LSP <- LSP[2:n]
  SSP <- SSP[2:n]
  
  result <- data.frame(LSP=LSP,SSP=SSP,f=f)

  #sub sample for low-resolution request
  if(res<1)
  {
    SEQ <- seq(1/res,length(f),1/res)
    SEQ <- round(SEQ)
    result <- result[SEQ,] 
  }

  return(result)
}

# slow periodogram code
periodogram.slow <- function(data,T=NULL,dt=NULL,res=1)
{
  t <- data$t
  x <- data$x
  y <- data$y

  # move some of this basic stuff into the wrapper
    
  # Nyquist frequency
  F <- 1/(2*dt)
  
  # frequency resolution
  df <- 1/(2*(T+dt)*res)
 
  # frequency grid
  f <- seq(df,F-df,df)
  
  # double angle matrix
  theta <- (4 * pi) * (f %o% t)

  # LSP lag shifts
  tau <- atan( rowSums(sin(theta)) / rowSums(cos(theta)) ) / (4*pi*f)
  
  # lagged angle matrix
  theta <- (2 * pi) * ((f %o% t) - (f * tau))

  # trigonometric matrices
  COS <- cos(theta)
  SIN <- sin(theta)
  
  LSP <- ((COS %*% x)^2 + (COS %*% y)^2)/rowSums(COS^2) + ((SIN %*% x)^2 + (SIN %*% y)^2)/rowSums(SIN^2)
  LSP <- LSP/4
  
  # sampling schedule periodogram
  SSP <- rowSums(COS)^2/rowSums(COS^2) + rowSums(SIN)^2/rowSums(SIN^2)
  SSP <- SSP/2
  
  result <- data.frame(LSP=LSP,SSP=SSP,f=f)
  
  return(result)
}


# plot periodograms
plot.periodogram <- function(x,diagnostic=FALSE,col="black",transparency=0.25,grid=TRUE,...)
{
  # frequency in 1/days
  f <- x$f*24*60^2
  
  LSP <- x$LSP
  LSP <- log(LSP/max(LSP))
  
  SSP <- x$SSP
  SSP <- log(SSP/max(SSP))
  
  at=labels=gcol=c()
  
  # tick specification function
  ticker <- function(time,div,name)
  {
    # fundamental
    at <<- c(at,time)
    labels <<- c(labels,name)
    gcol <<- c(gcol,grDevices::rgb(0.5,0.5,0.5,1))
    # harmonics
    DIV <- 2:div
    at <<- c(at,time/DIV)
    labels <<- c(labels,rep(NA,length(DIV)))
    gcol <<- c(gcol,grDevices::rgb(0.5,0.5,0.5,1/DIV))
  }
  
  # yearly periods
  ticker(365.24,11,"year")
  
  # lunar periods
  ticker(29.53059,29,"month")
  
  # diurnal periods
  ticker(1,24,"day")  

  col <- scales::alpha(col,alpha=((f[1]/f)^transparency))
  plot(1/f,LSP,log="x",xaxt="n",xlab="Period",ylab="Log Spectral Density",col=col,...)

  if(diagnostic)
  {
    col <- scales::alpha("red",alpha=((f[1]/f)^transparency))
    graphics::points(1/f,SSP,col=col,...)
  }

  graphics::axis(1,at=at,labels=labels) # tck can't be a vector apparently... :(
  if(grid){ graphics::abline(v=at,col=gcol) }
}
#methods::setMethod("plot",signature(x="periodogram",y="missing"), function(x,y,...) plot.periodogram(x,...))
#methods::setMethod("plot",signature(x="periodogram"), function(x,...) plot.periodogram(x,...))
