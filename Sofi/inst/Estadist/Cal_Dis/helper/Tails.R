normTail <-
  function(m=0, s=1, L=NULL, U=NULL, M=NULL, df=1000, curveColor=1, border=1, col='#569BBD', xlim=NULL, ylim=NULL, xlab='', ylab='', digits=2, axes=1, detail=999, xLab=c('number', 'symbol'), cex.axis=1, xAxisIncr=1, ...){
    if(is.null(xlim)[1]){
      xlim <- m + c(-1,1)*3.5*s
    }
    temp <- diff(range(xlim))
    x    <- seq(xlim[1] - temp/4, xlim[2] + temp/4, length.out=detail)
    y    <- dt((x-m)/s, df)/s
    if(is.null(ylim)[1]){
      ylim <- range(c(0,y))
    }
    plot(x, y, type='l', xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, axes=FALSE, col=curveColor, ...)
    if(!is.null(L[1])){
      these <- (x <= L)
      X <- c(x[these][1], x[these], rev(x[these])[1])
      Y <- c(0, y[these], 0)
      polygon(X, Y, border=border, col=col)
    }
    if(!is.null(U[1])){
      these <- (x >= U)
      X <- c(x[these][1], x[these], rev(x[these])[1])
      Y <- c(0, y[these], 0)
      polygon(X, Y, border=border, col=col)
    }
    if(all(!is.null(M[1:2]))){
      these <- (x >= M[1] & x <= M[2])
      X <- c(x[these][1], x[these], rev(x[these])[1])
      Y <- c(0, y[these], 0)
      polygon(X, Y, border=border, col=col)
    }
    
    if(axes == 1 || axes > 2){
      if(xLab[1]=='symbol'){
        xAt  <- m + (-3:3)*s
        xLab <- expression(mu-3*sigma, mu-2*sigma,
                           mu-sigma, mu,	mu+sigma,
                           mu+2*sigma, mu+3*sigma)
      } else if(xLab[1] != 'number'){
        stop('Argument "xLab" not recognized.\n')
      } else {
        temp <- seq(xAxisIncr, max(abs(xlim-m))/s, xAxisIncr)*s
        xAt <- m + c(-temp, 0, temp)
        xLab <- round(xAt, digits=digits)
      }
    }
    if(axes > 2){
      axis(1, at=xAt, labels=xLab, cex.axis=cex.axis)
      buildAxis(2, c(y,0), n=3, nMax=3, cex.axis=cex.axis)
    } else if(axes > 1){
      buildAxis(2, c(y,0), n=3, nMax=3, cex.axis=cex.axis)
    } else if(axes > 0){
      axis(1, at=xAt, labels=xLab, cex.axis=cex.axis)
    }
    
    abline(h=0)
  }

chiTail <-
function(U=NULL, df = 10, curveColor=1, border=1, col="#569BBD", xlim=NULL, ylim=NULL, xlab='', ylab='', detail=999){
	#if(U <= 30){xlim <- c(0,30)}
	#if(U > 30){xlim <- c(0,U+0.01*U)}
	temp <- diff(range(xlim))
	x    <- seq(xlim[1] - temp/4, xlim[2] + temp/4, length.out=detail)
	y    <- dchisq(x, df)
	ylim <- range(c(0,y))
	plot(x, y, type='l', xlim=xlim, ylim=ylim, axes=FALSE, col=curveColor, xlab = "", ylab = "")
	these <- (x >= U)
	X <- c(x[these][1], x[these], rev(x[these])[1])
	Y <- c(0, y[these], 0)
	polygon(X, Y, border=border, col=col)
	abline(h=0)
	axis(1, at = c(0,U), label = c(NA,round(U,4)))
}

FTail <-
  function(U=NULL, df_n=100, df_d = 100, curveColor=1, border=1, col="#569BBD", xlim=NULL, ylim=NULL, xlab='', ylab='', detail=999){
    if(U <= 5){xlim <- c(0,5)}
    if(U > 5){xlim <- c(0,U+0.01*U)}
    temp <- diff(range(xlim))
    x    <- seq(xlim[1] - temp/4, xlim[2] + temp/4, length.out=detail)
    y    <- df(x, df_n, df_d)
    ylim <- range(c(0,y))
    plot(x, y, type='l', xlim=xlim, ylim=ylim, axes=FALSE, col=curveColor, xlab = "", ylab = "")
    these <- (x >= U)
    X <- c(x[these][1], x[these], rev(x[these])[1])
    Y <- c(0, y[these], 0)
    polygon(X, Y, border=border, col=col)
    abline(h=0)
    axis(1, at = c(0,U), label = c(NA,round(U,4)))
  }

Valor_Fin<-function(dist_CalDis, tail_CalDis, mu_CalDis, sd_CalDis, a_CalDis, b_CalDis, df_CalDis, df1_CalDis, df2_CalDis, 
                    n_CalDis, p_CalDis, lower_bound_CalDis, upper_bound_CalDis){
    if (is.null(tail_CalDis) | is.null(a_CalDis))
    {
      shiny:::flushReact()
      return()
    }
    
    L = a_CalDis
    U = NULL
    
    if (tail_CalDis %in% c("both","middle")) 
    {
      if (is.null(b_CalDis))
      {
        shiny:::flushReact()
        return()
      }
      
      U = b_CalDis
      
      error = FALSE
      if (L>U) error = TRUE
      if (error){
        return()
      }
    }
    
    f = function() NULL
    
    if (dist_CalDis == "rnorm")
    {
      if (is.null(mu_CalDis) | is.null(sd_CalDis))
      {
        shiny:::flushReact()
        return()
      }
      
      f = function(x) pnorm(x,mu_CalDis,sd_CalDis)
    }  
    else if (dist_CalDis == "rt")
    {
      if (is.null(df_CalDis))
      {
        shiny:::flushReact()
        return()
      }
      
      f = function(x) pt(x,df_CalDis)
    }
    else if (dist_CalDis == "rchisq"){
      if (is.null(df_CalDis))
      {
        shiny:::flushReact()
        return()
      }
      
      f = function(x) pchisq(x,df_CalDis)
    }
    else if (dist_CalDis == "rf"){
      if (is.null(df1_CalDis) | is.null(df2_CalDis))
      {
        shiny:::flushReact()
        return()
      }
      
      f = function(x) pf(x,df1_CalDis,df2_CalDis)
    }    
    else if (dist_CalDis == "rbinom")
    {
      if (is.null(n_CalDis) | is.null(p_CalDis) | is.null(lower_bound_CalDis))
      {
        shiny:::flushReact()
        return()
      }
      
      if (tail_CalDis == "equal")
      {
        f = function(x) dbinom(x,n_CalDis,p_CalDis)
      }
      else
      {
        f = function(x) pbinom(x,n_CalDis,p_CalDis)
        
        if (tail_CalDis %in% c("lower","both") & lower_bound_CalDis == "open") L = L-1
        if (tail_CalDis %in% c("upper")        & lower_bound_CalDis == "closed") L = L-1
        if (tail_CalDis %in% c("middle")       & lower_bound_CalDis == "closed") L = L-1
        
        if (tail_CalDis %in% c("both","middle")) 
        {
          if (is.null(upper_bound_CalDis))
          {
            shiny:::flushReact()
            return()
          }
          
          if (tail_CalDis == "both"   & upper_bound_CalDis == "closed") U = U-1
          if (tail_CalDis == "middle" & upper_bound_CalDis == "open") U = U-1
        } 
      }
    }
    
    val = NA
    if (tail_CalDis == "lower")
      val = f(L)
    else if (tail_CalDis == "upper")
      val = 1-f(L)
    else if (tail_CalDis == "equal")
      val = f(L)
    else if (tail_CalDis == "both")
      val = f(L) + (1-f(U))
    else if (tail_CalDis == "middle")
      val = f(U) - f(L)
    return(val)
}