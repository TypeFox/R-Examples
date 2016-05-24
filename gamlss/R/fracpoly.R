# this is creates a base for given power 
bfp <-  function(x, powers=c(1,2), shift = NULL, scale = NULL)
 {
 #--------
 fp.scale <-function(x)
  {
    if(min(x) <= 0) 
    {        z <- sort(x)[-1] - sort(x)[ - length(x)]
         shift <- min(z[z > 0]) - min(x)
    }
    else shift <- 0
    range <- max(x) - min(x)
    scale <- 10^(sign(log10(range)) * trunc(abs(log10(range))))
    list(shift=shift, scale=scale)
  }
 #------- 
    nobs  <- length(x)
    npoly <- length(powers)  
    X <- matrix(0, nrow = nobs, ncol = npoly)   #
 #    if any is not defined
        if(is.null(scale)|is.null(shift)) 
         {
            out <- fp.scale(x)
          shift <- out$shift
          scale <- out$scale
         }
        x <- x + shift
        x <- x/scale
        x1 <- ifelse(powers[1] != rep(0, nobs), x^powers[1], log(x))
        X[, 1] <- x1
if (npoly >= 2)
 {
     for(i in 2:npoly)
      {
        if(powers[i] == powers[(i-1)])  x2 <- log(x) * x1
        else x2 <- ifelse(powers[i] != rep(0,nobs), x^powers[i], log(x))
        X[,i] <- x2
        x1 <- x2
      }
 }
   X
}


#----------------------------------------
# define the  formula function
fp <-function(x, npoly=2, shift = NULL, scale = NULL) {
  fp.scale <-function(x)
   {
    if(min(x) <= 0) 
    {        z <- sort(x)[-1] - sort(x)[ - length(x)]
         shift <- min(z[z > 0]) - min(x)
    }
    else shift <- 0
    range <- max(x) - min(x)
    scale <- 10^(sign(log10(range)) * trunc(abs(log10(range))))
    list(shift=shift, scale=scale)
   }
   #----
    scall <- deparse(sys.call())
    if(ncol(as.matrix(x)) > 1)
      stop(paste("The default smoother is bivariate; you gave a matrix as an argument in ",
            scall, "\n"))
    if(!is.null(levels(x))) 
        {
        if(inherits(x, "ordered"))
            x <- as.numeric(x)
        else stop("unordered factors cannot be used as smoothing variables")
        }
    if(!any(npoly %in%c(1,2,3))) stop("the number of fractional polynomials can be 1,2 or 3")
    if(is.null(scale)|is.null(shift)) 
        {
           out <- fp.scale(x)
         shift <- out$shift
         scale <- out$scale
        }
    if(shift > 0) warning("The origing of the x variable has been shifted by ", shift, "\n" )
    x <- x + shift
    x <- x/scale 
    {
     if(npoly==1) df <- 2 
     if(npoly==2) df <- 4 
     if(npoly==3) df <- 6 
    }  
    xvar <- rep(0, length(x))
    attr(xvar, "values") <- x
    attr(xvar, "npoly") <- npoly
    attr(xvar, "shift") <- shift
    attr(xvar, "scale") <- scale
    attr(xvar, "df") <- df
    real.call <- substitute(gamlss.fp(data[[scall]], z, w, npoly))
    attr(xvar, "call") <- real.call
    attr(xvar, "class") <- "smooth"
    a <- is.na(x)
    if(any(a))
    attr(xvar, "NAs") <- seq(along = x)[a]
    xvar
}
# -------------------------------
# define the backfitting function
# -------------------------------
gamlss.fp <-function(x, y, w, npoly = 2, xeval = NULL )
{
    bfp <-  function(x, powers=c(1,2))
  { nobs  <- length(x)
    npoly <- length(powers)  
    X <- matrix(0, nrow = nobs, ncol = npoly)   
        x1 <- ifelse(powers[1] != rep(0, nobs), x^powers[1], log(x))
        X[, 1] <- x1
        if (npoly >= 2)
       {
           for(i in 2:npoly)
          {
          if(powers[i] == powers[(i-1)])  x2 <- log(x) * x1
          else x2 <- ifelse(powers[i] != rep(0,nobs), x^powers[i], log(x))
          X[,i] <- x2
          x1 <- x2
          }
       }
   X
  }
   xvar <-  if (is.null(xeval)) as.vector(attr(x,"values"))#ms Wednesday, July 21, for prediction
        else  as.vector(attr(x,"values"))[seq(1,length(y))]
        npoly <- as.vector(attr(x,"npoly"))
          df  <- as.vector(attr(x,"df"))
       shift  <- as.vector(attr(x,"shift"))
       scale  <- as.vector(attr(x,"scale")) 
    if(!any(npoly %in%c(1,2,3))) stop("the number of fractional polynomials can be 1, 2 or 3")
       powers <- c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)
      npowers <- length(powers)
       devold <-  1000000000000000    #
 if (npoly==1)
 {
   for(i in 1:npowers) 
   {
      x.fp <- bfp(xvar, powers[i])
      one  <- rep(1,length(xvar))
      x.fp <- cbind(one,x.fp)
       fit <- lm.wfit(x=x.fp,y=y,w=w,method="qr")
      residuals <- fit$residuals
            dev <- sum(w*residuals^2)
            if(dev < devold) 
              {
                   devold <- dev
                powerbest <- powers[i] 
              }
   } 
 }
if (npoly==2)  
 {
 for(i in 1:npowers) 
 {       
           j <- i 
            while(j <= npowers) 
                    {                           
         x.fp <- bfp(xvar, c(powers[i], powers[j]))
         one  <- rep(1,length(xvar))
         x.fp <- cbind(one,x.fp)
         fit <- lm.wfit(x=x.fp,y=y,w=w,method="qr") # 
    residuals <- fit$residuals
          dev <- sum(w*residuals^2)  
                   if(dev < devold) 
                     {
                     devold <- dev
                  powerbest <- c(powers[i], powers[j])
                     }
                  j <- j + 1
                    }
  }
 }
if (npoly==3)
 {
 for(i in 1:npowers) 
  {             
           j <- i 
                  while(j <= npowers) 
                    {  
                     k <- j                  
                               while(k <= npowers) 
                            {    
                        x.fp <- bfp(xvar, c(powers[i], powers[j], powers[k] ))
                        one  <- rep(1,length(xvar))
                        x.fp <- cbind(one,x.fp)
                         fit <- lm.wfit(x=x.fp,y=y,w=w,method="qr")# 
                   residuals <- fit$residuals
                         dev <- sum(w*residuals^2)
                   if(dev < devold) 
                     {
                      devold <- dev
                   powerbest <- c(powers[i], powers[j], powers[k])
                     }                
                          k <- k + 1
                            }
                  j <- j + 1
                    }
  }
 }
## final fit
    power <- powerbest
     x.fp <- bfp(xvar, c(power))
    # one  <- rep(1,length(xvar))
    # x.fp <- cbind(one,x.fp)
      fit <- lm(y~x.fp,weights=w)# 
residuals <- fit$residuals
      cov <-  diag(chol2inv(fit$qr$qr[1:fit$rank, 1:fit$rank, drop = FALSE]))
      lev <- (hat(fit$qr) -.hat.WX(w,rep(1,length(w))))
      var <- lev/w # MS Tuesday, June 22, 2004 at 20:58
fit$power <- power
#fit$varcoeff <- cov
if (is.null(xeval)) # if no prediction  
 {
list(x = xvar, fitted.values = fitted(fit), residuals = resid(fit),
     var = var, nl.df = df, lambda = power, 
     coefSmo = fit)    
 }
else 
 {# if prediction
    ll <- length(longx <- as.vector(attr(x,"values")))
 nxvar <- as.vector(attr(x,"values"))[seq(length(y)+1,ll)]
 nx.fp <- bfp(nxvar, c(power))
 none  <- rep(1,length(nxvar))
 nx.fp <- cbind(none,nx.fp)
  pred <- drop(nx.fp %*% coef(fit))
 pred
 }
}
