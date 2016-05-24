# Checked at 10-08-12 MS
# 
# This is an smoothing cubic spline function 
# originaly created  Monday, April 28, 2008 at 08:42 
# Mikis Stasinopoulos
# TO DO 
# 1) all argument apart from df and spar should be put under control OK
# 2) Can we use ML? the answer is probably no sinse we do not have access to lambda
# for example while sigma_e can be calulated as
# sum((y-fitted(fit))^2)/ (length(y)-fit$df)
# and sigma_b as
# (fit$fit$coef^2)/m1$df
# the ratio is not input for smooth.spline() but spar
# spar is lambda = r * 256^(3*spar - 1) where r is in the code but not in the output
# 3) cs() can be use the same interface and fix df's YES
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
# this is the scs() function
scs <-function(x, df = NULL, spar = NULL, 
        control=cs.control(...),...) 
{   
  scall <- deparse(sys.call()) 
     if(ncol(as.matrix(x)) > 1)
    stop(paste("The default smoother is bivariate; you gave a matrix as an argument in ", scall, "\n"))
   # len <- if(is.null(dim(x))) length(x) else dim(x)[1]
    if(!is.null(levels(x))) 
      {
        if(inherits(x, "ordered"))
            x <- as.numeric(x)
        else stop("unordered factors cannot be used as smoothing variables")
      }
                     a <- is.na(x)
            real.call  <- substitute(gamlss.cs(data[[scall]], z, w, spar = spar, df = df))
     attr(x,"control") <- control
       attr(x, "call") <- real.call
      attr(x, "class") <- "smooth"
        if(any(a))
        attr(x, "NAs") <- seq(along = x)[a]
    x
}
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
# this is the cs() function
cs <-function(x, df = 3, spar = NULL,  c.spar=NULL,
        control=cs.control(...),...)
{   
  scall <- deparse(sys.call()) 
     if(ncol(as.matrix(x)) > 1)
    stop(paste("The default smoother is bivariate; you gave a matrix as an argument in ", scall, "\n"))
   # len <- if(is.null(dim(x))) length(x) else dim(x)[1]
    if(!is.null(levels(x))) 
      {
        if(inherits(x, "ordered"))
            x <- as.numeric(x)
        else stop("unordered factors cannot be used as smoothing variables")
      }
      if(is.null(c.spar)) control$control.spar <- list(low=-1.5,high=2) 
    else
      { control$control.spar <- if (is.list(c.spar) & length(c.spar)==2)  
                      {  list(low=c.spar[[1]],high=c.spar[[2]])}
                  else if (is.vector(c.spar) & length(c.spar)==2)   
                      {  list(low=c.spar[1],high=c.spar[2])}
        else stop("c.spar is not defined properly") 
      }   
    if(!is.null(df)&&df < 0) {df <- 3; warning("the df are set to 3")}
        df <- if (is.null(df))   NULL
          else df
                     a <- is.na(x)
            real.call  <- substitute(gamlss.cs(data[[scall]], z, w, spar = spar, df = df))
     attr(x,"control") <- control
       attr(x, "call") <- real.call
      attr(x, "class") <- "smooth"
        if(any(a))
        attr(x, "NAs") <- seq(along = x)[a]
    x
}
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
# control for cs() and scs() functions
cs.control <- function( cv = FALSE, all.knots = TRUE, nknots = NULL,
                       keep.data = TRUE, df.offset = 0, penalty = 1.4,# note that his is 1 in new 
                       control.spar = list(),  ...) # smooth splines() 
{                   
        if(df.offset < 0) {
warning("the value of df.offset supplied is negative the default value of 0 was used instead")
                df.offset <- 0} 
         if(penalty < 0) {
warning("the value of penalty supplied is negative the default value of 1 was used instead")
                penalty <- 1}                   
        list( cv = cv, all.knots = all.knots, nknots = nknots,
              keep.data = keep.data, df.offset = df.offset, penalty = penalty,
              control.spar = control.spar)#
}
#----------------------------------------------------------------------------------------
# the definition of the backfitting additive function
# fitting cubic splines
gamlss.cs <-function(x, y, w, df = NULL, spar = NULL, xeval = NULL, ...)
{      
             x <- signif(x, 6)
           pox <- order(x)
          freq <- table(x)
            df <- if (!is.null(df))  df+2
       control <- as.list(attr(x, "control")) 
    if (is.null(df)&&is.null(spar))
      {
       fit <- smooth.spline(y=y, x=x, w=w, 
                             cv= control$cv,
                             all.knots=control$all.knots, nknots=control$nknots, 
                             control.spar=control$control.spar,  penalty = control$penalty )
      }
    else if (is.null(df))
      {  fit <- smooth.spline(y=y, x=x, w=w, spar=spar, 
                              all.knots=control$all.knots, nknots=control$nknots, 
                              control.spar=control$control.spar, penalty = control$penalty)
      }
     else 
     {  fit <- smooth.spline(y=y, x=x, w=w, df=df, 
                             all.knots=control$all.knots, nknots=control$nknots, 
                             control.spar=control$control.spar,  penalty = control$penalty)
      }
       longfv  <- rep(fit$y,freq)
   longfv[pox] <- longfv # OK
          llev <- rep(fit$lev,freq) # the leverage   calculations    
     llev[pox] <- llev           
          sumw <- rep(fit$w,freq)     
     sumw[pox] <- sumw
            wt <-  (w * sum(w > 0))/sum(w)           
       longlev <- llev*wt/sumw        
          levl <- (longlev-.hat.WX(wt,x))
           var <- levl/w # MS Tuesday, June 22, 2004 at 20:58
if (is.null(xeval)) # if no prediction  
    {
       obj.out <- list(residuals=y-longfv, fitted.values = longfv, 
                   var = var,  nl.df = fit$df-2, lambda = fit$lambda, 
                   coefSmo = list(knot = fit$fit[["knot"]], #
                                    nk = fit$fit[["nk"]], 
                                   min = fit$fit[["ux[1]"]], 
                                 range = fit$fit[["r.ux"]], 
                                  coef = fit$fit[["coef"]], 
                              pen.crit = fit$pen.crit, 
                                   lev = levl,
                   lambda1 = fit$lambda))
class(obj.out) <- "smooth.spline"
    obj.out
    }
else 
   { # if  prediction  
                pred <- predict(fit,x = xeval)
                pred$y  
    }    
}







 
