#----------------------------------------------------------------------------------------
mrf<-function(x, precision=NULL, neighbour=NULL, polys= NULL, area=NULL, adj.weight=1000, control=mrf.control(...), ...)  
{
# ---------------------------------------------------
            scall <- deparse(sys.call(), width.cutoff = 500L)
#----------------------------------------------------
            rexpr <- regexpr("gamlss",sys.calls())
          for (i in length(rexpr):1)
          { 
          position <- i 
           if (rexpr[i]==1) break
          }
        gamlss.env <- sys.frame(position)
#if (is.null(precision)) stop("the function needs the penalty matrix")
if (!is(x, "factor")) stop("x must be a factor")
k <- area
if (is.null(k)) k <- factor(levels(x),levels=levels(x)) # default knots = all regions are in 
else{
  if (class(area)=="character") k <- as.factor(k)
  if (!(class(k)=="character"||class(k)=="factor")) 
    stop("area must be a factor or a chacacter vector")
}
if (length(levels(x))>length(levels(k))) 
  stop("MRF basis dimension set too high")
if (sum(!levels(x)%in%levels(k))) 
  stop("data contain regions that are not contained in the area specification")
x <- factor(x,levels=levels(k))
#X <- model.matrix(~x-1,) 
nlevx <- nlevels(x)
if (is.null(precision)&&is.null(neighbour)&&is.null(polys))
  stop("precision matrix, boundary polygons and/or neighbours list must be supplied")
if (!is.null(precision))
{ 
  if (!is.matrix(precision)||dim(precision)[1]!=nlevx||dim(precision)[2]!=nlevx) 
    stop("the precision matrix is not suitable")
  precision <- precision 
} 
# check the precision matrix
if (!is.null(neighbour)&&is.null(precision))
{ # if neighbour exits then calculate the precision
  precision <- nb2prec(neighbour=neighbour,x=x,area=area)
}
if (!is.null(polys)&&is.null(neighbour)&&is.null(precision))  
{ # if polys exits then calculate the precision
  a.name <- names(polys)
  d.name <- unique(a.name[duplicated(a.name)])
  if (length(d.name)) 
  {
    for (i in 1:length(d.name)) 
    {
      ind <- (1:length(a.name))[a.name == d.name[i]]
      for (j in 2:length(ind)) 
        polys[[ind[1]]] <- rbind(polys[[ind[1]]], c(NA, NA), polys[[ind[j]]])
    }
    #now delete the un-wanted duplicates
    ind <- (1:length(a.name))[duplicated(a.name)]
    if (length(ind) > 0) 
      for (i in length(ind):1) polys[[ind[i]]] <- NULL
  }#polygon list in correct format
  neighbour <- polys2nb(polys)        ####??
  precision   <- nb2prec(neighbour=neighbour,x=x,area=area) ###??
}
#--------
## get a random name to use it in the gamlss() environment
#--------
                sl <- sample(letters, 4)
       fourLetters <- paste(paste(paste(sl[1], sl[2], sep=""), sl[3], sep=""),sl[4], sep="")
   startLambdaName <- paste("start.Lambda",fourLetters, sep=".")
## put the starting values in the gamlss()environment
      weights.Adj  <- FALSE
            sigmas <- c(control$sig2e, control$sig2b, weights.Adj)
            assign(startLambdaName, sigmas, envir = gamlss.env)
#--------
                             xvar <- rep(0,length(x)) #  set x to zero
                           #  xvar <- seq(1,length(x))
      attr(xvar, "control")       <- control
      attr(xvar, "call")          <- substitute(gamlss.mrf(data[[scall]], z, w)) 
      attr(xvar, "gamlss.env")    <- gamlss.env
      attr(xvar, "x")             <- x
      attr(xvar, "adj.weight")    <- adj.weight
      attr(xvar, "area")          <- area
      attr(xvar, "precision")     <- precision
      attr(xvar, "NameForLambda") <- startLambdaName
      attr(xvar, "class")         <- "smooth"
      xvar
}
#----------------------------------------------------------------------------------------
##---------------------------------------------------------------------------------------
mrf.control <- function(sig2e=1, sig2b=1,  
                        sig2e.fix=FALSE, sig2b.fix = FALSE, 
                        penalty = FALSE, 
                        delta = c(0.01, 0.01), 
                        shift = c(0, 0), ...) # 
{ 
                        
         if(sig2e < 0) {
warning("the value of sig2e supplied is zero or negative the default value of 1 was used instead")
                sig2e <- 1} 
         if(sig2b < 0) {
warning("the value of sig2b supplied is zero or negative the default value of 1 was used instead")
                sig2b <- 1}    
          if(any(delta <= 0)||length(delta)!=2)
          {
warning("delta should be positive and of length 2, deltat is set at default values c(0.01,0.01)")           
         delta <- c(0.01, 0.01) 	
          }
         if(length(shift)!=2)
          {
warning("shift  length should be 2, is it set at default values c(0,0)")           
         shift <- c(0, 0) 	
          }
                        
        list( sig2e=sig2e, sig2b=sig2b, sig2e.fix=sig2e.fix, sig2b.fix=sig2b.fix,
              penalty=penalty, delta=delta, shift=shift)#
}
#----------------------------------------------------------------------------------------

gamlss.mrf <- function(x, y, w, xeval = NULL, ...)
{
# -------------------------------------------------- 
# the main function starts here
# get the attributes

          xvar <-  if (is.null(xeval)) as.factor(attr(x, "x"))
                   else as.factor(attr(x, "x"))[seq(1,length(y))]
          area <- attr(x, "area")
               if (!is.null(area)) area <- as.factor(attr(x, "area"))
      precision <- as.matrix(attr(x, "precision"))  
     adj.weight <- as.numeric(attr(x, "adj.weight"))
        control <- as.list(attr(x, "control")) 
     gamlss.env <- as.environment(attr(x, "gamlss.env"))
startLambdaName <- as.character(attr(x, "NameForLambda")) 
         # order <- control$order # the order of the penalty matrix
              N <- length(y) # the no of observations
          tau2  <- sig2 <- NULL
# now the action depends on the values of lambda and df
#--------------------------------------------------------------------      
sigmasWEIGHT <- get(startLambdaName, envir=gamlss.env) ## geting the starting value
  #cat("sigmas", sigmas, "\n")
sigmas <- sigmasWEIGHT[1:2]
ifWeiAdj <- sigmasWEIGHT[3]

if (!ifWeiAdj) 
{
  fit <- MRF(y = y, x=xvar, precision=precision,  weights = w,  
             area=area, sig2e = 1, sig2b = 1, 
             sig2e.fix = control$sig2e.fix, sig2b.fix = control$sig2b.fix,  
             penalty = control$penalty,  delta = control$delta, shift = control$shift)
  if(any(abs(coef(fit))>=20))
  {
    fit <- MRF(y = y, x=xvar, precision=precision,  weights = w*adj.weight,  
               area=area, sig2e = 1, sig2b = 1, 
               sig2e.fix = control$sig2e.fix, sig2b.fix = control$sig2b.fix,  
               penalty = control$penalty,  delta = control$delta, shift = control$shift)
    ifWeiAdj <- TRUE 
  }
}
else 
{
  fit <- MRF(y = y, x=xvar, precision=precision,  weights = w*adj.weight,  
             area=area, sig2e = 1, sig2b = 1, 
             sig2e.fix = control$sig2e.fix, sig2b.fix = control$sig2b.fix,  
             penalty = control$penalty,  delta = control$delta, shift = control$shift)
}
#cat(coef(fit),"\n")
sigmas <- c(fit$sig2e, fit$sig2b,ifWeiAdj)
        assign(startLambdaName, sigmas, envir=gamlss.env)
  # var <- fit$var[fit$var!=0]
    if (is.null(xeval)) # if no prediction 
    {
     list(fitted.values=fitted(fit), residuals=y-fitted(fit), var=fit$var, nl.df =fit$df-2,
          lambda=fit$sig2e/fit$sig2b, coefSmo=fit )
    }                            
    else
    {
      nx <- as.factor(attr(x, "x"))[seq(length(y)+1, length(x) )]
    pred <- predict(fit, newdata=nx)
    pred  
    }
}
