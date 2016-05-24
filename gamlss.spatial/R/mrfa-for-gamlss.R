#----------------------------------------------------------------------------------
# this is the mrf method using alternating estimation for sigmas
# FDB - 
#-------------------------------------------------------------------------------
mrfa<-function(x, precision=NULL, neighbour=NULL, polys= NULL, area=NULL, start=10, df=NULL, adj.weight=1000, ...) 
{
  #require(spam)
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
  ## this code is similar to the starting code in MRFA()
  ## and it takes a polys or neighbour informationand creates 
  ## the precision matrix
  ##---------------------
  if (!is(x, "factor")) stop("x must be a factor")
       k <- area
  if (is.null(k)) k <- factor(levels(x),levels=levels(x)) # 
  #the line above can change the order of the precision matrix
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
    neighbour <- polys2nb(polys)      
    precision   <- nb2prec(neighbour=neighbour,x=x,area=area) 
  }
  ## finish precision ----------
  ##
  ## get a random name to use it in the gamlss() environment
  #--------
               sl <- sample(letters, 4)
      fourLetters <- paste(paste(paste(sl[1], sl[2], sep=""), sl[3], sep=""),sl[4], sep="")
  startLambdaName <- paste("start.Lambda",fourLetters, sep=".")
  ## put the starting values in the gamlss()environment
    weights.Adj  <- FALSE
  lambda <- c(start, weights.Adj)
  assign(startLambdaName, lambda, envir = gamlss.env)
  #--------
  xvar <- rep(0,length(x)) #  set x to zero
  #  attr(xvar, "control")       <- control
  attr(xvar, "call")          <- substitute(gamlss.mrfa(data[[scall]], z, w)) 
  attr(xvar, "gamlss.env")    <- gamlss.env
  attr(xvar, "x")             <- x
  attr(xvar, "df")            <- df
  attr(xvar, "area")          <- area
  attr(xvar, "adj.weight")    <- adj.weight
  attr(xvar, "precision")     <- precision
  attr(xvar, "NameForLambda") <- startLambdaName
  attr(xvar, "class")         <- "smooth"
  xvar
}
#----------------------------------------------------------------------------------------
gamlss.mrfa <- function(x, y, w, xeval = NULL, ...)
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
  gamlss.env <- as.environment(attr(x, "gamlss.env"))
          df <- attr(x, "df")
  startLambdaName <- as.character(attr(x, "NameForLambda")) 
           N <- length(y) # the no of observations
  #     tau2  <- sig2 <- NULL
  # now the action depends on the values of lambda and df
  #--------------------------------------------------------------------  
  lambdaWEIGHT <- get(startLambdaName, envir=gamlss.env) ## geting the starting val
        lambda <- lambdaWEIGHT[1]
      ifWeiAdj <- lambdaWEIGHT[2]
      if (!ifWeiAdj) 
       {
        fit <- MRFA( y, xvar, precision=precision, area=area, weights=w, start=lambda)
        if (any(coef(fit)<1e-8))
         {
          fit <- MRFA( y, xvar, precision=precision, area=area, weights=w*adj.weight, start=lambda)
          ifWeiAdj <- TRUE 
         }
       }
      else 
       {
  fit <- MRFA( y, xvar, precision=precision, area=area, weights=w*adj.weight, start=lambda)
       }
lambda <- c(fit$lambda, ifWeiAdj) 
assign(startLambdaName, lambda, envir=gamlss.env)
 # var <- fit$var[fit$var!=0]
if (is.null(xeval)) # if no prediction 
{
  list(fitted.values=fitted(fit), residuals=y-fitted(fit), var=fit$var, nl.df =fit$df-2, lambda=fit$lambda, coefSmo=fit )
}                            
else
{
  lx <- length(attr(x,"x"))
  nx <- as.factor(attr(x, "x"))[seq(length(y)+1, lx )]
  pred <- predict(fit, newdata=nx)
  pred  
}
}
#-------------------------------------------------------------------------------------------
#
