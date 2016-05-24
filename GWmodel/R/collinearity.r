#
#
#########################################################################
#########################################################################
# gw_ridge_regression
# Author: MC, PH, BL


#######################################################################
#######################################################################
# ########## ridge.gwr ########## main ridge routine
#
# y               response variable
# x               predictor variables
# locs            X,Y for locations
# weight          =g: fixed radius Gaussian kernel
#                 =b: variable radius Bisquare kernel
#                 =i: fixed radius bisquare kernel
#                 =x: fixed radius boxcar kernel
# bandwidth       =g|i|x: in locs units
#                 =b: in nn
#
# lambda           ridge parameter for gwr ridge model
# cv              =T: return cv score
# [returns]       cv score only
#
# lambda.adjust    =F standard gwr (i.e. ridge=0) with local ridge adjustment
#                 =T use cn.thresh as maximum CN
# cn.thresh       desired maximum condition number for local adjustment
# [returns]       parameter esitmates, standard errors, goodness of fit, CN, ridge
#
# non.sample      =F standard gwr location fitting
#                 =T use fitlocs as fitting locations
# fitlocs         X,Y for fitting locations
# [returns]       parameter estimates, standard errors
#
gwr.lcr <-function(formula, data, regression.points, bw, kernel="bisquare",
                    lambda=0,lambda.adjust=FALSE,cn.thresh=NA,
                    adaptive=FALSE, p=2, theta=0, longlat=F,cv=T,dMat)
{
  ##Record the start time
  timings <- list()
  timings[["start"]] <- Sys.time()
  ###################################macth the variables
  this.call <- match.call()
  p4s <- as.character(NA)
  #####Check the given data frame and regression points
  #####Regression points
  if (missing(regression.points))
  {
  	rp.given <- FALSE
    regression.points <- data
    rp.locat <- coordinates(data)
    hatmatrix<-T
  }
  else
  {
    rp.given <- TRUE
    hatmatrix<-F
    if (is(regression.points, "Spatial"))
    {
       rp.locat<-coordinates(regression.points)
    }
    else if (is.numeric(regression.points) && dim(regression.points)[2] == 2)
       rp.locat<-regression.points
    else
      {
        warning("Output loactions are not packed in a Spatial object,and it has to be a two-column numeric vector")
        rp.locat<-dp.locat
      }
  }
  ##Data points{
  if (is(data, "Spatial"))
  {
    p4s <- proj4string(data)
    dp.locat<-coordinates(data)
    data <- as(data, "data.frame")
  }
  else
  {
    stop("Given regression data must be Spatial*DataFrame")
  }

  ####################
  ######Extract the data frame
  ####Refer to the function lm
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data"), names(mf), 0L)

    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    y <- model.extract(mf, "response")
    x <- model.matrix(mt, mf)
    idx1 <- match("(Intercept)", colnames(x))
    if(!is.na(idx1))
      colnames(x)[idx1]<-"Intercept" 
	  var.n<-ncol(x)
	  rp.n<-nrow(rp.locat)
    dp.n<-nrow(data)
    if (missing(dMat))
    {
      DM.given<-F
      DM1.given<-F
      if(dp.n + rp.n <= 10000)
      {
        dMat <- gw.dist(dp.locat=dp.locat, rp.locat=rp.locat, p=p, theta=theta, longlat=longlat)
        DM.given<-T
      }
    }
    else
    {
      DM.given<-T
      DM1.given<-T
      dim.dMat<-dim(dMat)
      if (dim.dMat[1]!=dp.n||dim.dMat[2]!=rp.n)
         stop("Dimensions of dMat are not correct")
    }

   # arrays for the results
   betas <-matrix(nrow=rp.n, ncol=var.n)
   local.cn <- rep(0,rp.n)
   local.lambda <- rep(0,rp.n)
   hatrow <- rep(0,rp.n)
   trs <- 0
   trsts <- 0
   colnames(betas) <- colnames(x)
   #colnames(betas)[1]<-"Intercept"

   #######################
   ###### Main Loop ######
   #######################
   for (i in 1:rp.n)
   {
      if (DM.given)
         dist.vi<-dMat[,i]
      else
      {
        if (rp.given)
          dist.vi<-gw.dist(dp.locat, rp.locat, focus=i, p, theta, longlat)
        else
          dist.vi<-gw.dist(dp.locat=dp.locat, focus=i, p=p, theta=theta, longlat=longlat)
      }
      W.i<-gw.weight(dist.vi,bw,kernel,adaptive)
      # CV=TRUE: apply crossvalidation the weight for the focus point is zero
      #
      #if (cv) {W.i[i] <- 0}
      #yw <- y * wgt
      xw <- as.matrix(sweep(x[,-1],1,W.i,"*"))
      #
      # Condition number check
      # Intercept
      x1w <- as.matrix(sweep(x,1,W.i,"*"))                    # Weight it
      svd.x <- svd(sweep(x1w, 2,sqrt(colSums(x1w^2)), "/"))$d  # SVD
      local.cn[i] <- svd.x[1] / svd.x[var.n]

      # this is the currently set global value of lambda
      local.lambda[i] <- lambda

       # Lambda adjustment for the locally compensated model
       if (lambda.adjust)
       {
          if (local.cn[i] > cn.thresh)
          {
             local.lambda[i] <- (svd.x[1] - cn.thresh*svd.x[var.n]) / (cn.thresh - 1)
          }
       }
      # Now, fit the ridge regression, and save the coefficients for this run
      ridge.fit <- ridge.lm(y,x,W.i,local.lambda[i],add.int=F)      # Compute the coefficients
      betas[i,] <- ridge.fit$coef                     # store adjusted coefficients for focus
      if (!rp.given)
      {
         xm <- x                                 # unweighted design matrix
         xtwx <- t(xm) %*% diag(W.i) %*% xm                  # x'wx
         xtwxinv <- solve(xtwx)   # (x'wx)^-1
         hatrow <- x1w[i,] %*% xtwxinv %*% t(x1w)        # focus'th row of the hat matrix
         trs <- trs + hatrow[i]                          # running sum for tr(S)
         trsts <- trsts + sum(hatrow^2)                      # running sum for tr(S'S)
      }
   }
   #################################
   ####### End of main loop ########
   #################################

   if (!rp.given)
   {
      # Fitted values, residuals, and rss
      yhat <- rowSums(x * betas)   # fitted values
      residual <- y - yhat                     # residuals
      rss <- sum(residual^2)                   # residual sum of squares

      enp <- 2*trs - trsts        # effective number of parameters
      edf <- dp.n - enp              # effective degrees of freedom of the residual

      s2 <- rss/(dp.n - enp)         # sigma-squared

      # Hurvich CM, Tsai C-L, 1991, Biometrike, 78(3), 499-509
      aic <-  dp.n*(log(2*pi*s2) + 1) + 2*(enp + 1)
      aicc <- dp.n*log(2*pi*s2) + dp.n*( (1+enp/dp.n) / (1-(enp+2)/dp.n) )
      CV<-numeric(dp.n)
      if(cv)
        CV<-gwr.lcr.cv.contrib(bw,x,y,dp.locat,kernel,lambda,lambda.adjust,cn.thresh,adaptive, p, theta, longlat,dMat)
      GW.diagnostic<-list(AIC=aic,AICc=aicc,enp=enp, edf=edf,RSS=rss)
   }
   ####encapsulate the GWR results
   GW.arguments<-list(formula=formula,rp.given=rp.given,bw=bw,
                       kernel=kernel,adaptive=adaptive, p=p, theta=theta, longlat=longlat,DM.given=DM1.given,
                       lambda=lambda,lambda.adjust=lambda.adjust,cn.thresh=cn.thresh)
   if (!rp.given)
    {
      gwres.df<-data.frame(betas,y,yhat,residual,local.cn, local.lambda,CV)
      colnames(gwres.df)<-c(colnames(betas),c("y","yhat","residual","Local_CN","Local_Lambda","CV_Score"))
    }
    else
    {
      gwres.df<-data.frame(betas,local.cn, local.lambda)
      colnames(gwres.df)<-c(colnames(betas),c("Local_CN","Local_Lambda"))
    }
    rownames(rp.locat)<-rownames(gwres.df)
    griddedObj <- F
    if (is(regression.points, "Spatial"))
    { 
        if (is(regression.points, "SpatialPolygonsDataFrame"))
        {
           polygons<-polygons(regression.points)
           #SpatialPolygons(regression.points)
           #rownames(gwres.df) <- sapply(slot(polygons, "polygons"),
                              #  function(i) slot(i, "ID"))
           SDF <-SpatialPolygonsDataFrame(Sr=polygons, data=gwres.df)
        }
        else
        {
           griddedObj <- gridded(regression.points)
           SDF <- SpatialPointsDataFrame(coords=rp.locat, data=gwres.df, proj4string=CRS(p4s))
           gridded(SDF) <- griddedObj 
        }
    }
    else
        SDF <- SpatialPointsDataFrame(coords=rp.locat, data=gwres.df, proj4string=CRS(p4s))
    timings[["stop"]] <- Sys.time()

     if(rp.given)
     {
       res <- list(SDF=SDF,GW.arguments=GW.arguments,timings=timings,this.call=this.call)
     }
     else
       res <- list(SDF=SDF,GW.arguments=GW.arguments, GW.diagnostic=GW.diagnostic,timings=timings,this.call=this.call)
     class(res) <-"gwrlcr"
     invisible(res)
}

############################Layout function for outputing the GWR results
##Author: BL
print.gwrlcr<-function(x, ...)
{
  if(class(x) != "gwrlcr") stop("It's not a gwrlcr object")
  cat("   ***********************************************************************\n")
  cat("   *                       Package   GWmodel                             *\n")
  cat("   ***********************************************************************\n")
  cat("   Program starts at:", as.character(x$timings$start), "\n")
  cat("   Call:\n")
  cat("   ")
  print(x$this.call)
  vars<-all.vars(x$GW.arguments$formula)
  var.n<-length(vars)
	cat("\n   Dependent (y) variable: ",vars[1])
	cat("\n   Independent variables: ",vars[-1])
	cat("\n   ***********************************************************************\n")
	  cat("   *        Results of Ridge Geographically Weighted Regression           *\n")
	cat("   ***********************************************************************\n")
	cat("\n   *********************Model calibration information*********************\n")
	cat("   Kernel function:", x$GW.arguments$kernel, "\n")
	if(x$GW.arguments$adaptive)
	   cat("   Adaptive bandwidth: ", x$GW.arguments$bw, " (number of nearest neighbours)\n", sep="")
  else
     cat("   Fixed bandwidth:", x$GW.arguments$bw, "\n")
	if(x$GW.arguments$rp.given)
     cat("   Regression points: A seperate set of regression points is used.\n")
  else
     cat("   Regression points: the same locations as observations are used.\n")
	if (x$GW.arguments$DM.given)
     cat("   Distance metric: A distance matrix is specified for this model calibration.\n")
  else
   {
     if (x$GW.arguments$longlat)
        cat("   Distance metric: Great Circle distance metric is used.\n")
     else if (x$GW.arguments$p==2)
        cat("   Distance metric: Euclidean distance metric is used.\n")
     else if (x$GW.arguments$p==1)
        cat("   Distance metric: Manhattan distance metric is used.\n")
     else if (is.infinite(x$GW.arguments$p))
        cat("   Distance metric: Chebyshev distance metric is used.\n")
     else
        cat("   Distance metric: A generalized Minkowski distance metric is used with p=",x$GW.arguments$p,".\n")
     if (x$GW.arguments$theta!=0&&x$GW.arguments$p!=2&&!x$GW.arguments$longlat)
        cat("   Coordinate rotation: The coordinate system is rotated by an angle", x$GW.arguments$theta, "in radian.\n")
   }
   cat("   Lambda(ridge parameter for gwr ridge model):", x$GW.arguments$lambda)
   if(x$GW.arguments$lambda.adjust)
   {
     cat("   The threshold of condition number is used for adjusting local Lambda:", x$GW.arguments$cn.thresh)
   }

	cat("\n   **************Summary of Ridge GWR coefficient estimates:***************\n")
		df0 <- as(x$SDF, "data.frame")[,1:var.n, drop=FALSE]
        if (any(is.na(df0))) {
            df0 <- na.omit(df0)
            warning("NAs in coefficients dropped")
        }
    dp.n <- nrow(df0)
	CM <- t(apply(df0, 2, summary))[,c(1:3,5,6)]
	rnames<-rownames(CM)
		for (i in 1:length(rnames))
			 rnames[i]<-paste("   ",rnames[i],sep="")
	rownames(CM) <-rnames
	printCoefmat(CM)
	cat("   ************************Diagnostic information*************************\n")
	if (!x$GW.arguments$rp.given)
  {
		cat("   Number of data points:", dp.n, "\n")
		cat("   Effective number of parameters (2trace(S) - trace(S'S)):", x$GW.diagnostic$enp, "\n")
		cat("   Effective degrees of freedom (n-2trace(S) + trace(S'S)):", x$GW.diagnostic$edf, "\n")
		cat("   AICc (GWR book, Fotheringham, et al. 2002, p. 61, eq 2.33):",
                    x$GW.diagnostic$AICc, "\n")
		cat("   AIC (GWR book, Fotheringham, et al. 2002,GWR p. 96, eq. 4.22):", x$GW.diagnostic$AIC, "\n")
		cat("   Residual sum of squares:", x$GW.diagnostic$RSS, "\n")
  }
	cat("\n   ***********************************************************************\n")
	cat("   Program stops at:", as.character(x$timings$stop), "\n")
	invisible(x)
}
#######################################################################
############### ridge.lm ######## fitting routine #####################
#######################################################################
ridge.lm <- function(y,X,w=rep(1,length(y)),lambda=0,add.int=TRUE) {
	if (add.int) X <- cbind(1,X)  # Add a column of 1s if required
	Xw <- sqrt(w) * X             # Multiply by sqrt weights
	yw <- sqrt(w) * y             # Multiply by sqrt weights
  Xsd <- c(1,apply(X[,-1],2,sd))
  #if (!add.int) # Compute sd's for X - set the scale factor to 1 for intercepts,  if they are there...
	#	{Xsd <- c(1,apply(X[,-1],2,sd))}
	#else
	#	{Xsd <- apply(X,2,sd)}
	Xws <- as.matrix(sweep(Xw,2,Xsd,"/")) # Divide the X columns by the sd's
    ysd <- sd(y)                              # Compute sd for y
    yws <- yw / ysd                           # divide y by its sd
    # Final line does the ridge regression and back-transforms the coefs
    b <- t(solve(crossprod(Xws)+lambda*diag(ncol(Xws)),crossprod(Xws,yws))*ysd / Xsd)

    if (add.int) colnames(b)[1] <- "Intercept"
    return(list(coef=b))
    }
#######################################################################
############### ridge gwr bandwidth #####################
#######################################################################
bw.gwr.lcr <-function(formula, data, kernel="bisquare",
                    lambda=0,lambda.adjust=FALSE,cn.thresh=NA,
                    adaptive=FALSE, p=2, theta=0, longlat=F,dMat)
{
  ##Record the start time
  timings <- list()
  timings[["start"]] <- Sys.time()
  ###################################macth the variables
  this.call <- match.call()
  p4s <- as.character(NA)
  #####Check the given data frame and regression points
  #####Regression points
  ##Data points{
  if (is(data, "Spatial"))
  {
    p4s <- proj4string(data)
    dp.locat<-coordinates(data)
    rp.locat<-dp.locat
    data <- as(data, "data.frame")
  }
  else
  {
    stop("Given regression data must be Spatial*DataFrame")
  }

  ####################
  ######Extract the data frame
  ####Refer to the function lm
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data"), names(mf), 0L)

    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    y <- model.extract(mf, "response")
    x <- model.matrix(mt, mf)
	  var.n<-ncol(x)
	  rp.n<-nrow(rp.locat)
    dp.n<-nrow(data)
    if (missing(dMat))
    {
      DM.given<-F 
      if(dp.n + rp.n <= 10000)
      {
        dMat <- gw.dist(dp.locat=dp.locat, rp.locat=rp.locat, p=p, theta=theta, longlat=longlat)
        DM.given<-T
      }
    }
    else
    {
      DM.given<-T
      dim.dMat<-dim(dMat)
      if (dim.dMat[1]!=dp.n||dim.dMat[2]!=rp.n)
         stop("Dimensions of dMat are not correct")
    }
   #########Find the range for the bandwidth selection
  if(adaptive)
  {
    upper<-dp.n
    lower<-20
  }
  else
  {
    if(DM.given)
    {
      upper<-range(dMat)[2]
      lower<-upper/5000
    }
    ###!!!!!!!! Important note: if the distance matrix is not specified, the running time will be consuming very much by choosing the range of fixed bandwidth when the p is not 2;
    ### because the range can't be decided by boundary box
    else
    {
      dMat<-NULL
      if (p==2)
      {
        b.box<-bbox(dp.locat)
        upper<-sqrt((b.box[1,2]-b.box[1,1])^2+(b.box[2,2]-b.box[2,1])^2)
        lower<-upper/5000
      }
      else
      {
        upper<-0
        for (i in 1:dp.n)
        {
          dist.vi<-gw.dist(dp.locat=dp.locat, focus=i, p=p, theta=theta, longlat=longlat)
          upper<-max(upper, range(dist.vi)[2])
        }
        lower<-upper/5000
      }
    }
  }
  #Select the bandwidth by golden selection
  bw<-NA
  bw <- gold(gwr.lcr.cv,lower,upper,adapt.bw=adaptive,x,y,dp.locat,kernel,lambda,
             lambda.adjust,cn.thresh,adaptive, p, theta, longlat,dMat)
  bw
}

gwr.lcr.cv <- function(bw,X,Y,locs,kernel="bisquare",
                    lambda=0,lambda.adjust=FALSE,cn.thresh=NA,
                    adaptive=FALSE, p=2, theta=0, longlat=F,dMat)
{
   # get the dimensions of the problem
   n <- dim(X)[1]
   m <- dim(X)[2]
   fitlocs <- locs
   # arrays for the results
   betas <- array(0, dim=c(n,m))
   local.cn <- rep(0,n)
   local.lambda <- rep(0,n)
   if (missing(dMat))
    DM.given<-F
   else
    DM.given<-T
   #######################
   ###### Main Loop ######
   #######################
   for (focus in 1:n)
   {

      if (DM.given)
        dist.vi<-dMat[,focus]
      else
      {
        dist.vi<-gw.dist(locs, focus=focus, p=p, theta=theta, longlat=longlat)
      }
	  	wgt <-gw.weight(dist.vi,bw,kernel,adaptive)
      wgt[focus] <- 0
      #yw <- y * wgt
      xw <- as.matrix(sweep(X[,-1],1,wgt,"*"))
      #
      # Condition number check
      # Intercept
      x1w <- as.matrix(sweep(X,1,wgt,"*"))                    # Weight it
      svd.x <- svd(sweep(x1w, 2,sqrt(colSums(x1w^2)), "/"))$d  # SVD
      local.cn[focus] <- svd.x[1] / svd.x[m]

      # this is the currently set global value of lambda
      local.lambda[focus] <- lambda
         # Lambda adjustment for the locally compensated model
         if (lambda.adjust) {
            if (local.cn[focus] > cn.thresh) {
               local.lambda[focus] <- (svd.x[1] - cn.thresh*svd.x[m]) / (cn.thresh - 1)
               #print(paste("Adjust",lambda,lambda.local,local.cn[focus]))
            }
         }

      # Now, fit the ridge regression, and save the coefficients for this run
      ridge.fit <- ridge.lm(Y,X,wgt,local.lambda[focus],add.int=F)      # Compute the coefficients
      betas[focus,] <- ridge.fit$coef                     # store adjusted coefficients for focus
   }
   yhat <- rowSums(X * betas)
   residual <- Y - yhat
   cv <- sum(residual^2)
   if(adaptive)
    cat("Adaptive bandwidth(number of nearest neighbours):", bw, "CV score:", cv, "\n")
  else
    cat("Fixed bandwidth:", bw, "CV score:", cv, "\n")
   cv
}
 #Contribution of each observation to the score statistic used in cross-validation for lcr.gwr
gwr.lcr.cv.contrib <- function(bw,X,Y,locs,kernel="bisquare",
                    lambda=0,lambda.adjust=FALSE,cn.thresh=NA,
                    adaptive=FALSE, p=2, theta=0, longlat=F,dMat)
{
   # get the dimensions of the problem
   n <- dim(X)[1]
   m <- dim(X)[2]
   fitlocs <- locs
   # arrays for the results
   betas <- array(0, dim=c(n,m))
   local.cn <- rep(0,n)
   local.lambda <- rep(0,n)
   if (missing(dMat))
    DM.given<-F
   else
    DM.given<-T
   #######################
   ###### Main Loop ######
   #######################
   for (focus in 1:n)
   {

      if (DM.given)
        dist.vi<-dMat[,focus]
      else
      {
        dist.vi<-gw.dist(locs, focus=focus, p=p, theta=theta, longlat=longlat)
      }
	  	wgt <-gw.weight(dist.vi,bw,kernel,adaptive)
      wgt[focus] <- 0
      #yw <- y * wgt
      xw <- as.matrix(sweep(X[,-1],1,wgt,"*"))
      #
      # Condition number check
      # Intercept
      x1w <- as.matrix(sweep(X,1,wgt,"*"))                    # Weight it
      svd.x <- svd(sweep(x1w, 2,sqrt(colSums(x1w^2)), "/"))$d  # SVD
      local.cn[focus] <- svd.x[1] / svd.x[m]

      # this is the currently set global value of lambda
      local.lambda[focus] <- lambda
         # Lambda adjustment for the locally compensated model
         if (lambda.adjust) {
            if (local.cn[focus] > cn.thresh) {
               local.lambda[focus] <- (svd.x[1] - cn.thresh*svd.x[m]) / (cn.thresh - 1)
               #print(paste("Adjust",lambda,lambda.local,local.cn[focus]))
            }
         }

      # Now, fit the ridge regression, and save the coefficients for this run
      ridge.fit <- ridge.lm(Y,X,wgt,local.lambda[focus],add.int=F)      # Compute the coefficients
      betas[focus,] <- ridge.fit$coef                     # store adjusted coefficients for focus
   }
   yhat <- rowSums(X * betas)
   cv <- Y - yhat
   cv
}
