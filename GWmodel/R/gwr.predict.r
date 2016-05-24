#########################
###Basic GWR  as a predictor
# Author: BL, PH
# dMat1: distance matrix between data points and prediction locations
# dMat2: sysmetric distance matrix between data points 
gwr.predict<-function(formula, data, predictdata, bw, kernel="bisquare",adaptive=FALSE, p=2, theta=0, longlat=F,dMat1, dMat2)
{
   ##Record the start time
  timings <- list()
  timings[["start"]] <- Sys.time()
  ###################################macth the variables
  #####Check the given data frame and prediction points
  this.call <- match.call()
  p4s <- as.character(NA)
  ##Data points for fitting GWR model
  if (is(data, "Spatial"))
  {
    p4s <- proj4string(data)
    fd.locat<-coordinates(data)
    predict.SPDF <- data
    fd.n <- nrow(fd.locat)
    data <- as(data, "data.frame")
  }
  else
  {
       stop("Given data for GWR calibration must be Spatial*DataFrame")
  }
  #########################################fit x and y
  ######Extract the data frame
  ####Refer to the function lm
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data"), names(mf), 0L)

    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    y <- model.extract(mf, "response")
    x <- model.matrix(mt, mf)
	  var.n<-ncol(x)
	  idx1 <- match("(Intercept)", colnames(x))
    if(!is.na(idx1))
      colnames(x)[idx1]<-"Intercept" 
	  inde_vars <- colnames(x)[-1]
  #######data for prediction
  #####Prediction points
  if (missing(predictdata))
  {
    warning("Data for prediction is not given and the same data as fitting data will be used!")
    predictdata <- data
    pd.locat <- fd.locat
    pd.given <- F
  }
  else
  {
    if (is(predictdata, "Spatial"))
    {
      p4s <- proj4string(predictdata)
      pd.locat<-coordinates(predictdata)
      predict.SPDF <- predictdata 
      predictdata <- as(predictdata, "data.frame")
      pd.given <- T
      if(any((inde_vars %in% names(predictdata))==F))
        stop("All the independent variables should be included in the predictdata")
    }
  }
  pd.n <- nrow(pd.locat)
  x.p <-  predictdata[,inde_vars]
  x.p <- cbind(rep(1,pd.n),x.p)
  x.p  <- as.matrix(x.p)
  #####################Distance metric
  DM1.given <- T
  DM2.given <- T
  if (missing(dMat1))
      DM1.given<-F
  if (missing(dMat2))
      DM2.given<-F
  if (!DM1.given)
  {
    if(fd.n + pd.n <= 10000)
    {
      dMat1 <- gw.dist(dp.locat=fd.locat, rp.locat=pd.locat, p=p, theta=theta, longlat=longlat)
      DM3.given <- F
      DM1.given<-T
    }
  }
  else
  {
     DM3.given <- T
     dim.dMat1<-dim(dMat1)
     if (dim.dMat1[1]!=fd.n||dim.dMat1[2]!=pd.n)
       stop("Dimensions of dMat1 (prediction matrix) are not correct")
  }
  
  if (!DM2.given)
  {
    if(fd.n <= 5000)
    {
      dMat2 <- gw.dist(dp.locat=fd.locat, rp.locat=fd.locat, p=p, theta=theta, longlat=longlat)
      DM2.given<-T
    }
  }
  else
  {
     dim.dMat2<-dim(dMat2)
     if (dim.dMat2[1]!=fd.n||dim.dMat2[2]!=fd.n)
       stop("Dimensions of dMat1 (fitting matrix) are not correct")
  }
  
  #####################GWR Fit
  #############Calibration the model
  wt <-matrix(nrow=fd.n, ncol=pd.n) 
  xtxinv <- as.double(rep(0,pd.n*var.n*var.n))
  dim(xtxinv) <- c(pd.n,var.n,var.n)
  betas1 <- matrix(nrow=pd.n, ncol=var.n)
  #predict the model
  for (i in 1:pd.n)
  {
    if (DM1.given)
       dist.vi<-dMat1[,i]
    else
    {
        dist.vi<-gw.dist(fd.locat, pd.locat, focus=i, p, theta, longlat)
    }
    W.i<-gw.weight(dist.vi,bw,kernel,adaptive)
    if (!pd.given)
      W.i[i] <- 0 
    wt[,i] <- W.i
    gw.resi<-gw.reg1(x,y,W.i,i)
    betas1[i,]<-gw.resi[[1]]
    xtxinv[i,,] <-gw.resi[[2]]
  }
  gw.predict <- gwr.fitted(x.p, betas1)
  ###### fit the model
  betas2 <- matrix(nrow=fd.n, ncol=var.n)
  S <- matrix(nrow=fd.n, ncol=fd.n)
  for (i in 1:fd.n)
  {
    if (DM2.given)
       dist.vi<-dMat2[,i]
    else
    {
        dist.vi<-gw.dist(fd.locat, focus=i, p, theta, longlat)
    }
    W.i<-gw.weight(dist.vi,bw,kernel,adaptive)
    gw.resi<-gw.reg(x,y,W.i,T,i)
    betas2[i,]<-gw.resi[[1]] ######See function by IG
    S[i,]<-gw.resi[[2]]
    #Ci<-gw.resi[[3]]
  }
  tr.S<-sum(diag(S))
  tr.StS<-sum(S^2)
  Q<-t(diag(fd.n)-S)%*%(diag(fd.n)-S)
  RSS.gw<-t(y)%*%Q%*%y
  sigma.hat<-RSS.gw/(fd.n-2*tr.S+tr.StS)
  ###Prediction variance
  predict.var <- numeric(pd.n)
  for(i in 1:pd.n)
  {
    w2 <- wt[,i]*wt[,i]
    w2x<- x * w2
    xtw2x <- t(x)%*%w2x
    # NEED TO CHOOSE ONE ACCORDING TO NUMBER OF EXPLANATORY VARIABLES (& INTERCEPT).....
    xtxinvp <- c()
    for (j in 1:var.n)
       xtxinvp <- rbind(xtxinvp, xtxinv[i,,j])
    s0 <- xtxinvp%*%xtw2x%*%xtxinvp
    x.pi <-x.p[i,]
    dim(x.pi) <- c(1, length(x.pi)) 
    s1 <- x.pi%*%s0%*%t(x.pi)
    pse <- (sqrt(sigma.hat))*(sqrt(1+s1))
    pvar <- (pse)^2
    predict.var[i] <- pvar
  }
  gwr.pred.df <- data.frame(betas1, gw.predict, predict.var)
  colnames(gwr.pred.df) <- c(paste(colnames(x), "coef", sep = "_"), "prediction", "prediction_var")
  rownames(pd.locat)<-rownames(gwr.pred.df)
  griddedObj <- F
  if (is(predict.SPDF, "Spatial"))
  { 
      if (is(predict.SPDF, "SpatialPolygonsDataFrame"))
      {
         polygons<-polygons(predict.SPDF)
         #SpatialPolygons(regression.points)
         #rownames(gwres.df) <- sapply(slot(polygons, "polygons"),
                            #  function(i) slot(i, "ID"))
         SDF <-SpatialPolygonsDataFrame(Sr=polygons, data=gwr.pred.df)
      }
      else
      {
         griddedObj <- gridded(predict.SPDF)
         SDF <- SpatialPointsDataFrame(coords=pd.locat, data=gwr.pred.df, proj4string=CRS(p4s), match.ID=F)
         gridded(SDF) <- griddedObj 
      }
  }
  else
      SDF <- SpatialPointsDataFrame(coords=pd.locat, data=gwr.pred.df, proj4string=CRS(p4s), match.ID=F)
  
 # if (is(predict.SPDF, "SpatialPolygonsDataFrame"))
#    {
#       polygons<-polygons(predict.SPDF)
#       #SpatialPolygons(regression.points)
#       rownames(gwres.df) <- sapply(slot(polygons, "polygons"),
#                            function(i) slot(i, "ID"))
#       SDF <-SpatialPolygonsDataFrame(Sr=polygons, data=gwr.pred.df, match.ID=F)
#    }
#    else
#       SDF <- SpatialPointsDataFrame(coords=pd.locat, data=gwr.pred.df, proj4string=CRS(p4s))
       
    timings[["stop"]] <- Sys.time()
  ##############
  GW.arguments<-list(formula=formula,bw=bw, kernel=kernel,adaptive=adaptive, p=p, theta=theta, longlat=longlat, fd.n=fd.n, DM.given=DM3.given)
  res<-list(GW.arguments=GW.arguments,SDF=SDF,timings=timings,this.call=this.call)
  class(res) <-"gwrm.pred"
  invisible(res)
}

gw.reg1<-function(X,Y,W.i,focus)
###GWR for prediction
{
    xtxinv <- solve(t(X*W.i)%*%X)
    xty <- t(X*W.i)
    betai<-xtxinv%*%xty%*%Y
    res<-list(betai,xtxinv,xty) 
    res    
}

############################Layout function for outputing the GWR results
##Author: BL	
print.gwrm.pred<-function(x, ...)
{
  if(class(x) != "gwrm.pred") stop("It's not a gwm object")
  cat("   ***********************************************************************\n")
  cat("   *                       Package   GWmodel                             *\n")
  cat("   ***********************************************************************\n")
  cat("   Program starts at:", as.character(x$timings$start), "\n")
  cat("   Call:\n")
  cat("   ")
  print(x$this.call)
  vars<-all.vars(x$GW.arguments$formula)
  var.n<-length(vars)
	cat("\n   Dependent (y) variable for prediction: ",vars[1])
	cat("\n   Independent variables: ",vars[-1])
	fd.n<-x$GW.arguments$fd.n
	cat("\n   Number of data points:",fd.n)	
	#########################################################################
	cat("\n   ***********************************************************************\n")
	  cat("   *     Results of Geographically Weighted Regression for prediction    *\n")
	cat("   ***********************************************************************\n")
	cat("\n   *********************Model calibration information*********************\n")
	cat("   Kernel function:", x$GW.arguments$kernel, "\n")
	if(x$GW.arguments$adaptive)
	   cat("   Adaptive bandwidth: ", x$GW.arguments$bw, " (number of nearest neighbours)\n", sep="") 
  else
     cat("   Fixed bandwidth:", x$GW.arguments$bw, "\n")
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
	SDF.df <-as(x$SDF, "data.frame")
	cat("\n   ****************Summary of GWR coefficient estimates:******************\n")       
		df0 <- SDF.df[,1:var.n, drop=FALSE]
        if (any(is.na(df0))) {
            df0 <- na.omit(df0)
            warning("NAs in coefficients dropped")
        }
	CM <- t(apply(df0, 2, summary))[,c(1:3,5,6)]
	rnames<-rownames(CM)
		for (i in 1:length(rnames))
			 rnames[i]<-paste("   ",rnames[i],sep="")
	rownames(CM) <-rnames 
	printCoefmat(CM)
	cat("\n   ****************       Results of GW prediction       ******************\n")
  #m <- ncol(SDF.df)       
		df1 <- SDF.df[,c("prediction", "prediction_var"), drop=FALSE]
        if (any(is.na(df0))) {
            df0 <- na.omit(df0)
            warning("NAs in coefficients dropped")
        }
	CM1 <- t(apply(df1, 2, summary))[,c(1:3,5,6)]
	rnames<-rownames(CM1)
		for (i in 1:length(rnames))
			 rnames[i]<-paste("   ",rnames[i],sep="")
	rownames(CM1) <-rnames 
	printCoefmat(CM1)
	cat("\n   ***********************************************************************\n")
	cat("   Program stops at:", as.character(x$timings$stop), "\n")
	invisible(x)
}