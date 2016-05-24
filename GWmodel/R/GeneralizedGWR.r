##Generalized GWR functions
#Author: Binbin Lu
#Referenced to Chris' R code
gwr.fitted <- function(x,b) apply(x*b,1,sum)

gwr.generalised<-function(formula, data, regression.points, bw, family ="poisson", kernel="bisquare",
              adaptive=FALSE, p=2, theta=0, longlat=F, dMat, cv=T,tol=1.0e-5, maxiter=20)
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
    hatmatrix<-T
  }
  else 
  {
    rp.given <- TRUE
    hatmatrix<-F
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
    ############################################
    ##Generalized linear regression
    glm.res<-glm.fit(x, y, family = poisson())  
    ############################################
	  var.n<-ncol(x)
	  rp.locat<-coordinates(regression.points)
	  rp.n<-nrow(rp.locat)
    dp.n<-nrow(data)
    betas <-matrix(nrow=rp.n, ncol=var.n)
    betas1<- betas
    betas.SE <-matrix(nrow=rp.n, ncol=var.n)
    betas.TV <-matrix(nrow=rp.n, ncol=var.n)
    ##S: hatmatrix
    S<-matrix(nrow=dp.n,ncol=dp.n)
    #C.M<-matrix(nrow=dp.n,ncol=dp.n)
    idx1 <- match("(Intercept)", colnames(x))
    if(!is.na(idx1))
      colnames(x)[idx1]<-"Intercept" 
    colnames(betas) <- colnames(x)
    #colnames(betas)[1]<-"Intercept"
    ####################################################GWR
	  #########Distance matrix is given or not

    if (missing(dMat))
    {
      DM1.given<-F
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
      DM1.given<-T
      dim.dMat<-dim(dMat)
      if (dim.dMat[1]!=dp.n||dim.dMat[2]!=rp.n)
         stop("Dimensions of dMat are not correct")
    }
    ####Generate the weighting matrix
     #############Calibration the model
     W1.mat<-matrix(numeric(dp.n*dp.n),ncol=dp.n)
     W2.mat<-matrix(numeric(dp.n*rp.n),ncol=rp.n)
    for (i in 1:dp.n)
    {
      if (DM.given)
         dist.vi<-dMat[,i]
      else
      {
          dist.vi<-gw.dist(dp.locat=dp.locat, focus=i, p=p, theta=theta, longlat=longlat) 
      }
      W.i<-gw.weight(dist.vi,bw,kernel,adaptive)
      W1.mat[,i]<-W.i
    }
    if (rp.given)
    {
      for (i in 1:dp.n)
      {
        if (DM.given)
           dist.vi<-dMat[,i]
        else
        {
            dist.vi<-gw.dist(dp.locat, rp.locat, focus=i, p, theta, longlat) 
        }
        W.i<-gw.weight(dist.vi,bw,kernel,adaptive)
        W2.mat[,i]<-W.i
      }
    }
    else
        W2.mat<-W1.mat
    
    ##model calibration
    if(family=="poisson")
      res1<-gwr.poisson(y,x,regression.points,W1.mat,W2.mat,hatmatrix,tol, maxiter)
    if(family=="binomial")  
      res1<-gwr.binomial(y,x,regression.points,W1.mat,W2.mat,hatmatrix,tol, maxiter)
    ####################################
    CV <- numeric(dp.n)
    if(hatmatrix && cv)
    {
       CV <- ggwr.cv.contrib(bw, x, y,family, kernel,adaptive, dp.locat, p, theta, longlat,dMat)
    }
      ####encapsulate the GWR results
      GW.arguments<-list()
     GW.arguments<-list(formula=formula,rp.given=rp.given,hatmatrix=hatmatrix,bw=bw, family=family,
                       kernel=kernel,adaptive=adaptive, p=p, theta=theta, longlat=longlat,DM.given=DM1.given)    
    
    timings[["stop"]] <- Sys.time()
   ##############
    res<-list(GW.arguments=GW.arguments,GW.diagnostic=res1$GW.diagnostic,glm.res=res1$glm.res,SDF=res1$SDF,CV=CV,timings=timings,this.call=this.call)
    class(res) <-"ggwrm"
    invisible(res) 
}



############ Possipon GWGLM
gwr.poisson<-function(y,x,regression.points,W1.mat,W2.mat,hatmatrix,tol=1.0e-5, maxiter=20)
{
    p4s <- as.character(NA)
    if (is(regression.points, "Spatial"))
    {
      p4s <- proj4string(regression.points)
    }
    ############################################
    ##Generalized linear regression
    glm.res<-glm.fit(x, y, family = poisson())  
    ############################################
	  var.n<-ncol(x)
	  rp.locat<-coordinates(regression.points)
	  rp.n<-nrow(rp.locat)
    dp.n<-nrow(x)
    betas <-matrix(nrow=rp.n, ncol=var.n)
    betas1<- betas
    betas.SE <-matrix(nrow=rp.n, ncol=var.n)
    betas.TV <-matrix(nrow=rp.n, ncol=var.n)
    ##S: hatmatrix
    S<-matrix(nrow=dp.n,ncol=dp.n)
    #C.M<-matrix(nrow=dp.n,ncol=dp.n)
    colnames(betas) <- colnames(x)
   # colnames(betas)[1]<-"Intercept" 
    ####################################
    ##model calibration
    it.count <- 0
    llik <- 0.0
    mu <- y + 0.1
    nu <- log(mu)
    cat(" Iteration    Log-Likelihood\n=========================\n")
    wt2 <- rep(1,dp.n)
    repeat {
     y.adj <- nu + (y - mu)/mu
     for (i in 1:dp.n)
     {
        W.i<-W1.mat[,i]
        gw.resi<-gw.reg(x,y.adj,W.i*wt2,hatmatrix=F,i)
        betas1[i,]<-gw.resi[[1]]
     }
     nu <- gwr.fitted(x,betas1)
     mu <- exp(nu)
     old.llik <- llik
     llik <- sum(y*nu - mu - log(gamma(y+1)))
     cat(paste("   ",formatC(it.count,digits=4,width=4),"    ",formatC(llik,digits=4,width=7),"\n"))
     if (abs((old.llik - llik)/llik) < tol) break
     wt2 <- mu
     it.count <- it.count+1
     if (it.count == maxiter) break}
     if(hatmatrix)
     {
        for (i in 1:rp.n)
        { 
          W.i<-W2.mat[,i]
          gw.resi<-gw.reg(x,y.adj,W.i*wt2,hatmatrix,i)
          betas[i,]<-gw.resi[[1]] ######See function by IG
          S[i,]<-gw.resi[[2]]
          Ci<-gw.resi[[3]]
          #betas.SE[i,]<-diag(Ci%*%t(Ci)) 
          invwt2 <- 1.0 /wt2
          betas.SE[i,] <- diag((Ci*invwt2) %*% t(Ci))# diag(Ci/wt2%*%t(Ci))  #see Nakaya et al. (2005) 
        }
        tr.S<-sum(diag(S))
        tr.StS<-sum(S^2)
        ###edf is different from the definition in Chris' code
        edf<-dp.n-2*tr.S+tr.StS
        yhat<-gwr.fitted(x, betas)
        residual<-y-exp(yhat)
        ########rss <- sum((y - gwr.fitted(x,b))^2)
        rss <- sum((y-exp(yhat))^2)
        sigma.hat <- rss/edf
        sigma.aic <- rss/dp.n
        for(i in 1:dp.n)
        {
           betas.SE[i,]<-sqrt(sigma.hat*betas.SE[i,])
           betas.TV[i,]<-betas[i,]/betas.SE[i,] 
        }
        #AICc <- -2*llik + 2*tr.S*dp.n/(dp.n-tr.S-2)
        AICc <- -2*llik + 2*tr.S + 2*tr.S*(tr.S+1)/(dp.n-tr.S-1)  # This is generic form of AICc (TN)
        yss.g <- sum((y - mean(y))^2)
        gw.R2<-1-rss/yss.g; ##R Square valeu
        gwR2.adj<-1-(1-gw.R2)*(dp.n-1)/(edf-1) #Adjusted R squared value
        GW.diagnostic<-list(rss=rss,AICc=AICc,edf=edf,gw.R2=gw.R2,gwR2.adj=gwR2.adj)        
     }
     else
     {
        for (i in 1:rp.n)
        { 
          W.i<-W2.mat[,i]
          gw.resi<-gw.reg(x,y.adj,W.i*wt2,hatmatrix,i)
          betas[i,]<-gw.resi[[1]] ######See function by IG
        }
     }
    if (hatmatrix)                                         
    {
      gwres.df<-data.frame(betas,y,exp(yhat),residual,betas.SE,betas.TV)
      colnames(gwres.df)<-c(c(c(colnames(betas),c("y","yhat","residual")),paste(colnames(betas), "SE", sep="_")),paste(colnames(betas), "TV", sep="_"))
    }
    else
    {
      gwres.df<-data.frame(betas)
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
           SDF <- SpatialPointsDataFrame(coords=rp.locat, data=gwres.df, proj4string=CRS(p4s), match.ID=F)
           gridded(SDF) <- griddedObj 
        }
    }
    else
        SDF <- SpatialPointsDataFrame(coords=rp.locat, data=gwres.df, proj4string=CRS(p4s), match.ID=F)
   ##############
    res<-list(GW.diagnostic=GW.diagnostic,glm.res=glm.res,SDF=SDF) 
}

############ Binomial GWGLM
gwr.binomial<-function(y,x,regression.points,W1.mat,W2.mat,hatmatrix,tol=1.0e-5, maxiter=20)
{
    p4s <- as.character(NA)
    if (is(regression.points, "Spatial"))
    {
      p4s <- proj4string(regression.points)
    }
    
    ############################################
    ##Generalized linear regression
    glm.res<-glm.fit(x, y, family = binomial())  
    ############################################
	  var.n<-ncol(x)
	  rp.locat<-coordinates(regression.points)
	  rp.n<-nrow(rp.locat)
    dp.n<-nrow(x)
    betas <-matrix(nrow=rp.n, ncol=var.n)
    betas1<- betas
    betas.SE <-matrix(nrow=rp.n, ncol=var.n)
    betas.TV <-matrix(nrow=rp.n, ncol=var.n)
    ##S: hatmatrix
    S<-matrix(nrow=dp.n,ncol=dp.n)
    #C.M<-matrix(nrow=dp.n,ncol=dp.n)
    colnames(betas) <- colnames(x)
    #colnames(betas)[1]<-"Intercept" 
    ####################################
    ##model calibration
    n=rep(1,length(y))
    it.count <- 0
    llik <- 0.0
    mu <- 0.5
    nu <- 0
    cat(" Iteration    Log-Likelihood\n=========================\n")
    wt2 <- rep(1,dp.n)
    repeat {
     y.adj <- nu + (y - n*mu)/(n*mu*(1 - mu))
     for (i in 1:dp.n)
     {
        W.i<-W1.mat[,i]
        gw.resi<-gw.reg(x,y.adj,W.i*wt2,hatmatrix=F,i)
        betas1[i,]<-gw.resi[[1]]
     }
     nu <- gwr.fitted(x,betas1)
     mu <- exp(nu)/(1 + exp(nu))  
     old.llik <- llik
     llik <- sum(lchoose(n,y) + (n-y)*log(1 - mu/dp.n) + y*log(mu/n))
     if(is.na(llik)) llik <-old.llik
     cat(paste("   ",formatC(it.count,digits=4,width=4),"    ",formatC(llik,digits=4,width=7),"\n"))
     if (abs((old.llik - llik)/llik) < tol) break
     wt2 <- n*mu*(1-mu) 
     it.count <- it.count+1
     if (it.count == maxiter) break}
     
     if(hatmatrix)
     {
        for (i in 1:rp.n)
        { 
          W.i<-W1.mat[,i]
          gw.resi<-gw.reg(x,y.adj,W.i*wt2,hatmatrix,i)
          betas[i,]<-gw.resi[[1]] ######See function by IG
          S[i,]<-gw.resi[[2]]
          Ci<-gw.resi[[3]]
          #betas.SE[i,]<-diag(Ci%*%t(Ci))
          invwt2 <- 1.0 /wt2
          betas.SE[i,] <- diag((Ci*invwt2) %*% t(Ci))   #see Nakaya et al. (2005)
        }
        tr.S<-sum(diag(S))
        tr.StS<-sum(S^2)
        ###edf is different from the definition in Chris' code
        edf<-dp.n-2*tr.S+tr.StS
        yhat<-gwr.fitted(x, betas)
        residual<-y-exp(yhat)/(1+exp(yhat))
        ########rss <- sum((y - gwr.fitted(x,b))^2)
        rss <- sum(residual^2)
        sigma.hat <- rss/edf
        sigma.aic <- rss/dp.n   ### can be omitted? (TN)
        for(i in 1:dp.n)
        {
           betas.SE[i,]<-sqrt(sigma.hat*betas.SE[i,])
           betas.TV[i,]<-betas[i,]/betas.SE[i,] 
        }
        #AICc <- -2*llik + 2*tr.S*dp.n/(dp.n-tr.S-2)
        AICc <- -2*llik + 2*tr.S + 2*tr.S*(tr.S+1)/(dp.n-tr.S-1)
        yss.g <- sum((y - mean(y))^2)
        gw.R2<-1-rss/yss.g; ##R Square valeu  ### is R2 needed? (TN)
        gwR2.adj<-1-(1-gw.R2)*(dp.n-1)/(edf-1) #Adjusted R squared value
        GW.diagnostic<-list(rss=rss,AICc=AICc,edf=edf,gw.R2=gw.R2,gwR2.adj=gwR2.adj)        
     }
     else
     {
        for (i in 1:rp.n)
        { 
          W.i<-W2.mat[,i]
          gw.resi<-gw.reg(x,y.adj,W.i*wt2,hatmatrix,i)
          betas[i,]<-gw.resi[[1]]
        }
     }
    if (hatmatrix)                                         
    {
      gwres.df<-data.frame(betas,y,exp(yhat)/(1+exp(yhat)),residual,betas.SE,betas.TV)
      colnames(gwres.df)<-c(c(c(colnames(betas),c("y","yhat","residual")),paste(colnames(betas), "SE", sep="_")),paste(colnames(betas), "TV", sep="_"))
    }
    else
    {
      gwres.df<-data.frame(betas)
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
           SDF <- SpatialPointsDataFrame(coords=rp.locat, data=gwres.df, proj4string=CRS(p4s), match.ID=F)
           gridded(SDF) <- griddedObj 
        }
    }
    else
        SDF <- SpatialPointsDataFrame(coords=rp.locat, data=gwres.df, proj4string=CRS(p4s), match.ID=F)
   ##############
    res<-list(GW.diagnostic=GW.diagnostic,glm.res=glm.res,SDF=SDF) 
}


############################Layout function for outputing the GWR results
##Author: BL	
print.ggwrm<-function(x, ...)
{
  if(class(x) != "ggwrm") stop("It's not a gwm object")
  cat("   ***********************************************************************\n")
  cat("   *                       Package   GWmodel                             *\n")
  cat("   ***********************************************************************\n")
  cat("   Program starts at:", as.character(x$timings$start), "\n")
  cat("   Call:\n")
  cat("   ")
  print(x$this.call)
  vars<-all.vars(x$GW.arguments$formula)
  var.n<-length(x$glm.res$coefficients)
	cat("\n   Dependent (y) variable: ",vars[1])
	cat("\n   Independent variables: ",vars[-1])
	dp.n<-length(x$glm.res$residuals)
	cat("\n   Number of data points:",dp.n)
	cat("\n   Used family:",x$GW.arguments$family)
	################################################################ Print Linear
 	cat("\n   ***********************************************************************\n")
	  cat("   *              Results of Generalized linear Regression               *\n")
	  cat("   ***********************************************************************\n")
	print(summary.glm(x$glm.res))	
	#########################################################################
	cat("\n   ***********************************************************************\n")
	  cat("   *          Results of Geographically Weighted Regression              *\n")
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
	
	cat("\n   ************Summary of Generalized GWR coefficient estimates:**********\n")      
		df0 <- as(x$SDF, "data.frame")[,1:var.n, drop=FALSE]
        if (any(is.na(df0))) {
            df0 <- na.omit(df0)
            warning("NAs in coefficients dropped")
        }
	CM <- t(apply(df0, 2, summary))[,c(1:3,5,6)]
	if(var.n==1) 
    { 
      CM <- matrix(CM, nrow=1)
      colnames(CM) <- c("Min.", "1st Qu.", "Median", "3rd Qu.", "Max.")
      rownames(CM) <- names(x$SDF)[1]
    }
	rnames<-rownames(CM)
		for (i in 1:length(rnames))
			 rnames[i]<-paste("   ",rnames[i],sep="")
	rownames(CM) <-rnames 
	printCoefmat(CM)
	cat("   ************************Diagnostic information*************************\n")

	if (x$GW.arguments$hatmatrix) 
  {	
		cat("   Number of data points:", dp.n, "\n")
		cat("   Effective degrees of freedom (n-2trace(S) + trace(S'S)):", x$GW.diagnostic$edf, "\n")
		cat("   AICc :",
                    x$GW.diagnostic$AICc, "\n")
		cat("   Residual sum of squares:", x$GW.diagnostic$rss, "\n")
    cat("   R-square value: ",x$GW.diagnostic$gw.R2,"\n")
		cat("   Adjusted R-square value: ",x$GW.diagnostic$gwR2.adj,"\n")	
  }
	cat("\n   ***********************************************************************\n")
	cat("   Program stops at:", as.character(x$timings$stop), "\n")
	invisible(x)
}

