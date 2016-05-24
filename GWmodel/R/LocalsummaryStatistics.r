###Calculate the local summary statistics
##GW means, variances (standard deviations), skew
##Using different distance metrics
##Author: Binbin Lu, Isabella Gollini
#Binbin - Could you do a bw.gwss() function just for the means and medians?
#It wouldn't take much to do - i.e. the local means and medians are the predictions for the cross-validation approach.
gwss <- function (data, summary.locat, vars, kernel = "bisquare", adaptive = FALSE, 
          bw, p = 2, theta = 0, longlat = F, dMat, quantile = FALSE) 
{
  findq <- function(x, w, p = c(0.25, 0.5, 0.75)) {
    lw <- length(w)
    lp <- length(p)
    q <- rep(0, lp)
    xo <- sort(x)
    wo <- w[order(x)]
    for (i in 1:lp) {
      cond <- max({
        cumsum(wo) <= p[i]
      } * seq(1:lw))
      if (cond == 0) 
        cond <- 1
      q[i] <- xo[cond]
    }
    q
  }
  S_Rho <- function(x, y, w) {
    n <- length(w)
    xr <- rank(x)
    yr <- rank(y)
    scorr <- 1 - {
      6 * w %*% {
       {xr - yr}^2
      }
    }/{
      n * {
        n^2 - 1
      }
    }
  }
  if (is(data, "Spatial")) {
    p4s <- proj4string(data)
    dp.locat <- coordinates(data)
  }
  else if (is(data, "data.frame") && (!missing(dMat))) 
    data <- data
  else stop("Given data must be a Spatial*DataFrame or data.frame object")
  if (missing(summary.locat)) {
    sp.given <- FALSE
    summary.locat <- data
    sp.locat <- coordinates(summary.locat)
  }
  else {
    sp.given <- T
    if (is(summary.locat, "Spatial")) 
      sp.locat <- coordinates(summary.locat)
    else {
      warning("Output loactions are not packed in a Spatial object,and it has to be a two-column numeric vector")
      summary.locat <- sp.locat
    }
  }
  data <- as(data, "data.frame")
  dp.n <- nrow(data)
  sp.n <- nrow(sp.locat)
  if (missing(dMat)) 
    DM.given <- F
  else {
    DM.given <- T
    dim.dMat <- dim(dMat)
    if (dim.dMat[1] != dp.n || dim.dMat[2] != sp.n) 
      stop("Dimensions of dMat are not correct")
  }
  if (missing(vars)) 
    stop("Variables input error")
  if (missing(bw) || bw <= 0) 
    stop("Bandwidth is not specified incorrectly")
  len.var <- length(vars)
  col.nm <- colnames(data)
  var.idx <- match(vars, col.nm)[!is.na(match(vars, col.nm))]
  if (length(var.idx) == 0) 
    stop("Variables input doesn't match with data")
  x <- data[, var.idx]
  x <- as.matrix(x)
  var.nms <- names(data)[var.idx]
  var.n <- ncol(x)
  if(len.var > var.n)
     warning("Invalid variables have been specified, please check them again!")
  
  local.mean <- matrix(numeric(var.n * sp.n), ncol = var.n)
  standard.deviation <- matrix(numeric(var.n * sp.n), ncol = var.n)
  local.skewness <- matrix(numeric(var.n * sp.n), ncol = var.n)
  LCV <- matrix(numeric(var.n * sp.n), ncol = var.n)
  LVar <- matrix(numeric(var.n * sp.n), ncol = var.n)
  if (quantile == TRUE) {
    local.median <- matrix(numeric(var.n * sp.n), ncol = var.n)
    IQR <- matrix(numeric(var.n * sp.n), ncol = var.n)
    QI <- matrix(numeric(var.n * sp.n), ncol = var.n)
  }
  cov.nms <- c()
  corr.nms <- c()
  cov.mat <- c()
  corr.mat <- c()
  s.corr.mat <- c()
  if (var.n > 1)
  {
     cov.mat <- matrix(numeric((var.n - 1) * var.n * sp.n/2), 
                    nrow = sp.n)
     corr.mat <- matrix(numeric((var.n - 1) * var.n * sp.n/2), 
                     nrow = sp.n)
      s.corr.nms <- c()
      s.corr.mat <- matrix(numeric((var.n - 1) * var.n * sp.n/2), 
                       nrow = sp.n)
  }
  for (i in 1:sp.n) {
    if (DM.given) 
      dist.vi <- dMat[, i]
    else {
      if (sp.given) 
        dist.vi <- gw.dist(dp.locat, sp.locat, focus = i, 
                           p, theta, longlat)
      else dist.vi <- gw.dist(dp.locat = dp.locat, focus = i, 
                              p = p, theta = theta, longlat = longlat)
    }
    W.i <- matrix(gw.weight(dist.vi, bw, kernel, adaptive), 
                  nrow = 1)
    sum.w <- sum(W.i)
    Wi <- W.i/sum.w
    local.mean[i, ] <- Wi %*% x
    if (quantile == TRUE) {
      kecor <- cor(x, method = "kendall")
      quant <- apply(x, 2, findq, w = c(Wi))
      local.median[i, ] <- quant[2, ]
      IQR[i, ] <- quant[3, ] - quant[1, ]
      QI[i, ] <- {
        2 * quant[2, ] - quant[3, ] - quant[1, ]   # quant[2,] changed to quant[1,] in this line...
      }/IQR[i, ]
    }
    for (j in 1:var.n) {
      LVar[i, j] <- Wi %*% ((x[, j] - local.mean[i, j])^2)
      standard.deviation[i, j] <- sqrt(LVar[i, j])
      local.skewness[i, j] <- (Wi %*% ((x[, j] - local.mean[i, 
                                                            j])^3))/(standard.deviation[i, j]^3) # Exponent in denominator is 3, not 1.5
      LCV[i, j] <- standard.deviation[i, j]/local.mean[i, 
                                                       j]
    }
    if (var.n >= 2) {
      tag <- 0
      for (j in 1:(var.n - 1)) for (k in (j + 1):var.n) {
        tag <- tag + 1
        cov.mat[i,tag] <- cov.wt(cbind(x[, j],x[, k]), wt=Wi[1,])$cov[1,2] # Replace the old code for something faster
        corr.mat[i, tag] <- cov.wt(cbind(x[, j],x[, k]), wt=Wi[1,],cor=TRUE)$cor[1,2] # Also replaced this as this way is a bit faster
        s.corr.mat[i, tag] <- cov.wt(cbind(rank(x[, j]), rank(x[, k])), wt=Wi[1,],cor=TRUE)$cor[1,2] # Had to replace Spearmans Rho function - supplied one not correct for weighted version
      }
    }
  }
  colnames(local.mean) <- paste(var.nms, "LM", sep = "_")
  colnames(standard.deviation) <- paste(var.nms, "LSD", sep = "_")
  colnames(LVar) <- paste(var.nms, "LVar", sep = "_")
  colnames(local.skewness) <- paste(var.nms, "LSKe", sep = "_")
  colnames(LCV) <- paste(var.nms, "LCV", sep = "_")
  if (quantile == TRUE) {
    colnames(local.median) <- paste(var.nms, "Median", sep = "_")
    colnames(IQR) <- paste(var.nms, "IQR", sep = "_")
    colnames(QI) <- paste(var.nms, "QI", sep = "_")
  }
  if (var.n>1)
  {
    for (i in 1:(var.n - 1)) {
      for (j in (i + 1):var.n) {
        cov.v1v2 <- paste("Cov", paste(var.nms[i], var.nms[j], 
                                       sep = "."), sep = "_")
        corr.v1v2 <- paste("Corr", paste(var.nms[i], var.nms[j], 
                                         sep = "."), sep = "_")
        cov.nms <- c(cov.nms, cov.v1v2)
        corr.nms <- c(corr.nms, corr.v1v2)
        s.corr.v1v2 <- paste("Spearman_rho", paste(var.nms[i], 
                                                   var.nms[j], sep = "."), sep = "_")
        s.corr.nms <- c(s.corr.nms, s.corr.v1v2)
      }
    }
    colnames(cov.mat) <- cov.nms
    colnames(corr.mat) <- corr.nms
    colnames(s.corr.mat) <- s.corr.nms
  }
  if (quantile == TRUE)
  {
    if(var.n >1) 
       res.df <- data.frame(local.mean, standard.deviation, 
                         LVar, local.skewness, LCV, cov.mat, corr.mat, s.corr.mat, 
                         local.median, IQR, QI)
    else 
       res.df <- data.frame(local.mean, standard.deviation, 
                         LVar, local.skewness, LCV, 
                         local.median, IQR, QI)
  }
  else 
  {
    if(var.n >1) 
       res.df <- data.frame(local.mean, standard.deviation, 
                            LVar, local.skewness, LCV, cov.mat, corr.mat, s.corr.mat)
    else
       res.df <- data.frame(local.mean, standard.deviation, 
                            LVar, local.skewness, LCV)
  }
  rownames(res.df) <- rownames(sp.locat)
  griddedObj <- F
  if (is(summary.locat, "Spatial"))
  { 
    if (is(summary.locat, "SpatialPolygonsDataFrame")) 
    {
      polygons <- polygons(summary.locat)
      SDF <- SpatialPolygonsDataFrame(Sr = polygons, data = res.df, 
                                      match.ID = F)
    }
    else
    {
      griddedObj <- gridded(summary.locat)
      SDF <- SpatialPointsDataFrame(coords = sp.locat, data = res.df, 
                                       proj4string = CRS(p4s), match.ID=F)
      gridded(SDF) <- griddedObj
    }
  }
  else
    SDF <- SpatialPointsDataFrame(coords = sp.locat, data = res.df, 
                                       proj4string = CRS(p4s), match.ID=F) 
  res <- list(SDF = SDF, vars = vars, kernel = kernel, adaptive = adaptive, 
              bw = bw, p = p, theta = theta, longlat = longlat, DM.given = DM.given, 
              sp.given = sp.given, quantile = quantile)
  class(res) <- "gwss"
  invisible(res)
}

###print the summary information for local summary statistics
#Author: Binbin Lu, Isabella Gollini
print.gwss<-function(x, ...)
{
    if (class(x) != "gwss") 
        stop("It's not a lss object")
    cat("   ***********************************************************************\n")
    cat("   *                       Package   GWmodel                             *\n")
    cat("   ***********************************************************************\n")
    cat("\n   ***********************Calibration information*************************\n")
    vars <- x$vars
    var.n <- length(vars)
    cat("\n   Local summary statistics calculated for variables:")
    cat("\n   ", vars)
    dp.n <- nrow(data.frame(x$SDF))
    cat("\n   Number of summary points:", dp.n)
    cat("\n   Kernel function:", x$kernel, "\n")
    if (x$sp.given) 
        cat("   Summary points: A seperate set of summary points is used.\n")
    else cat("   Summary points: the same locations as observations are used.\n")
    if (x$adaptive) 
        cat("   Adaptive bandwidth: ", x$bw, " (number of nearest neighbours)\n", 
            sep = "")
    else cat("   Fixed bandwidth:", x$bw, "\n")
    if (x$DM.given) 
        cat("   Distance metric: A distance matrix is specified for this model calibration.\n")
    else {
        if (x$longlat) 
            cat("   Distance metric: Great Circle distance metric is used.\n")
        else if (x$p == 2) 
            cat("   Distance metric: Euclidean distance metric is used.\n")
        else if (x$p == 1) 
            cat("   Distance metric: Manhattan distance metric is used.\n")
        else if (is.infinite(x$p)) 
            cat("   Distance metric: Chebyshev distance metric is used.\n")
        else cat("   Distance metric: A generalized Minkowski distance metric is used with p=", 
            x$p, ".\n")
        if (x$theta != 0 && x$GW.agruments$p != 2 && !x$longlat) 
            cat("   Coordinate rotation: The coordinate system is rotated by an angle", 
                x$theta, "in radian.\n")
    }
    cat("\n   ************************Local Summary Statistics:**********************\n")
    df0 <- as(x$SDF, "data.frame")
    cat("   Summary information for Local means:\n")
    nms.LM <- paste(vars, "LM", sep = "_")
    df.lm <- df0[, nms.LM]
    if (var.n ==1)
      dim(df.lm) <- c(dp.n,var.n)
    #LM <- t(apply(df.lm, 2, summary))[, c(1:3, 5, 6)]
    LM <- t(apply(df.lm, 2, summary))[, c(1:3, 5, 6)]
    if (var.n ==1 )
    {
      cat(nms.LM[1],"\n")
      print(LM)
    }
    else
    {
      #dim(LM) <- c(var.n, 5)
      rnames <- rownames(LM)
      for (i in 1:length(rnames)) rnames[i] <- paste("   ", rnames[i], 
          sep = "")
      rownames(LM) <- rnames
      printCoefmat(LM)
    }
    cat("   Summary information for local standard deviation :\n")
    nms.sd <- paste(vars, "LSD", sep = "_")
    df.sd <- df0[, nms.sd]
    if (var.n ==1 )
        dim(df.sd) <- c(dp.n,var.n)
    SD <- t(apply(df.sd, 2, summary))[, c(1:3, 5, 6)]
    if (var.n ==1 )
    {
      cat(nms.sd[1],"\n")
      print(SD)
    }
    else
    {
      rnames <- rownames(SD)
      for (i in 1:length(rnames)) rnames[i] <- paste("   ", rnames[i], 
          sep = "")
      rownames(SD) <- rnames
      printCoefmat(SD)
    }
    cat("   Summary information for local variance :\n")
    nms.lvar <- paste(vars, "LVar", sep = "_")
    df.lvar <- df0[, nms.lvar]
    if (var.n ==1 )
       dim(df.lvar) <- c(dp.n,var.n)
    LVar <- t(apply(df.lvar, 2, summary))[, c(1:3, 5, 6)]
    if (var.n ==1 )
    {
      cat(nms.lvar[1],"\n")
      print(LVar)
    }
    else
    {
      rnames <- rownames(LVar)
      for (i in 1:length(rnames)) rnames[i] <- paste("   ", rnames[i], 
          sep = "")
      rownames(LVar) <- rnames
      printCoefmat(LVar)
    }
    cat("   Summary information for Local skewness:\n")
    nms.ske <- paste(vars, "LSKe", sep = "_")
    df.ske <- df0[, nms.ske]
    if (var.n ==1 )
       dim(df.ske) <- c(dp.n,var.n)
    SKE <- t(apply(df.ske, 2, summary))[, c(1:3, 5, 6)]
    if (var.n ==1 )
    {
      cat(nms.ske[1],"\n")
      print(SKE)
    }
    else
    {
      rnames <- rownames(SKE)
      for (i in 1:length(rnames)) rnames[i] <- paste("   ", rnames[i], 
          sep = "")
      rownames(SKE) <- rnames
      printCoefmat(SKE)
    }
    cat("   Summary information for localized coefficient of variation:\n")
    nms.lcv <- paste(vars, "LCV", sep = "_")
    df.lcv <- df0[, nms.lcv]
    if (var.n ==1 )
       dim(df.lcv) <- c(dp.n,var.n)
    LCV <- t(apply(df.lcv, 2, summary))[, c(1:3, 5, 6)]
    if (var.n ==1 )
    {
      cat(nms.lcv[1],"\n")
      print(LCV)
    }
    else
    {
      rnames <- rownames(LCV)
      for (i in 1:length(rnames)) rnames[i] <- paste("   ", rnames[i], 
          sep = "")
      rownames(LCV) <- rnames
      printCoefmat(LCV)
    }
    if (var.n >= 2) {
        cov.nms <- c()
        corr.nms <- c()
        s.corr.nms <- c()
        #k.corr.nms <- c()
        for (i in 1:(var.n - 1)) {
            for (j in (i + 1):var.n) {
                cov.v1v2 <- paste("Cov", paste(vars[i], vars[j], 
                  sep = "."), sep = "_")
                corr.v1v2 <- paste("Corr", paste(vars[i], vars[j], 
                  sep = "."), sep = "_")
                 s.corr.v1v2 <- paste("Spearman_rho", paste(vars[i], vars[j], 
                  sep = "."), sep = "_")
                # k.corr.v1v2 <- paste("Kendall_tau", paste(vars[i], vars[j], sep = "."), sep = "_")
                cov.nms <- c(cov.nms, cov.v1v2)
                corr.nms <- c(corr.nms, corr.v1v2)
           		s.corr.nms <- c(s.corr.nms, s.corr.v1v2)
               # k.corr.nms <- c(k.corr.nms, k.corr.v1v2)
                
            }
        }
        cat("   Summary information for localized Covariance and Correlation between these variables:\n")
        
        df.covar <- data.frame(df0[, cov.nms], df0[, corr.nms],df0[, s.corr.nms])
        Covar <- t(apply(df.covar, 2, summary))[, c(1:3, 5, 6)]
        rownames(Covar) <- c(cov.nms, corr.nms,s.corr.nms)
        
        rnames <- rownames(Covar)
        for (i in 1:length(rnames)) rnames[i] <- paste("   ", 
            rnames[i], sep = "")
        rownames(Covar) <- rnames
        printCoefmat(Covar)
   # cat("   Summary information for localized Correlation between these variables:\n")
#    
#    df.corr<-
#    Corr <- t(apply(df.corr, 2, summary))[,c(1:3,5,6)]
#    rnames<-rownames(Corr)
#    	for (i in 1:length(rnames))
#    		 rnames[i]<-paste("   ",rnames[i],sep="")
#    rownames(Corr) <-rnames 
#    printCoefmat(Corr)
     
     if(x$quantile==TRUE){ 
	cat("   Summary information for Local median:\n")
    nms.LMed <- paste(vars, "Median", sep = "_")
    df.lmed <- df0[, nms.LMed]
    LMed <- t(apply(df.lmed, 2, summary))[, c(1:3, 5, 6)]
    rnames <- rownames(LMed)
    for (i in 1:length(rnames)) rnames[i] <- paste("   ", rnames[i], 
        sep = "")
    rownames(LMed) <- rnames
    printCoefmat(LMed)
        	cat("   Summary information for Interquartile range:\n")
    nms.IQR <- paste(vars, "IQR", sep = "_")
    df.IQR <- df0[, nms.IQR]
    sIQR <- t(apply(df.IQR, 2, summary))[, c(1:3, 5, 6)]
    rnames <- rownames(sIQR)
    for (i in 1:length(rnames)) rnames[i] <- paste("   ", rnames[i], 
        sep = "")
    rownames(sIQR) <- rnames
    printCoefmat(sIQR)
        	cat("   Summary information for Quantile imbalance:\n")
        nms.QI <- paste(vars, "QI", sep = "_")
    df.QI <- df0[, nms.QI]
    sQI <- t(apply(df.QI, 2, summary))[, c(1:3, 5, 6)]
    rnames <- rownames(sQI)
    for (i in 1:length(rnames)) rnames[i] <- paste("   ", rnames[i], 
        sep = "")
    rownames(sQI) <- rnames
    printCoefmat(sQI)
        
       }

  }
	cat("\n   ************************************************************************\n")
	invisible(x)
}

#  Local correlation calculator - can do adaptive and non-adaptive bandwidths.
#Author: Binbin Lu
local.corr<-function(x,dp.locat, sp.locat,SD)
{
  var.nms<-colnames(x)
  var.n<-length(var.nms)
  dp.n<-nrow(dp.locat)
  sp.n<-nrow(sp.locat)
  
  for (i in 1:(var.n-1))
  {
    for(j in (i+1):var.n)
    {
      cov.v1v2<-paste("Cov",paste(var.nms[i],var.nms[j],sep="."),sep="_")
      corr.v1v2<-paste("Corr",paste(var.nms[i],var.nms[j],sep="."),sep="_")
      cov.nms<-c(cov.nms,cov.v1v2)
      corr.nms<-c(corr.nms, corr.v1v2)
      
    }
  }
}     


# Randomisation Tests for GWSS
montecarlo.gwss<-function(data, vars, kernel = "bisquare", 
                 adaptive = FALSE, bw, p = 2, theta = 0, longlat = F, 
                 dMat, quantile=FALSE,nsim=99) 
{
	if (is(data, "Spatial")) {
        p4s <- proj4string(data)
        dp.locat <- coordinates(data)
    }
    else if (is(data, "data.frame") && (!missing(dMat))) 
        data <- data
    else stop("Given data must be a Spatial*DataFrame or data.frame object")
    data <- as(data, "data.frame")
    dp.n <- nrow(data)
    if (missing(dMat))
    { 
        DM.given <- F
        dMat <- gw.dist(dp.locat=dp.locat, p=p, theta=theta, longlat=longlat)
    }
    else {
        DM.given <- T
        dim.dMat <- dim(dMat)
        if (dim.dMat[1] != dp.n || dim.dMat[2] != dp.n) 
            stop("Dimensions of dMat are not correct")
    }
    if (missing(vars)) 
        stop("Variables input error")
    if (missing(bw) || bw <= 0) 
        stop("Bandwidth is not specified incorrectly")
    col.nm <- colnames(data)
    var.idx <- match(vars, col.nm)[!is.na(match(vars, col.nm))]
    if (length(var.idx) == 0) 
        stop("Variables input doesn't match with data")
    x <- data[, var.idx]
    x <- as.matrix(x)
    var.nms <- colnames(x)
    var.n <- ncol(x)
    colnames(x) <- vars

 
	dataI<-SpatialPointsDataFrame(dp.locat,data.frame(x))

	resi<-(gwss(data=dataI, vars=vars, kernel=kernel, adaptive=adaptive, bw=bw, p=p, theta=theta, 
                      longlat=longlat, dMat=dMat, quantile=quantile)$SDF)@data
	
	nstats<-ncol(resi)
	
	lss.i<-array(0,dim=c(dp.n,nstats,nsim+1))

	lss.i[,,1]<- as.matrix(resi)
	dMat1<- dMat
	for(i in 1:nsim)
	{
		mcs <- sample(dp.n)   	
		dataI<-SpatialPointsDataFrame(dp.locat[mcs,],data.frame(x))
		dMat1[mcs,]<-dMat[1:dp.n,]
		dMat1[,mcs]<-dMat[,1:dp.n]
		resi<-(gwss(data=dataI, vars=vars, kernel=kernel, adaptive=adaptive, 
                       bw=bw, p=p, theta=theta, longlat=longlat, dMat=dMat1, quantile=quantile)$SDF)@data		
		lss.i[,,i+1]<-as.matrix(resi)
	}
	
	dimnames(lss.i)[[1]]<-seq(0,dp.n-1)
	dimnames(lss.i)[[2]]<-names(resi)
	dimnames(lss.i)[[3]]<-c('Original_Data',paste('Sample',seq(1:nsim),sep='_'))	
	test<-apply(lss.i,c(1,2),function(x,n) 1 - rank(x,ties.method='first')[1]/n,n=nsim+1)	
	test
}
                                    