####GW Discriminant Analysis
#The GWDA technique proposed here exploits the fact that linear and quadratic
#discriminant analyses (LDA and QDA) rely only on the mean vector and covariance
#matrix of {x} for each population (the former assumes that the covariance is the same for all m populations).

gwda <- function(formula, data, predict.data,validation = T, COV.gw=T, 
                 mean.gw=T, prior.gw=T, prior=NULL, wqda =F,
                kernel = "bisquare", adaptive = FALSE, bw,
                 p = 2, theta = 0, longlat = F,dMat)
{
  ##########################
  this.call <- match.call()
  p4s <- as.character(NA)
  ############################training data
  #data must be given as training data 
  if (is(data, "Spatial")) {
    p4s <- proj4string(data)
    dp.locat <- coordinates(data)
  }
  else
     stop("Given training data must be a Spatial*DataFrame or data.frame object")
  
  
  ##########################prediction data
  #if it is not given, the training data will be also predicted using cross-validation appraoch
  if(missing(predict.data))
  {
     cv.predict <- T
     pr.locat <- dp.locat
     predict.data <- data
     pr.data <- as(predict.data, "data.frame")
  }
  else
  {
     cv.predict <- F
     if (is(predict.data, "Spatial")) 
     {
        pr.locat <- coordinates(predict.data)
        pr.data <- as(predict.data, "data.frame")
     }
     else
         stop("Prediction data must be a Spatial*DataFrame or data.frame object")
  }
  data <- as(data, "data.frame")
  ##########################validation data
  ##As default, the training data will be also used as validation data
#  if (is(validat.data, "Spatial")) {
#    vp.locat <- coordinates(validat.data)
#    va.data <- as(validat.data, "data.frame")
#  }
#  else
#     stop("Validation data must be a Spatial*DataFrame or data.frame object")
  ######## variables from formula
  vars <- all.vars(formula)
  grouping.nm <- vars[1]
  expl.vars <- vars[-1]
  m <- length(expl.vars)
  if (m<2)
     stop("Two or more variables shoule be specfied for analysis")
  #x, y from training data
  res1 <-grouping.xy(data, grouping.nm, expl.vars)
  x<- res1$x
  grouping <- res1$y 
  lev <- levels(grouping)
  #x for prediction data
  x.pr <- grouping.xy(data=pr.data, expl.vars=expl.vars)$x
  #######Distance matrix for training
  dp.n <- nrow(dp.locat)
  pr.n <- nrow(pr.locat)
  if (missing(dMat))
  {
     dMat <- gw.dist(dp.locat=dp.locat, rp.locat=pr.locat, p=p, theta=theta, longlat=longlat)
  }
  else
  {
    dim.dMat<-dim(dMat)
    if (dim.dMat[1]!=dp.n||dim.dMat[2]!=pr.n)
    stop ("Dimensions of dMat are not correct")
  }
  
   ##Weighting matrix for 
   wt <- gw.weight(dMat,bw,kernel,adaptive)
  ##########  prior probility
  if(!is.null(prior))
  {
     prior.given <- T
     if(any(prior < 0) || round(sum(prior), 5) != 1) stop("invalid 'prior'")
  }
  else
     prior.given <- F
  ####if the prediction data is not given, cross-validation appraoch will be used for prediction
  if(cv.predict) diag(wt) <- 0
  #####Weighted qda
  if(wqda) 
     res.df <- wqda(x, grouping, x.pr, wt, COV.gw, 
                 mean.gw, prior.gw, prior)
  else
     res.df <- wlda(x, grouping, x.pr, wt, COV.gw, 
                 mean.gw, prior.gw, prior)
  ###For validation
  correct.ratio <- 0
  if(validation)
  {
    obs.grouping <- grouping.xy(pr.data, grouping.nm, expl.vars)$y
    n.correct <- length(which(obs.grouping == res.df[,"group.predicted"]))
    correct.ratio <- n.correct/nrow(pr.locat)
  }
  #####output in a spatial*dataframe#
    res.df <- data.frame(res.df)
    rownames(res.df) <- rownames(pr.locat)
    griddedObj <- F
    if (is(predict.data, "Spatial"))
    { 
        if (is(predict.data, "SpatialPolygonsDataFrame"))
        {
          polygons <- polygons(predict.data)
          SDF <- SpatialPolygonsDataFrame(Sr = polygons, data = res.df, 
                                    match.ID = F)
        }
        else
        {
           griddedObj <- gridded(predict.data)
           SDF <- SpatialPointsDataFrame(coords = pr.locat, data = res.df, 
                                     proj4string = CRS(p4s), match.ID=F)
           gridded(SDF) <- griddedObj 
        }
    }
    else
        SDF <- SpatialPointsDataFrame(coords = pr.locat, data = res.df, 
                                     proj4string = CRS(p4s), match.ID=F)   
#   if (is(predict.data, "SpatialPolygonsDataFrame")) {
#    polygons <- polygons(predict.data)
#    SDF <- SpatialPolygonsDataFrame(Sr = polygons, data = res.df, 
#                                    match.ID = F)
#  }
#  else SDF <- SpatialPointsDataFrame(coords = pr.locat, data = res.df, 
#                                     proj4string = CRS(p4s), match.ID=F)                                     
  res<-list()
  res$SDF<- SDF
  res$this.call <- this.call
  res$grouping.nm <- grouping.nm
  res$lev <- lev
  res$expl.vars <- expl.vars
  res$cv.predict <- cv.predict
  res$mean.gw <- mean.gw
  res$COV.gw <- COV.gw
  res$prior.given <- prior.given
  res$prior.gw <- prior.gw        
  res$kernel <- kernel
  res$adaptive <- adaptive
  res$bw <- bw
  res$p = p
  res$theta = theta 
  res$longlat = longlat
  res$wqda <- wqda 
  res$validation <- validation
  res$pr.n <- nrow(pr.locat)
  res$correct.ratio <- correct.ratio
  res
  class(res) <-"gwda"
  invisible(res) 
}
############################Layout function for outputing the GWDA results
##Author: BL	
print.gwda<-function(x, ...)
{
  if(class(x) != "gwda") stop("It's not a gwda object")
  cat("   ***********************************************************************\n")
  cat("   *                       Package   GWmodel                             *\n")
  cat("   ***********************************************************************\n")
  cat("   Call:\n")
  cat("   ")
  print(x$this.call)
	cat("\n   Grouping factor: ",x$grouping.nm, " with the following groups: \n")
	cat("\n   ",x$lev)
	cat("\n    Discriminators: ",x$expl.vars)
	if(x$cv.predict) cat("\n    Prediction: No prediction data is given and leave-one-out cross-validation will be applied")
	else  cat("\n    Prediction: Ordinary prediction is made with given prediction data")
	if(x$mean.gw) cat("\n    Meams: Localised mean is used for GW discriminant analysis")
	else  cat("\n    Meams: Global means is used for GW discriminant analysis")
	if(x$COV.gw) cat("\n    Variance-covariance: Localised variance-covariance matrix is used for GW discriminant analysis")
	else  cat("\n    Variance-covariance: Global variance-covariance matrix is used for GW discriminant analysis")
	if(x$prior.given) cat("\n    Prior probability is pre-set for GW discriminant analysis")
  else
  {
    if(x$prior.gw) cat("\n    Localised prior probability is used for GW discriminant analysis")
	  else  cat("\n    Global prior probability is used for GW discriminant analysis")
  }
  if(x$adaptive)
	   cat("\n    Adaptive bandwidth: ", x$bw, " (number of nearest neighbours)\n", sep="") 
  else
     cat("\n    Fixed bandwidth:", x$bw, "\n")
  if (x$longlat)
    cat("    Distance metric: Great Circle distance metric is used.\n")
  else if (x$p==2)
    cat("    Distance metric: Euclidean distance metric is used.\n")
  else if (x$p==1)
    cat("    Distance metric: Manhattan distance metric is used.\n") 
  else if (is.infinite(x$p))
    cat("    Distance metric: Chebyshev distance metric is used.\n")
  else 
    cat("    Distance metric: A generalized Minkowski distance metric is used with p=",x$p,".\n")
  if (x$theta!=0&&x$p!=2&&!x$longlat)
    cat("    Coordinate rotation: The coordinate system is rotated by an angle", x$theta, "in radian.\n")
  if(x$validation)
  {
    cat("    The correct ratio is validated as ", x$correct.ratio) 
    cat("\n    The number of points for prediction is ", x$pr.n)        
  }  
	cat("\n   ***********************************************************************\n")
	invisible(x)
}
#Return the x (and y) data for grouping data
grouping.xy <- function(data, grouping.nm, expl.vars)
{
  y <- factor()
  col.nm <- colnames(data)
  if(!missing(grouping.nm)) 
  {
    idx <- match(grouping.nm, col.nm)
    y <- as.factor(data[,idx])
  }
  idx <- match(expl.vars, col.nm)
  n <- length(idx)
  x <- as.matrix(data[,idx],ncol=n)
  res <- list()
  res$x <- x
  res$y <- y
  res 
}
###Weighted qda function
wqda <- function(x, grouping, x.pr, wt, COV.gw=T, 
                 mean.gw=T, prior.gw=T, prior=NULL)
{
  lev <- levels(grouping)
  m <- length(lev)
  dp.n <- nrow(x)
  pr.n <- nrow(x.pr)
  
  wt.ones <- matrix(rep(1, dp.n*pr.n),ncol=pr.n)
  x.g <- splitx(x, grouping)
  local.mean <- list()
  ##variance-covariance matrix
  sigma.gw <- list()
  ##prior proportion
  prior.given <- F
  if(is.null(prior))
    prior <- list()
  else
  {
    prior.given <- T
    if(is.numeric(prior) && length(prior)==m && sum(prior)==1)
    {
      prior1 <- list()
      for(i in 1:m)
        prior1[[lev[i]]] <- rep(prior[i], pr.n)
    }
    prior <- prior1    
  }
  ##################################################################
  #Training
  for (i in 1:m)
  {
     xi <- x.g[[lev[i]]]
     idx <- as.numeric(rownames(xi))
     if(mean.gw)
       wti <- wt[idx,]
     else
       wti <- wt.ones[idx,]
     ####weighted means for each population i
      local.mean[[lev[i]]] <- wmean(xi, wti)
     ####weighted variances for each population i
     if(COV.gw)
       wti <- wt[idx,]
     else
       wti <- wt.ones[idx,]
     sigma.gw[[lev[i]]] <- wvarcov(xi, wti)
     if(!prior.given)
     {
        if(prior.gw)
        {
          wti <- wt[idx,]
          sum.w <- sum(wt)
        }
        else
        {
          wti <- wt.ones[idx,]
          sum.w <- sum(wt.ones)
        }
        prior[[lev[i]]] <- wprior(wti, sum.w)
     }
  }
  ##################################################################
  #prediction
  log.pf <- matrix(numeric(pr.n*m), ncol=m) 
  for(i in 1:m)
    for(j in 1:pr.n)
    {
      x.prj <- as.matrix(x.pr[j,],nrow=1)
      meani <- as.matrix(local.mean[[lev[i]]][j,],nrow=1)
      cov.matj <- sigma.gw[[lev[i]]][j,,]
      
      log.pf[j,i] <- (m/2)*log(norm(cov.matj)) + 0.5 *t(x.prj - meani)%*%solve(cov.matj)%*% (x.prj - meani) - log(prior[[lev[i]]][j])
    }
  colnames(log.pf) <- paste(lev, "logp", sep="_")
  group.pr <- vector("character", pr.n)
  for(i in 1:pr.n)
    group.pr[i] <- lev[which.min(log.pf[i,])[1]]
  res.df <- cbind(log.pf, group.pr)
  colnames(res.df) <-c(colnames(log.pf), "group.predicted") 
  res.df
}

###Weighted lda function
wlda <- function(x, grouping, x.pr, wt, COV.gw=T, 
                 mean.gw=T, prior.gw=T, prior=NULL)
{
  lev <- levels(grouping)
  m <- length(lev)
  dp.n <- nrow(x)
  pr.n <- nrow(x.pr)
  var.n <- ncol(x)
  wt.ones <- matrix(rep(1, dp.n*pr.n),ncol=pr.n) 
  x.g <- splitx(x, grouping)
  local.mean <- list()
  ##variance-covariance matrix
  sigma.gw <- list()
  ##prior proportion
  prior.given <- F
  if(is.null(prior))
    prior <- list()
  else
  {
    prior.given <- T
    if(is.numeric(prior) && length(prior)==m && sum(prior)==1)
    {
      prior1 <- list()
      for(i in 1:m)
        prior1[[lev[i]]] <- rep(prior[i], pr.n)
    }
    prior <- prior1    
  }
  ##################################################################
  #Training
  for (i in 1:m)
  {
     xi <- x.g[[lev[i]]]
     idx <- as.numeric(rownames(xi))
     if(mean.gw)
       wti <- wt[idx,]
     else
       wti <- wt.ones[idx,]
     ####weighted means for each population i
     local.mean[[lev[i]]] <- wmean(xi, wti)
     ####weighted variances for each population i
     if(COV.gw)
       wti <- wt[idx,]
     else
       wti <- wt.ones[idx,]
     sigma.gw[[lev[i]]] <- wvarcov(xi, wti)
     if(!prior.given)
     {
        if(prior.gw)
        {
          wti <- wt[idx,]
          sum.w <- sum(wt)
        }
        else
        {
          wti <- wt.ones[idx,]
          sum.w <- sum(wt.ones)
        }
        prior[[lev[i]]] <- wprior(wti, sum.w)
     }
  }
  sigma1.gw <- array(0, dim=c(pr.n,var.n,var.n)) 
  counts <- as.vector(table(grouping))
  for (i in 1:pr.n)
  {
    sigmai <- array(0, dim=c(var.n,var.n))
    for(j in 1:m)
    {
      sigmai <- sigmai + counts[j]*sigma.gw[[lev[j]]][i,,]
    }
    sigma1.gw[i,,] <- sigmai/sum(counts)
  }
  ##################################################################
  #prediction
  log.pf <- matrix(numeric(pr.n*m), ncol=m) 
  for(i in 1:m)
    for(j in 1:pr.n)
    {
      x.prj <- as.matrix(x.pr[j,],nrow=1)
      meani <- as.matrix(local.mean[[lev[i]]][j,],nrow=1)
      cov.matj <- sigma1.gw[j,,]
      log.pf[j,i] <- (m/2)*log(norm(cov.matj)) + 0.5 *t(x.prj - meani)%*%solve(cov.matj)%*%(x.prj - meani) - log(prior[[lev[i]]][j])
    }
  colnames(log.pf) <- paste(lev, "logp", sep="_")
  group.pr <- vector("character", pr.n)
  for(i in 1:pr.n)
    group.pr[i] <- lev[which.min(log.pf[i,])[1]]
  res.df <- cbind(log.pf, group.pr)
  colnames(res.df) <-c(colnames(log.pf), "group.predicted") 
  res.df
}

##Split X into groupes
splitx <- function(x, grouping)
{
  lev <- levels(grouping)
  p <- length(lev)
  counts <- as.vector(table(grouping))
  names(counts) <- lev
  if(any(counts < p+1)) stop("Some group is too small for training")
  res <- list()
  for(i in 1:p)
  {
    idx <- which(grouping==lev[i])
    xi <- x[idx,]
    rownames(xi) <- as.character(idx)
    res[[lev[i]]] <- xi
  }
  res
}
##Weighted means
wmean <- function(x, wt)
{
  var.n <- ncol(x)
  dp.n <- nrow(x)
  pr.n <- ncol(wt)
  local.mean <- matrix(numeric(var.n * pr.n), ncol = var.n)
  for (i in 1:pr.n)
  {
    w.i <- matrix(wt[,i],nrow = 1)
    sum.w <- sum(w.i)
    w.i <- w.i / sum.w
    local.mean[i, ] <- w.i %*% x
  }
  local.mean
}
##Weighted variance-covariance matrix
wvarcov <- function(x, wt)
{
  var.n <- ncol(x)
  dp.n <- nrow(x)
  pr.n <- ncol(wt)
  cov.mat <- array(0, dim=c(pr.n,var.n,var.n))
  for (i in 1:pr.n)
  {
    w.i <- wt[,i]
    sum.w <- sum(w.i)
    w.i <- w.i / sum.w
    cov.mat[i, ,] <-cov.wt(x, wt=w.i)$cov
  }
  cov.mat
}
##Weighted means
wprior <- function(wt, sum.w)
{
  pr.n <- ncol(wt)
  local.prior <- numeric(pr.n)
  for (i in 1:pr.n)
    local.prior[i] <- sum(wt[,i])/sum.w
  local.prior
}


##Confusion matrix
confusion.matrix <- function(original, classified)
{
  classes <- levels(original)
  n <- length(classes)
  cf.mat <- matrix(nrow=n+1, ncol = n+1)
  total <- 0
  for (i in 1:n)
  {
    tag1 <- 0
    for (j in 1:n)
    {
       cf.mat[i,j] <- length(which(is.na(match(which(original == classes[j]), which(classified == classes[i])))==FALSE))
       total <- total + cf.mat[i,j]
       tag1 <- tag1 + cf.mat[i,j]
    }
    cf.mat[i,n+1] <- tag1 
  }
  cf.mat[n+1,n+1] <- total
  for (i in 1:n)
     cf.mat[n+1,i] <- sum(cf.mat[1:n,i]) 
  rownames(cf.mat) <- colnames(cf.mat) <- c(classes, "Total")
  cf.mat
}