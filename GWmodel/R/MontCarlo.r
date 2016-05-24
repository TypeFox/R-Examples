#####Mont Carlo simulation
##Change the positions of a sequence randomly
#Author: Binbin Lu

montecarlo.gwr<-function(formula, data = list(),nsims=99, kernel="bisquare",adaptive=F, bw,
                         p=2, theta=0, longlat=F,dMat)
{
  ##Extract the model data frame
  this.call <- match.call() 
  if (!is.null(data))
  {
    if (is(data, "Spatial"))
    {
       dp.locat <- coordinates(data)
       data <- as(data, "data.frame")
       dp.n<-nrow(dp.locat)
    }
    else
    {
      if (!is(data, "data.frame"))
         stop("Given regression data must be data.frame or Spatial*DataFrame")
    }
  }
  else stop("No regression data frame is avaiable!")
  if (missing(dMat))
  {
      dMat <- gw.dist(dp.locat=dp.locat, p=p, theta=theta, longlat=longlat)
  }
  else
  {
    dim.dMat<-dim(dMat)
    if (!is.numeric(dMat)||!is(dMat, "matrix"))
      stop("Distance matrix(dMat) has to be specified correctly")
    if (dim.dMat[1]!=dp.n||dim.dMat[2]!=dp.n)
       stop("Dimensions of dMat are not correct")
  }
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data"), names(mf), 0)

    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    dp.n <- length(model.extract(mf, "response"))
    y <- model.extract(mf, "response")
    x <- model.matrix(mt, mf)
  var.n<-ncol(x)
  
  if (missing(bw))
      stop("Bandwidth must be given for non-adaptive weights")
  if (adaptive)
  {
    stopifnot(is.numeric(bw))
    stopifnot((bw > 0))
    stopifnot((bw <= dp.n))
  }
  else
  {
    stopifnot(is.numeric(bw))
    stopifnot((bw > min(dMat)))
  }
  bandwidth<-bw
 #####Calibrate the original GWR model
  betas <- matrix(nrow=dp.n, ncol=var.n)
  Var_betas<-matrix(nrow=nsims+1, ncol=var.n)
  for (i in 1:dp.n)
  {
    dist.vi<-dMat[,i]
    W.i<-gw.weight(dist.vi,bw,kernel,adaptive)
    gw.resi<-gw.reg(x,y,W.i,hatmatrix=F,i)
    betas[i,]<-gw.resi[[1]]
  }
  for (j in 1:var.n)
  {
    Var_betas[1,j]<-var(betas[,j])
  }
 ###################### Random sequence
  sq<-1:length(y)
  for(k in 1:nsims)
  {
    #randsq<-randomSQ(sq)
    #randx<-x[randsq,]
    #randy<-y[randsq]
    mcs <- sample(dp.n)
		dMat[mcs,]<-dMat[1:dp.n,]
		dMat[,mcs]<-dMat[,1:dp.n]
    betas <- matrix(nrow=dp.n, ncol=var.n)
    for (i in 1:dp.n)
    {
      dist.vi<-dMat[,i]
      W.i<-gw.weight(dist.vi,bw,kernel,adaptive)
      gw.resi<-gw.reg(x,y,W.i,hatmatrix=F,i)
      betas[i,]<-gw.resi[[1]]
    }
    for (j in 1:var.n)
    {
      Var_betas[k+1,j]<-var(betas[,j])
    }
  }
  ######Compute the p-values
  p.values<-numeric(var.n)
  for (j in 1:var.n)
  {
     var_betaj<-Var_betas[,j]
     indx<-sort(var_betaj, index.return=T)$ix
     p.values[j]=1-(which(indx==1))/(nsims+1)
  }
  pmat<-matrix(p.values, ncol=1)
  colnames(pmat)<-c("p-value")
  rownames(pmat)<-colnames(x)
  cat("\nTests based on the Monte Carlo significance test\n\n")
  printCoefmat(pmat)
	invisible(pmat)
}
