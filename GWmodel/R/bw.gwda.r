###################################################
#Bandwidth selection using cross-validation
bw.gwda <- function(formula, data, COV.gw = T, prior.gw = T, mean.gw = T,
                 prior = NULL, wqda = F, kernel = "bisquare", adaptive
                 = FALSE, p = 2, theta = 0, longlat = F, dMat)
{
  #data must be given as training data
  if (is(data, "Spatial")) {
    p4s <- proj4string(data)
    dp.locat <- coordinates(data)
    data <- as(data, "data.frame")
  }
  else
    stop("Given training data must be a Spatial*DataFrame or data.frame object")
  ######## variables from formula
  vars <- all.vars(formula)
  grouping.nm <- vars[1]
  expl.vars <- vars[-1]
  m <- length(expl.vars)
  if (m < 2)
     stop("Two or more variables shoule be specfied for analysis")
  #x, y from training data
  res1 <- grouping.xy(data, grouping.nm, expl.vars)
  x<- res1$x
  grouping <- res1$y
  lev <- levels(grouping)
  dp.n<-nrow(data)
  if (missing(dMat))
  {
    dMat <- gw.dist(dp.locat=dp.locat, p=p, theta=theta, longlat=longlat) 
  }
  else
  {
    dim.dMat<-dim(dMat)
    if (dim.dMat[1]!=dp.n||dim.dMat[2]!=dp.n)
    stop ("Dimensions of dMat are not correct")
  }
  if(adaptive)
  {
    upper<-dp.n
    lower<-20
  }
  else
  {
    upper<-range(dMat)[2]
    lower<-upper/5000
  }
  bw<-NA
  if(wqda)
     bw <- optimize(wqda.cr,lower=lower,upper=upper,maximum=T,
                     x=x, grouping=grouping, dMat=dMat, COV.gw=COV.gw,
                 mean.gw=mean.gw, prior.gw=prior.gw, prior=prior,
                 kernel = kernel, adaptive =adaptive)
  else
     bw <- optimize(wlda.cr,lower=lower,upper=upper,maximum=T,
                     x=x, grouping=grouping, dMat=dMat, COV.gw=COV.gw,
                 mean.gw=mean.gw, prior.gw=prior.gw, prior=prior,
                 kernel = kernel, adaptive =adaptive)
  if(adaptive)
    bw <- round(bw$maximum)
  else
    bw <- bw$maximum
  bw
}

###Correct ratio
wqda.cr <- function(bw, x, grouping, dMat, COV.gw=T,
                 mean.gw=T, prior.gw=T, prior=NULL,
                 kernel = "bisquare", adaptive = FALSE)
{
  if(adaptive) bw <- round(bw)
  wt <- gw.weight(dMat,bw,kernel,adaptive)
  diag(wt) <- 0
  res.df <- try(wqda(x, grouping, x, wt, COV.gw,
                 mean.gw, prior.gw, prior))
  if(!inherits(res.df, "try-error"))
  {
     n.correct <- length(which(grouping == res.df[,"group.predicted"]))
     correct.ratio <- n.correct/nrow(x)
  }
  else
     correct.ratio <- 0
  if(adaptive)
      cat("Adaptive bandwidth:", bw, "Correctly predicted proportion:", correct.ratio, "\n")
    else
      cat("Fixed bandwidth:", bw, "Correctly predicted proportion:", correct.ratio, "\n")
  correct.ratio
}

wlda.cr <- function(bw, x, grouping, dMat, COV.gw=T,
                 mean.gw=T, prior.gw=T, prior=NULL,kernel = "bisquare", adaptive = FALSE)
{
  if(adaptive) bw <- round(bw)
  wt <- gw.weight(dMat,bw,kernel,adaptive)
  diag(wt) <- 0
  res.df <- try(wlda(x, grouping, x, wt, COV.gw,
                 mean.gw, prior.gw, prior))
  if(!inherits(res.df, "try-error"))
  {
     n.correct <- length(which(grouping == res.df[,"group.predicted"]))
     correct.ratio <- n.correct/nrow(x)
  }
  else
     correct.ratio <- 0
  if(adaptive)
      cat("Adaptive bandwidth:", bw, "Correctly predicted proportion:", correct.ratio, "\n")
    else
      cat("Fixed bandwidth:", bw, "Correctly predicted proportion:", correct.ratio, "\n")
  correct.ratio
}