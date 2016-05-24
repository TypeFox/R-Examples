integrate.intensity=function(x, dimyx=NULL, eps=NULL, se=FALSE, od=FALSE,
                            reps=100, silent=FALSE, J.inv=NULL, showplot=TRUE)
################################################################################
# Calculate b(theta) vector for Expected count delta method calculation
#
# Arguments:
#
#   x         - fitted dspat object
#   dimyx     - number of y,x pixels
#   eps       - height and width of pixels
#   se        - if FALSE only returns mu.B and lambda
#   od        - if FALSE, does not attempt over-dispersion correction
#   reps      - number of reps for MC integration
#   silent    - if FALSE, show progress of MC reps
#   J.inv     - var-cov matrix of fit
#   showplot  - if TRUE show plot of Poisson vs empirical & lgcp fitted K functions
#
# Value: list with the following elements:
#
#   Abundance     - estimate of expected total abundance
#   distribution  - dataframe containing N (predicted number of points in the cell),
#                    x,y (x,y coordinates of cell) and covariates used in the model}
#   precision     - List containing se, lcl.95, ucl.95, J.inv, and b.vec
#   precision.od  - For over-dispersion estimate a list containing se, lcl.95, ucl.95, J.inv, and b.vec
#   lambda        - image of predicted intensity
#   W             - mask of study area
#
# Devin Johnson & Jeff Laake
# 9 April 2008
################################################################################
{
log.normal.ci=function(est,se)
{
  ll = exp(log(est) - 1.96*se)
  ul = exp(log(est) + 1.96*se)
  se = est*se
  return(list(se=se,lcl.95=ll,ucl.95=ul))
}
# Check to make sure x is of class dspat
  if(class(x)[1]!="dspat") stop("\n Argument x must be of class dspat\n")
#
# Create mask for study area
#
  if(is.null(dimyx))
  {
     if(is.null(eps))
     {
        B.mask=as.mask(x$study.area,xy=list(x=x$covariate.im[[1]]$xcol,y=x$covariate.im[[1]]$yrow))
        Area=x$covariate.im[[1]]$xstep*x$covariate.im[[1]]$ystep
     }
     else
     {
        B.mask=as.mask(x$study.area,eps=eps)
        Area=eps^2
     }
  }
  else
  {
     B.mask=as.mask(x$study.area,dimyx=dimyx)
     Area=B.mask$xstep*B.mask$ystep
  }
#
# if distance included in model, then append an image with 0 distance everywhere
#
  varnames=names(coef(x$model))
  if(!is.null(x$covariate.im) & length(x$covariate.im)>0)
  {
     anyimage=x$covariate.im[[1]]
     covariates=x$covariate.im
  }
  else
  {
     anyimage=NULL
     covariates=vector("list",0)
  }
  if (length(grep("distance",varnames)!=0))
  {
     if(is.null(anyimage))
     {
        x.pg=unique(B.mask$xcol)
        y.pg=unique(B.mask$yrow)
        covariates=list(distance=im(matrix(0,ncol=length(x.pg),nrow=length(y.pg)),xcol=x.pg,yrow=y.pg))
     }
     else
        covariates=append(covariates, list(distance=eval.im(0*anyimage)))
  }
#
# compute intensity over grid (lambda)
#
  x$model$Q$data$marks <- x$model$Q$dummy$marks <- NULL
  locs = as.ppp(expand.grid(x=B.mask$xcol,y=B.mask$yrow)[,2:1],W=boundingbox(B.mask))
  locs=locs[B.mask]
  lambda <- predict(x$model,locations=locs,covariates=covariates,type="lambda")
# create dataframe of covariate values for grid
  nms <- names(covariates)
  covData <-   lapply(nms,
                function(nm, covariates, locs)
                {
                  covariates[[nm]][locs]
                }, covariates=covariates, locs=locs
              )
  names(covData) <- nms
  covData = as.data.frame(covData)
#
#  include x and y
#
  if(dim(covData)[1]==0)
     covData=data.frame(x=locs$x,y=locs$y)
  else
  {
     covData$x=locs$x
     covData$y=locs$y
  }
  distribution=cbind(N=lambda*Area,covData)
  if("distance" %in% nms) distribution$distance=NULL
# compute abundance over grid
  mu.B = sum(lambda*Area)
  if(!se) return(list(abundance=mu.B,distribution=distribution,
            lambda=im.clipped(lambda,B.mask), W=B.mask))
#
# If se=TRUE compute standard error and confidence interval
#
  if(x$use.gam)
       mm = predict(x$model$internal$glmfit,type="lpmatrix",newdata=covData)
  else
       mm <- model.matrix(x$model$trend, covData)
  mmNms <- dimnames(mm)[[2]]
  b.vec <- crossprod(mm,lambda*Area)/mu.B
  if(missing(J.inv)) J.inv <- vcov(x)
  se.ln.mu.B = sqrt(t(b.vec)%*%J.inv%*%b.vec)
  abundance.intervals=log.normal.ci(mu.B,se.ln.mu.B)
  precision=log.normal.ci(mu.B,se.ln.mu.B)
  precision$J.inv=J.inv
  precision$b.vec=b.vec

#
# If suppose to do over-dispersion correction factor continue
#
  if(od)
  {
    fit.ppm=x$model
    lam.obs <- as.vector(fitted(fit.ppm)[is.data(fit.ppm$Q)])
# Calculate empirical K function for IPP
# max distance, 'r', should be half width of a transect strip
    rmax <- min(x$lines.psp$width/2)
    K.hat <- Kinhom(X=x$model$Q$data, r=seq(0, rmax, rmax/100),
                 lambda=lam.obs, correction="translate" )
# Estimate LGCP spatial parameters vis 2-stage estimation
    fit.lgcp <- lgcp.estK( X=K.hat,startpar=c(sigma2=1,alpha=1),lambda=NULL,
                            q= 1/4, p = 2,rmin = 0, rmax = rmax )
# Show fit. Red=empirical, Black=LGCP fit, Green=IPP theoretical
# If black line is above green, there is some evidence of overdispersion.
    if(showplot) plot(fit.lgcp)
# Perform Monte Carlo calculations to adjust J.inv for overdispersion
# reps = number of independent replications for MC estimate
# Use silent=TRUE to avoid iteration updates
# !!! This is computationally intensive and will take a few minutes.
    od.corr <- lgcp.correction(fit.ppm=fit.ppm, fit.lgcp=fit.lgcp, reps=reps,
                J.inv=J.inv, lines.psp=x$lines.psp, silent=silent)
# Adjust initial estimates with MC overdispersion correction
    J.alt <- as.matrix(od.corr$J.inv.corr)
    b.vec <- as.matrix(b.vec)
    se.ln.mu.B.corr <- sqrt(t(b.vec)%*%J.alt%*%b.vec)
    precision.od=log.normal.ci(mu.B,se.ln.mu.B.corr)
    precision.od$J.inv=J.alt
    precision.od$b.vec=b.vec
    return(list( abundance=mu.B,distribution=distribution,precision=precision,
            precision.od=precision.od,
            lambda=im.clipped(lambda,B.mask), W=B.mask))
  }
  else
   return(list( abundance=mu.B,distribution=distribution,precision=precision,
            lambda=im.clipped(lambda,B.mask), W=B.mask))
}
