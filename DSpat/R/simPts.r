simPts=function(covariates,int.formula=~1, int.par=c(1), EN=100,
                model=NULL, showplot=FALSE, showpts=FALSE, ...)
################################################################################
#
# Simulates point process (homogeneous Poisson, inhomogeneous Poisson and
# both with additional clustering of points) over a rectangular grid.
# The matrix covariates defines the window and the covariate values. The formula
# for the covariates and int.par determine the intensity surface except that the
# intercept is set to provide a certain expected number of points (EN).
#
# Arguments:
#
#   covariates   - a matrix with columns x,y and any number of covariates
#                  x and y are the mid points of the grid cells; the order
#                  of the rows must match the formulation for function im
#
#   int.formula  - formula for deriving expected intensity from covariates
#
#   int.par      - parameters for intensity formula
#
#   EN           - expecte number of points
#
#   model        - a model for RFsimulate()
#
#   cor.par      - parameters controlling clustering of points
#                    cor.par[1] sigma^2 cor.par[2]=alpha
#                  where cov(y1,y2)=sigma^2*exp(-d^p/alpha) and
#                  d is the distance between y1 and y2 and p=1 for exp and
#                  p=2 for gauss
#
#  showplot      - if true plot intensity and point process
#  showpts       -if TRUE show points on the plot
#   ...          - parameters passed to plot
#
# Value: ppp object of point locations
#
#  Author: Devin Johnson; Jeff Laake
#  10 April 2008
#############################################################################
{
# Get x and y covariates
  x=covariates[,"x"]
  y=covariates[,"y"]
# Create a design matrix and compute the intensity image which is
# adjusted so the expected number of points is EN
  design.matrix=model.matrix(int.formula,data=as.data.frame(covariates))
  if(dim(design.matrix)[2]!=length(int.par))
     stop("\nMismatch between formula and length of int.par")
  I.x = as.vector(exp(design.matrix%*%int.par))
  mu = log(EN)-log(sum(I.x))
  lamIm = im(exp(mu)*I.x,xcol=unique(x),yrow=unique(y))
# Add clustering to the intensity surface using the specified exp or gauss model
# and the sigma^2 and alpha parameters.Skip this if cor.par has not been assigned.
  if(!is.null(model))
  {
     #if(!model %in% c("exp","gauss"))stop("\n model must be either exp or gauss")
     z=RFsimulate(model=model, x=x, y=y, grid=FALSE)@data[,1]
     Zim = im(exp(z),xcol=unique(x),yrow=unique(y))
     lamIm = eval.im(lamIm*Zim)
  }
# Simulate and return a point process with the specified intensity image
  ppp=rpoispp(lamIm)
  if(showplot)
  {
     plot(lamIm,main="True Intensity",...)
     if(showpts) plot(ppp,add=TRUE)
  }
  return(ppp)
}
