`lgcp.correction` <-
function(fit.ppm, fit.lgcp, reps=100, J.inv, silent=FALSE,
                              lines.psp)
################################################################################
# Calculate Overdispersion factor via Monte Carlo Integration
# 12/20/2007
#
# Arguments:
#  fit.ppm   -fitted model from ppm of spatstat
#  fit.lgcp  -fitted model from lgcp.estK
#  reps      -number of replicates for approximation
#  J.inv     -variance-covariance matrix from fitted ppm model
#  silent    -if FALSE, shows counter for replicates
#  lines.psp -line segment process
#
# Value:
#  J.inv.corr -Adjusted var-cov matrix
#  u          -score matrix
#
#
# Devin Johnson
################################################################################
{
  cat('Be patient, this might take a few minutes,... \n')
  sigma2 = fit.lgcp$par[1]
  alpha = fit.lgcp$par[2]
  winPoly = fit.ppm$Q$dummy$window$bdry
  ntran = length(winPoly)
  Qsumm = summary(fit.ppm$Q)
# grid.pts are only the quadrature dummy points
  grid.pts = fit.ppm$Q$dummy
  dummy.labels=grid.pts$label
  labels=lines.psp$label
  angles=lines.psp$angles
# true/false vector to select quad dummy points
  grid.ind = !is.data(fit.ppm$Q)
# covariate design matrix for dummy points; retrieval depends on whether it
# is glm or gam
  if(!is.null(fit.ppm$internal$glmfit$gcv.ubre)&!is.null(fit.ppm$internal$glmfit$smooth))
     grid.cov = cbind(x=grid.pts$x, y=grid.pts$y,
                predict(fit.ppm$internal$glmfit,type="lpmatrix")[grid.ind,])
  else
     grid.cov = cbind(x=grid.pts$x, y=grid.pts$y,
                model.matrix(fit.ppm$trend, fit.ppm$covariates[grid.ind,]))
# compute adjusted lambda for dummy points
  lambda.grid = fitted(fit.ppm)[grid.ind]/exp(sigma2/2)
# in grid.cov rotate x,y for each line and restore with new coordinates but same covariate values
# Use that to create a list of masks and images (one for each line) of the adjusted lambda values
# at the rotated (to vertical) dummy points in each line
  i=0
  lambdaLst=vector("list",length=length(labels))
  maskLst=vector("list",length=length(labels))
  for (label in labels)
  {
    i=i+1
    newdummy=rotate(ppp(grid.pts$x[dummy.labels==label],grid.pts$y[dummy.labels==label],
                           window=owin(poly=winPoly[i])),angle=angles[i])
    newdummy$x=round(newdummy$x,8)
    newdummy$y=round(newdummy$y,8)
    grid.cov[,"x"][dummy.labels==label]=newdummy$x
    grid.cov[,"y"][dummy.labels==label]=newdummy$y
    rotated.transect=rotate(owin(poly=winPoly[[i]]),angle=angles[i])
    maskLst[[i]]=as.mask(rotated.transect,xy=list(x=newdummy$x, y=newdummy$y))
    lambdaLst[[i]]=im.clipped(rev_val(newdummy$x,newdummy$y,
                           lambda.grid[dummy.labels==labels[i]]),maskLst[[i]])
  }
  u <- matrix(NA,reps,dim(grid.cov)[[2]]-2)
# create function to compute new random points for a transect
  new.pts=function(k, lambdaLst, z, labels, dummy.labels)
  {
      Zim=im.clipped(rev_val(grid.cov[,"x"][dummy.labels==labels[k]],
                             grid.cov[,"y"][dummy.labels==labels[k]],
                             exp(z[dummy.labels==labels[k]])),
                     maskLst[[k]]
                    )
      lamIm = lambdaLst[[k]]
      lamIm = eval.im(lamIm*Zim)
      rpoispp(lamIm)
  }
# create function to get covariates for new points
  get.cov=function(k, XLst, maskLst, grid.cov)
  {
     grid.snp = as.data.frame(nearest.raster.point(XLst[[k]]$x, XLst[[k]]$y, maskLst[[k]],indices=FALSE))
     return(merge(grid.snp,grid.cov))
  }
  cat('Beginning Monte Carlo reps, patience please,... \n')
# Loop over each simulation rep
  for(i in 1:reps)
  {
#   Generate a random field for the original coordinate system
    #z = GaussRF(x=cbind(grid.pts$x,grid.pts$y), model="exp", grid=FALSE, param=c(0,sigma2,0,alpha))
    z = RFsimulate(model=RMexp(var=sigma2, scale=alpha), x=grid.pts$x, y=grid.pts$y, grid=FALSE)@data[,1]
#   For each line, generate a random set of points using the lambdaLst and correlated field
    XLst =  lapply(c(1:ntran),new.pts,lambdaLst=lambdaLst,z=z,labels=labels,dummy.labels=dummy.labels)
#   For the generated set of points, snap them to the mask to get the covariate values
    sim.cov=do.call('rbind',lapply(c(1:ntran),get.cov, XLst=XLst,
                             maskLst=maskLst, grid.cov=grid.cov))
#    Compute score statistic for the ith replicate
     u[i,] <- apply(sim.cov[,-c(1:2),drop=FALSE],2,sum)
     if(!silent) cat(paste('Rep',i,'of',reps,'completed.\n',collapse=' '))
  }
# return the adjusted v-c matrix
  return(list(J.inv.corr = J.inv%*%var(u)%*%J.inv, u=u))
}

