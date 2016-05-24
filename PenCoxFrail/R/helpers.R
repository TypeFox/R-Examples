## help functions

bs.design<-function(x, xl, xr, spline.degree, nbasis, comp = NULL)
{
  
  ## generate a B-Spline-Matrix with equidistant knots (code by Thomas Kneib & Andreas Groll):
  ## x are the positions where spline to be evaluated
  ## xl, xr intervall boundaries where spline functions are relevant
  xmin<-xl-(xr-xl)/100
  xmax<-xr+(xr-xl)/100
  dx<-(xmax-xmin)/(nbasis-3)
  knots<-seq(xmin-spline.degree*dx,xmax+spline.degree*dx,by=dx)
  B<-splines::spline.des(knots,x,spline.degree+1,outer.ok = TRUE)$design
  if(is.null(comp))
  {
    return(B)
  }else{
    return(B[,comp])
  }
}

########

mirror <- function(x,low.tri)
{
  x[low.tri] <- t(x)[low.tri]
  return(x)
}

#############

penal.fct <- function(x,K)
{
  sqrt(t(x) %*% (K %*% x))
}

#############

factor.test <- function(x)
{
  substr(x,1,9)=="as.factor"
}

#############

penal.fct.inv <- function(x,K,c.app)
{
  (as.numeric(t(x) %*% (K %*% x))+c.app)^(-0.5)
}

#############

eucl.norm <- function(x)
{
  sqrt(sum(x^2))
}

#############

count.levels <- function(x)
{
  nlevels(x)
}

#############

center.fct <- function(fit.vec,mean.vec,sd.vec,nbasis,m, standardize = TRUE, center = TRUE)
{
  if(!center)
    mean.vec <- rep(0,length(mean.vec))
  if(!standardize)
    sd.vec <- rep(1,length(sd.vec))
  
    Delta.base.retrans <- c(1,-mean.vec/sd.vec) %*% t(matrix(fit.vec[1:(nbasis*(m+1))],nbasis,m+1))
  return(Delta.base.retrans)
}

############

int.approx <- function(z,time.grid,B,nbasis,alpha)
{
  index <- time.grid<z[1]
  diff(time.grid[1:2]) * sum(exp(B[index,] %*% (alpha * rep(z[2:length(z)],each=nbasis))))
}

