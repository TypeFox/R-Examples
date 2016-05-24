pearson0moments <- function(mean,sd) {
  c(mean=mean,variance=sd^2,skewness=0,kurtosis=3)
}

pearsonImoments <- function(a,b,location,scale) {
  c(mean=scale*a/(a+b)+location,variance=scale^2*a*b/((a+b)^2*(a+b+1)),
    skewness=sign(scale)*2*(b-a)*sqrt(a+b+1)/((a+b+2)*sqrt(a*b)),
    kurtosis=3+6*(a^3-a^2*(2*b-1)+b^2*(b+1)-2*a*b*(b+2))/
             (a*b*(a+b+2)*(a+b+3)))
}

pearsonIImoments <- function(a,location,scale) {
  c(mean=scale/2+location,variance=scale^2*a^2/((2*a)^2*(2*a+1)),
    skewness=0,
    kurtosis=3+6*(-2*(a+1))/
             ((2*a+2)*(2*a+3)))
}

pearsonIIImoments <- function(shape,location,scale) {
  c(mean=scale*shape+location,variance=shape*scale^2,
    skewness=sign(scale)*2/sqrt(shape),kurtosis=3+6/shape)
}

pearsonIVmoments <- function(m,nu,location,scale,params) {
  if (!missing(params)) { m <- params[[1]]; nu <- params[[2]]; 
                          location <- params[[3]]; scale <- params[[4]] }
  r <- 2*(m-1)
  if (m>1) mmm <- location-scale*nu/r else mmm <- NaN                           # improvable?... +/- Inf ...
  if (m>3/2) vvv <- scale^2/(r^2*(r-1))*(r^2+nu^2) else vvv <- Inf
  if (m>2) sss <- -4*nu/(r-2)*sqrt((r-1)/(r^2+nu^2)) else sss <- NaN
  if (m>5/2) kkk <- 3*(r-1)*((r+6)*(r^2+nu^2)-8*r^2)/
                      ((r-2)*(r-3)*(r^2+nu^2)) else kkk <- NaN
  c(mean=mmm,variance=vvv,skewness=sss,kurtosis=kkk)
}

pearsonVmoments <- function(shape,location,scale) {                             
  c(mean=ifelse(shape>1,location+scale/(shape-1),NaN),
    variance=ifelse(shape>2,scale^2/((shape-1)^2*(shape-2)),NaN),
    skewness=ifelse(shape>3,4*sign(scale)*sqrt(shape-2)/(shape-3),NaN),
    kurtosis=ifelse(shape>4,3+(30*shape-66)/((shape-3)*(shape-4)),NaN))
}

pearsonVImoments <- function(a,b,location,scale) {
  mmm <- scale*a/(b-1)
  vvn <- scale^2*(a+1)*a/((b-1)*(b-2))
  vvv <- vvn-mmm^2
  ssn <- scale^3*(a+2)*(a+1)*a/((b-1)*(b-2)*(b-3))
  sss <- ssn-3*mmm*vvn+2*mmm^3
  kkn <- scale^4*(a+3)*(a+2)*(a+1)*a/((b-1)*(b-2)*(b-3)*(b-4))
  kkk <- kkn - 4*mmm*ssn+6*mmm^2*vvn-3*mmm^4
  c(mean=ifelse(b>1,mmm+location,NaN),
    variance=ifelse(b>2,scale^2*a*(a+b-1)/((b-2)*(b-1)^2),NaN),
    skewness=ifelse(b>3,sss/(vvv^(3/2)),NaN),
    kurtosis=ifelse(b>4,kkk/(vvv^2),NaN))
}

pearsonVIImoments <- function(df,location,scale) {
  c(mean=ifelse(df>1,location,NaN),
    variance=ifelse(df>2,scale^2*df/(df-2),Inf),
    skewness=ifelse(df>3,0,NaN),
    kurtosis=ifelse(df>4,3+6/(df-4),NaN))
}
