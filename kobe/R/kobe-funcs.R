##################################################################################
##                                                                              ##
##  A range of methods for                                                      ##
##     file IO,                                                                 ##
##     summarising time series and                                              ##    
##     plotting densities on kobe plots                                         ##                        
##                                                                              ##
##                                                                              ##
##################################################################################

#' @importFrom plyr ddply ldply .
#' @importFrom MASS kde2d bandwidth.nrd 
#' @importFrom emdbook HPDregionplot
#' @importFrom coda mcmc
#' @import methods


utils::globalVariables(c("xFac","yFac","freq","count"))

ac=as.character

getExt=function(file)
  tolower(substr(file,max(gregexpr("\\.", file)[[1]])+1,nchar(file)))

#' kobeFreq
#' @description 
#' Calculates the frequency of an obervation in a 2D cell
#'           
#' @aliases kobeFreq
#' 
#' @param x a vector holding a time series
#' @param y a vector holding a time series
#' @param x.n a numeric vector giving number of bins
#' @param y.n a numeric vector giving number of bins
#' @param na.rm logical; if true, any NA and NaN's are removed from x before calculations
#'
#'
#' @return a \code{data.frame} frequency by bin
#' @export
#' @docType methods
#' @rdname kobeFreq
#' 
#' @examples
#' \dontrun{
#'    y=rnorm(20)
#'    x=rnorm(20)
#'    kobeFreq(x,y)}
kobeFreq=function(x,y,x.n=11,y.n=x.n,na.rm=FALSE){
  
  if (na.rm){
    .na=is.na(x)|is.na(y)|is.nan(x)|is.nan(y)
    x=x[.na]
    y=y[.na]}
  
  df=data.frame(x=x,y=y)
  df=data.frame(df,xFac  =cut(df$x,  seq(min(df$x),  max(df$x),  length.out=x.n)),
                yFac=cut(df$y,seq(min(df$y),max(df$y),length.out=y.n)))
  
  c.=ddply(data.frame(df,count=1),.(xFac,yFac), function(x) count(x$count))[,c("xFac","yFac","freq")]
  
  p.=merge(df,c.,by=c("xFac","yFac"))[,c("x","y","freq","xFac","yFac")]
  p.=ddply(p.,.(xFac,yFac),with, sum(freq))
  names(p.)=c("x","y","freq")
  
  return(p.)}

## calculates density of points
#' kobeDens
#' @description 
#' Calculates the Densities of obervation in a 2D cell using Two-dimensional
#'  kernel density estimation with an axis-aligned bivariate normal kernel, 
#'  evaluated on a square grid.
#'           
#' @aliases 
#' kobeDens
#' 
#' @param x a vector 
#' @param y a vector
#' @param n Number of grid points in each direction. Can be scalar or a length-2 integer vector.
#' @param h  vector of bandwidths for x and y directions. Defaults to normal reference bandwidth 
#' (see bandwidth.nrd). A scalar value will be taken to apply to both directions.
#' @param lims The limits of the rectangle covered by the grid as c(xl, xu, yl, yu).
#' @param na.rm logical; if true, any NA and NaN's are removed from x before calculations
#'
#' @return a \code{data.frame} with three variables
#' \code{x, y} coordinates of the grid points, vectors of length n.
#' \code{z} An n[1] by n[2] matrix of the estimated density: rows correspond to the value of 
#' x, columns to the value of y.
#' 
#' @export
#' @docType methods
#' @rdname kobeDens
#' 
#' @examples
#' \dontrun{
#'    y=rnorm(20)
#'    x  =rnorm(20)
#'    kobeDens(x,y)}
kobeDens=function(x,y,h=c(bandwidth.nrd(x),bandwidth.nrd(y)),n=11,lims=c(range(x),range(y)),na.rm=FALSE){
  
  if (na.rm){
    .na=is.na(x)|is.na(y)|is.nan(x)|is.nan(y)
    x=x[.na]
    y=y[.na]}
  
  dat=data.frame(x=x,y=y,n=n)
  f1 =with(dat, kde2d(x,y,h=h,n=n,lims=lims)) 
  f2 =data.frame(expand.grid(x=f1$x, y=f1$y), z=as.vector(f1$z))
  
  return(f2)}

## calculates probabilities
#' kobeProb
#' @description 
#' Calculates the probability of an obervations occurring in a 2D cell using HPDregionplot 
#' Given a sample calculates the  bivariate region of highest marginal posterior density 
#' for two variables, using kde2d from MASS to calculate a bivariate density.
#'            
#' @aliases 
#' kobeProb
#' 
#' @param x a vector
#' @param y a vector
#' @param prob probability levels
#' @param n number of points at which to evaluate the density grid
#' @param h bandwidth of 2D kernel smoother (previous default value was c(1,1), which worked 
#' poorly with some plots with very small scales; if not specified, defaults to values in kde2d)
#' @param lims limits, specified as (x.lower,x.upper,y.lower,y.upper) (passed to kde2d)
#' @param na.rm logical; if true, any NA and NaN's are removed from x before calculations
#'
#' @return a \code{data.frame} with three variables
#' \code{x, y} coordinates of the grid points, vectors of length n.
#' \code{level} contours corresponding to \code{prob}
#' 
#' 
#' @export
#' @docType methods
#' @rdname kobeProb
#' 
#' @examples
#' \dontrun{
#'    y=rnorm(20)
#'    x  =rnorm(20)
#'    kobeProb(x,y)}
kobeProb=function(x,y,prob=c(0.5, 0.75,0.95),n=21,h=c(bandwidth.nrd(x),bandwidth.nrd(y)),lims=NULL,na.rm=FALSE){
  
  if (na.rm){
    .na=is.na(x)|is.na(y)|is.nan(x)|is.nan(y)
    x=x[.na]
    y=y[.na]}
  
  tmp=HPDregionplot(mcmc(data.frame(x,y)),prob=prob,h=h)
  
  prb=ldply(tmp, function(dat) data.frame(level=dat$level,x=dat$x, y=dat$y))
  
  return(prb)}

#Utility methods for summarising time series to create performance measures

#' iav
#' @description 
#' Calculates the inter-annual variation in a time series, i.e. \code{(x[t+1]-x[t])/x[t]}
#' Used to show how variable a quantity like yield 
#' is under different management strategies within a Management Strategy Evaluation.
#'      
#' @aliases  iav
#' 
#' @param x a vector holding a time series
#' @return a \code{vector} with the inter-annual variation each time step
#' 
#' 
#' @export
#' @docType methods
#' @rdname iav
#' 
#' @examples
#' \dontrun{
#'    x=rnorm(2)
#'    iav(x)
#'    ## inter-annual average variation
#'    mean(iav(x),na.rm=T)}
iav=function(x) if (length(x)==1) return(NA) else c(NA,(x[-1]-x[-length(x)])/x[-length(x)])


#' incr
#' @description 
#' Is a quantity increasing from 1 time step to another \code{(x[t+1]-x[t])>0} 
#'      
#' @aliases 
#' incr
#' 
#' @param x a vector holding a time series
#' @return a \code{logical} indicating an increase
#' @export
#' @docType methods
#' @rdname incr
#' 
#' @examples
#' \dontrun{
#'    x=rnorm(2)
#'    incr(x)}
incr=function(x) c(NA,(x[-1]-x[-length(x)]>0))


#' recovered
#' @description 
#' Has the stock recovered yet? i.e. stock>=1 and harvest<=1 in the current or an earlier time step.
#' In other words has it been in the green Kobe quadrant.
#'      
#' @aliases 
#' recovered
#' 
#' @param stock a vector holding a time series
#' @param harvest a vector holding a time series
#' @return a \code{logical} indicating a recovered stock
#' 
#' @export
#' @docType methods
#' @rdname recovered
#' 
#' @examples
#' \dontrun{
#'    harvest=rlnorm(20)
#'    stock  =rlnorm(20)
#'    recovered(stock,harvest)}
recovered=function(stock,harvest) {
  
  x =pmin(floor(stock),1)
  y =1-pmin(floor(harvest),1)
  
  as.logical(cumsum(y*x))}

#' dRate
#' @description 
#' Discount rate reflects the risk (the higher the risk the higher the discount rate). Is used 
#' to discount all forecast future cash flows to calculate a present value:
#'      
#' @aliases 
#' dRate
#' 
#' @param x a vector holding a time series
#' @param r a numeric discount rate
#' @param wtAv a logical, FALSE by default, if TRUE then uses discount rate to calculate a weighted average
#' @return net present value
#' 
#' @export
#' @docType methods
#' @rdname dRate
#' 
#' @examples
#' \dontrun{
#'    x=rnorm(20)
#'    dRate(x,0.05)}
dRate=function(x,r,wtAv=FALSE) {
  if (wtAv) return( sum(x/(1+r)^(0:(length(x)-1)))/sum(1/(1+r)^(0:(length(x)-1))))
  return(sum(x/(1+r)^(0:(length(x)-1))))
}

#' incr
#' @description 
#' Break point for segmented regression stock recruitment relationship 
#' for a given steepness and virgin biomass
#'      
#' @aliases 
#' hinge
#' 
#' @param s steepness
#' @param v virgin biomass
#' @param rec a vector holding constant recruitment level
#' @return a \code{numeric} giving break point
#' @export
#' @docType methods
#' @rdname  hinge
#' 
#' @examples
#' \dontrun{
#'    hinge(.7,1000,2000)}
hinge=function(s,v,rec){
  #slope
  x=c(0,v*0.2)
  y=c(0,rec*s)
  
  #y=a+b*x
  #y1=a+b*x1
  #y2=a+b*x2
  #y1-y2=b*(x1-x2)
  
  b=(y[1]-y[2])/(x[1]-x[2])
  a=b*x[1]-y[1]
  
  (rec-a)/b}
