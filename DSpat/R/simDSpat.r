simDSpat=function(study.area=owin(xrange=c(0,100),yrange=c(0,100)),covariates,
                  angle=90,nlines=10,spacing=10,width=1,int.formula=~1,
                  int.par=1,model="exp",cor.par=NULL,EN=1000,detfct=hndetfct,
                  det.formula=~1,det.par=log(width/3),showplot=FALSE,
                  showlines=FALSE,showpts=FALSE,pts=NULL,...)
#################################################################################
#  This is a wrapper function that calls all of the functions needed to simulate
#  and sample a point process over a defined study area with a specified covariate
#  structure.  In sequence it calls create.lines, lines_to_strips, simPts,
#  and sample.points.
#
#  Arguments:
#  study.area    -owin class defining area
#  covariates    -a matrix with columns x,y and any number of covariates
#                     x and y are the mid points of the grid cells; the order
#                     of the rows must match the formulation for function im
#  angle          -angle of rotation in degrees anticlockwise from x-axis
#  nlines         -number of lines
#  spacing        -spacing distance between centerlines
#  width          -full transect width
#  int.formula    -formula for deriving expected intensity from covariates
#  int.par        -parameters for intensity formula
#  model          -either "exp" or "gauss" for exponential or Gaussian correlation
#  cor.par        -parameters controlling clustering of points
#                     cor.par[1] sigma^2 cor.par[2]=alpha
#                     where cov(y1,y2)=sigma^2*exp(-d^p/alpha) and
#                     d is the distance between y1 and y2 and p=1 for exp and
#                     p=2 for gauss; if it is not specified then no additional
#                     clustering is included.
#  EN             -expected number of points
#  detfct         -detection function name
#  det.formula    -formula of covariates to use for scale of distance
#                     if det.formula=~-1, uses a strip transect
#  det.par        -parameters for the detection function
#  showplot       -if TRUE show plot of the simulated points
#  showlines      -if TRUE show lines and transects on the plot
#  showpts        -if TRUE show points on the plot
#  pts            -if not NULL use these points rather than generating new ones;
#                     this allows generation of a single set of points and
#                     evaluation of different sampling designs or intensity}
# ...             - parameters for plot
#
#
#  Value: a list with elements
#
#    lines        - lines dataframe with label,x0,y0,x1,y1,width where x0,y0 is beginning
#                   and x1,y1 is end of the line
#    observations - an observation dataframe
#
#  Jeff Laake
#  23 April 2008
#################################################################################
{
#  Create sample of lines
   xlines=create.lines(study.area=study.area,nlines=nlines,width=width,angle=angle)
#  Create strips and line psp
   ls=lines_to_strips(xlines,study.area)
#  Simulate point process across study area unless pts are supplied
   if(is.null(pts))
      obs=simPts(covariates=covariates,int.formula=int.formula,
                     int.par=int.par,EN=EN,showplot=showplot, showpts=showpts, ...)
   else
      obs=pts
#  Sample points via lines and detection
   observations=sample.points(ls$transects,xlines,obs,detfct=hndetfct,
                det.par=det.par,det.formula=det.formula,covariates=covariates)
#  If showplot=TRUE, show points and lines
   if(showplot & showlines)
   {
      plot(ls$lines,lty=2,add=TRUE)
      plot(owin(poly=ls$transects),add=TRUE)
   }
   return(list(lines=xlines,observations=observations,pts=obs))
}



