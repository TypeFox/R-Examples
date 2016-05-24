dspat=function(int.formula=~1,det.formula=~1,study.area,
               obs,lines,covariates,epsvu=c(1,.01), width=NULL,
               use.gam=FALSE, show.warnings=FALSE, nclass=NULL)
##################################################################################
# Creates a dspat object by fitting model represented by formula to observations
# along line transects in a study area with covariates defined for a grid over the
# entire study area.
#
# Arguments:
#
#   formula      - formula for model
#   study.area   - owin class specifying study area
#   obs          - dataframe of observations
#   lines        - dataframe of lines
#   covariates   - dataframe of covariates with following structure
#                     x   - x coordinate of midpoint of grid cell
#                     y   - y coordinate of midpoint of grid cell
#                     ... - any number of covariate
#                     the data are ordered by column from left to right and
#                     from bottom to top such that y changes first from smallest
#                     to largest. Below are matrices showing y,x and their order
#                     3,1 3,2 3,3   3 6 9
#                     2,1 2,2 2,3   2 5 8
#                     1,1 1,2 1,3   1 4 7
#   epsvu        - vector of pixel sizes (v,u) or (height, width) after rotation
#   width        - full transect width; only needed if it is not specified in lines
#   use.gam      - if TRUE uses gam instead of glm for fitting; if formula contains
#                   s() use.gam will be set TRUE by default
#   show.warnings - if TRUE, show the warnings created in building the quadrature
#   nclass        - number of distance classes for expected/observed counts
#
# Value:  list of class "dspat" with elements
#
#   model        - output object from ppm
#   lines.psp    - psp line segment process for center lines
#   transects    - list of dataframes specifying polygonal transects
#   covariate.im - list of covariate images (class im)
#   study.area   - owin class of study area
#   use.gam      - TRUE/FALSE whether gam was used
#
# Jeff Laake
# 9 April 2008
##################################################################################
{
# adjust epsvu[2] such that width is a multiple
  if(is.null(width))
  {
     if("width" %in% colnames(lines))
        width=max(lines[,"width"])
     else
        stop("Missing either width specification or width column in lines matrix")
  }
    else
  {
       lines = cbind(lines,width=rep(width,dim(lines)[1]))
  }
  epsvu[2]=width/2/ceiling(width/2/epsvu[2])
# Create lt quadscheme and extract components
Q.list <- quadscheme.lt(study.area,obs,lines,epsvu=epsvu,
                         width=NULL,show.warnings=show.warnings)
Q.lt=Q.list$Q
obs.ppp=Q.list$Q$data
lines.psp=Q.list$lines.psp
transects=Q.list$transects
# Turn covariates into spatstat images; handle HPP with no distance
# effect (strip transects) by creating constant covariate
if(int.formula==~1) int.formula=~-1+constant
if(int.formula==~-1+constant)
   covariates$constant=1
# Create covariate dataframe with all covariates used in formula at
# data (observations) and dummy(quadrature) points.
LTdf = LTDataFrame(study.area, lines, lines.psp, int.formula,
                   det.formula, covariates, Q.lt)
# Set value of use.gam to TRUE if formula contains a smooth "s("
if(length(grep("s\\(",as.character(LTdf$formula)))!=0) use.gam=TRUE
# Fit Model
fit.ppm = ppm(Q.lt,LTdf$formula,covariates = LTdf$cov.df, use.gam=use.gam)
x=list(model=fit.ppm,lines.psp=lines.psp,transects=transects,
      covariate.im=LTdf$covariate.im,study.area=study.area,use.gam=use.gam)
# Compute observed and expected counts by transect
class(x)="dspat"
x=c(x,transect.intensity(x, epsvu=epsvu, obs.ppp=obs.ppp,
        covariates=LTdf$cov.df,nclass=nclass,width=width))
# Return object with class "dspat"
class(x)="dspat"
return(x)
}

