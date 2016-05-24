quadscheme.lt <-
function(study.area,observations,lines,width=NULL,epsvu=c(1,.01),show.warnings=FALSE)
################################################################################
# Creates a quadrature for spatstat from a study area, observations, and lines
# ppp process.
#
# Arguments:
#
#  study.area   - owin class giving the boundaries of the study area
#
#  observations - data frame of observations with the following structure
#                    label - label linking it to a unique line
#                    x     - x coordinate
#                    y     - y coordinate
#                    ...   - any number of covariates
#
#  lines        - data frame of lines with the following structure
#                    label - unique label
#                    x0    - x coordinate of beginning of line
#                    y0    - y coordinate of beginning of line
#                    x1    - x coordinate of end of line
#                    y1    - y coordinate of end of line
#                    width - full width of transect around line (optional)
#                    ...   - any number of covariates
#
#  width        - if no width field is given in lines then specified here
#                 as a constant width for all lines
#
#  epsvu        - pixel dimensions height = epsvu[1] in v direction and
#                 width = epsvu[2] in u direction (these are used
#                 after the line is rotated vertically.
#
#  show.warnings - if TRUE, warnings from quadrature construction will be shown
#
#
# Value: a list with the following elements
#   Q            - quadscheme
#   transects    - transect polygons
#   lines.psp    - line segment process
#
# Devin Johnson & Jeff Laake
# 4 April 2008
################################################################################
{
  method='gridweights'
  number.lines = dim(lines)[1]
#
# From the lines dataframe, construct a set of lines.psp, transect.polygons and
# rotation angles.
#
  transect.list=lines_to_strips(lines=lines, study.area=study.area, width=width)
  lines.psp=transect.list$lines
  angles=lines.psp$angles
  transects=transect.list$transects
  full.transects=transect.list$full.transects
  lengths=lengths.psp(lines.psp)
#
# Loop over each line
#
  dummyLst=vector("list",number.lines)
  wLst=vector("list",number.lines)
  wData=NULL
  wDummy=NULL
  DummyLabels=NULL
  nx=round(max(lines.psp$width)/epsvu[2]+.5)
  for (i in 1:number.lines)
  {
#   Get line label and observations
    line.label=lines$label[i]
    obs=observations[observations$label==line.label,]
#   Create observation point process and rotate it and the transect
    obs.ppp.rotate=rotate(ppp(obs$x,obs$y,window=owin(poly=transects[i])),angle=angles[i])
    win = rotate(owin(poly=transects[i]),angle=angles[i])
    win.full = rotate(owin(poly=full.transects[i]),angle=angles[i])
#   Compute number of pixels in x and y direction
    ny=round(lengths[i]/epsvu[1]+.5)
#   Compute dummy data point process
    dummy=data.frame(gridcenters(window=win.full, nx=nx, ny=ny))
    instrip=inside.owin(dummy$x, dummy$y, win)
    dummyLst[[i]]=ppp(dummy$x[instrip],dummy$y[instrip],window=win)
#   Concatenate data and dummy points
    X = superimpose(obs.ppp.rotate,dummyLst[[i]],W=win)
#   Calculate weights and then partition into data and dummy wts
    if(show.warnings)
       weightLst = do.call(method,list(X=X,ntile=c(nx,ny),window=win,exact=FALSE))
    else
       weightLst = suppressWarnings(do.call(method,list(X=X,ntile=c(nx,ny),
                                           window=win,exact=FALSE)))
    if(dim(obs)[1]>0)
    {
       wData=c(wData,weightLst[1:dim(obs)[1]])
       wDummy=c(wDummy,weightLst[(dim(obs)[1]+1):length(weightLst)])
    }
    else
       wDummy=c(wDummy,weightLst)
    DummyLabels=c(DummyLabels,rep(line.label,sum(as.numeric(instrip))))
#   Rotate dummy points back to original strip
    dummyLst[[i]]=rotate(dummyLst[[i]],angle=-angles[i])
   }
#  Create observation ppp and dummy ppp;gridwe a label is added to obs.ppp and dummy.ppp
   obs.ppp=ppp(observations$x,observations$y,window=owin(poly=transects))
   obs.ppp$label=observations$label
   dummy.ppp = do.call(superimpose, append(dummyLst, list(W=owin(poly=transects))))
#   dummy.ppp = superimpose(dummyLst,W=owin(poly=transects))
   dummy.ppp$label=DummyLabels
#  Order wts and then create quadrature
   wts=c(wData,wDummy)
   Q = list(data=obs.ppp, dummy=dummy.ppp, w=wts, param=NULL)
   class(Q) = "quad"
   return(list(Q=Q,transects=transects,lines.psp=lines.psp))
}

