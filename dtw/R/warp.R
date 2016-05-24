###############################################################
#                                                             #
#   (c) Toni Giorgino <toni.giorgino,gmail.com>           #
#       Istituto di Ingegneria Biomedica (ISIB-CNR)                 #
#       Consiglio Nazionale delle Ricerche                           #
#       www.isib.cnr.it                                    #
#                                                             #
#   $Id: warp.R 311 2013-06-03 17:23:17Z tonig $
#                                                             #
###############################################################


## Apply the warping curve to a given timeseries. If the curve is
## multi-valued, average multiple points, with a warning. Argument
## "inverse=T" warps the inverse of the curve (ie a reference into a
## query).

## Direct:
##  all points in query have at least one image (NO?)
##  only one image, if asymmetric
##  if partial, image <= reference
##  x should be as long as the range in ix
##  there could be gaps in ix, to be interpolated,
##  depending on the step pattern

## for each point in the reference space, do an interpolated lookup
## into the y->x mapping

## Inverse: considering reference as domain
##  if partial, some trailing points may have no image
##  otherwise, all points have one or more images


## sortedXyData

warp <- function(d,index.reference=FALSE) {

  if(!is.dtw(d))
    stop("dtw object required");
  

  if(!index.reference) {
    ## warp QUERY into reference space
    iset<-d$index1;
    jset<-d$index2;
  } else {
    ## warp REFERENCE into query
    iset<-d$index2;
    jset<-d$index1;
  }
  jmax<-max(jset);

  ## rebuild  index, interpolating holes
  ii<-approx(x=jset,y=iset,1:jmax);
  return(ii$y);    
     
}


## > ss<-warp.dtw(al,query)
## Warning message:
## In approx(x = jset, y = iset, 1:jmax) : collapsing to unique 'x' values
## > plot(reference);lines(query)
## > plot(reference);lines(ss)
## > st<-warp.dtw(al,reference,inverse=T)
## Warning message:
## In approx(x = jset, y = iset, 1:jmax) : collapsing to unique 'x' values
## > Cairo()
## > plot(query);lines(st)
## > plot(query,type="l");points(st)
## > 
