sample.points<-function(transects,lines,points.ppp,detfct=NULL,det.par=NULL,
                           det.formula=~1,covariates=NULL)
################################################################################
# Sample points within each transect and filter with specified detection
# function.
#
# Arguments:
#
#  transects  - list of transect polygons
#  lines      - dataframe of lines
#  points.ppp - point process
#  detfct     - detection function name
#  det.par    - parameters for the detection function
#  det.formula- formula of covariates to use for scale of distance
#               if det.formula=~-1, uses a strip transect
#  covariates - dataframe of covariates across study area
#
#  Value: observation dataframe
#
# Jeff Laake
# 10 April 2008
################################################################################
{
   cov.obs=function(pts, covariates, covnms)
   {
     if(pts$n>0)
     {
       covdata =  lapply(covnms,
               function(nm, covariates){assign(nm, covariates[[nm]][pts])},
                              covariates=covariates)
       covdata=as.data.frame(do.call("cbind",covdata))
       names(covdata) = covnms
       return(covdata)
     }
     else
         return(NULL)
   }
   sample.points.from.line<-function(k,transects,lines,points.ppp,
                                       detfct,det.par,scale)
   {
#     Extract points within the strip
      which.instrip=inside.owin(points.ppp$x,points.ppp$y,owin(poly=transects[[k]]))
      points.instrip=points.ppp[owin(poly=transects[[k]])]
#     Compute values of distance from line
      distance=dist2line(points.instrip,lines[k,c("x0","y0","x1","y1")])$distance
#     Determine which points are observed; if formula =~-1 then it is a
#     a strip transect so any in the strip are seen
      if(is.null(scale))
      {
         if(points.instrip$n>0)
            which.seen=rep(TRUE,points.instrip$n)
         else
            return(NULL)
      }
      else
      {
#     If there are any points in the strip continue otherwise assign 0 to count
         if(length(distance)>0)
            which.seen=runif(length(distance)) <
                      detfct(distance,scale[which.instrip,,drop=FALSE]%*%det.par)
         else
            return(NULL)
       }
#      Return matrix of points
       sampled.points=points.instrip[which.seen]
       nobs=sum(as.numeric(which.seen))
       if(nobs==0 | length(which.seen)==0)
           return(NULL)
       else
       {
           sampled.points=cbind(label=rep(lines[k,"label"],nobs),
                             x=sampled.points$x,y=sampled.points$y,distance=distance[which.seen])
           row.names(sampled.points)=NULL
           return(sampled.points)
       }
   }
   varnames=all.vars(det.formula)
   covnms=names(covariates)
   if(any(!varnames %in% covnms)) stop("\n variable in det.formula is not in covariates\n")
   if(det.formula==~-1)
      scale=NULL
   else
   {
      covariates.im <- create.covariate.images(covariates,varnames)
      if(det.formula==~1 )
         scale=matrix(1,ncol=1,nrow=points.ppp$n)
      else
      {
         cov.df=cov.obs(pts=points.ppp,covariates=covariates.im,covnms=varnames)
         scale=model.matrix(det.formula,cov.df)
         if(dim(scale)[2]!=length(det.par)) stop("\n inconsistent det.formula and par vector\n")
      }
   }
   return(as.data.frame(do.call("rbind",lapply(1:length(transects),sample.points.from.line,
             transects=transects,lines=lines,points.ppp=points.ppp,
             detfct=detfct,det.par=det.par,scale=scale))))
}
hndetfct=function(x,scale)
{
  return(exp(-(x^2/(2*exp(scale)^2))))
}

