transect.intensity=function(x, epsvu=NULL, obs.ppp, covariates, nclass=NULL, width)
###########################################################################
# Computes distribution of expected and observed numbers of objects for
# nclass equally spaced bins of distance for each transect.
#
# Arguments:
#
#    x       - dpsat object
#   epsvu    - epsvu setting for fitted model (only need epsvu[2] value for u
#  obs.ppp   - observation point process
# covariates - dataframe of covariates at quadrature points
#   nclass   - number of distance intervals
#   width    - maximum width of the transects
#
# Value:
#
#   exp.counts - matrix of expected counts in each distance bin (columns)
#                   for each transect (row)
#   obs.counts - matrix of observed counts in each distance bin (columns)
#                   for each transect (row)
#
# Jeff Laake
# 4/28/2008
###########################################################################
{
# Check to make sure x is of class dspat
  if(class(x)[1]!="dspat") stop("\n Argument x must be of class dspat\n")
#
# compute equally-spaced break points for distance bins; following code
# makes sure that they are equal multiples of width/2 and epsvu[2]
#
  n=obs.ppp$n
  if(is.null(nclass))
    nclass=max(2,ceiling(sqrt(n)))
  width=width/2
  max.class=width/epsvu[2]
  cm=width/((1:max.class)*epsvu[2])-floor(width/((1:max.class)*epsvu[2])+1/max.class)
  xm=(1:max.class)[abs(cm)<1/max.class]-nclass
  if(!any(xm>0) || abs(max(xm[xm<=0])) < abs(min(xm[xm>0])))
  {
    nclass=nclass+max(xm[xm<=0])
  } else
  {
    nclass=nclass+min(xm[xm>0])
  }
  breaks=(0:nclass)*width/nclass
#
# compute expected intensity over each transect (lambda*g(u)) using
# the dummy quadrature points that are equally spaced based on eps for u
# also compute the freq of observations in the same bins
#
  nlines= dim(x$lines.psp$ends)[1]
  distribution=vector("list",length=nlines)
  observed=vector("list",length=nlines)
  Area=prod(epsvu)
  for (i in 1:nlines)
  {
     inside=inside.owin(x$model$Q$dummy$x,x$model$Q$dummy$y,owin(poly=x$transects[[i]]))
     locs=x$model$Q$dummy[inside]
     distances=dist2line(locs,x$lines.psp$ends[i,])$distance
     lambda <- predict(x$model,locations=locs,
               covariates=covariates[(n+1):dim(covariates)[1],][inside,],type="lambda")
     distribution[[i]]=data.frame(label=rep(x$lines.psp$label[i],
                              length(distances)),N=lambda*Area,distance=distances)
     distances=dist2line(obs.ppp[owin(poly=x$transects[[i]])],x$lines.psp$ends[i,])$distance
     observed[[i]]=hist(distances,breaks=breaks,plot=FALSE)$counts
  }
  exp.counts=do.call("rbind",distribution)
  obs.counts=do.call("rbind",observed)
  exp.counts=tapply(exp.counts$N,list(exp.counts$label,cut(exp.counts$distance,breaks)),sum)
  exp.counts=exp.counts*sum(obs.counts)/sum(exp.counts)
  return(list(exp.counts=exp.counts,obs.counts=obs.counts))
}
