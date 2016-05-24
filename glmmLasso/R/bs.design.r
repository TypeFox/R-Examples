bs.design<-function(x, diff.ord, spline.degree, knots.no)
{
  
  ## generate a B-Spline-Matrix with equidistant knots (code by Thomas Kneib):
  n<-length(x)
  xl<-min(x)
  xr<-max(x)
  xmin<-xl-(xr-xl)/100
  xmax<-xr+(xr-xl)/100
  dx<-(xmax-xmin)/(knots.no-1)
  knots<-seq(xmin-spline.degree*dx,xmax+spline.degree*dx,by=dx)
  B<-splines::spline.des(knots,x,spline.degree+1)$design
  knots<-knots[1:(knots.no+spline.degree-1)]
  
  ## generate Penalization-Matrix
  D<-diag(length(knots))
  d<-min(diff.ord,spline.degree)
  if(d<diff.ord) warning(paste("order of differences > degree of
                               splines:\n new order of differences=",d,"\n"))
  if(d>0) {for(i in 1:d) D<-diff(D)}
  
  ## reparametrization: B_unpen = unpenalized part, B_pen=penalized part
  B.unpen.fact<-rep(1,length(knots))
  if(diff.ord>1) {for(i in 2:diff.ord)
    B.unpen.fact<-cbind(B.unpen.fact,knots^(i-1)) }
  
  B.unpen<-(B%*%B.unpen.fact)
  B.pen  <-B%*%t(D)%*%solve(D%*%t(D))
  
  return(list(B=B, X=B.unpen, Z=B.pen))
}
