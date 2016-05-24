
setClass("spatial",
         slots=list(
             name="character",
             p.start="numeric",
             p.min="function",p.max="function",
             vario="logical",
             cov.function="function"))

spatialWhittleMatern<-new("spatial",
                   name="WhittleMatern",
                   p.start=c(1,1,1),
                   p.min=function(d) c(0,0,0),
                   p.max=function(d) c(Inf,Inf,Inf),
                   cov.function=function(h,params) {
                       params[3]*ifelse(h==0,1,
                              (2^(1-params[2])/gamma(params[2]))*((h/params[1])^params[2])*
                              besselK((h/params[1]),params[2]))
                   },vario=FALSE
                   )

spatialCauchy<-new("spatial",
            name="Cauchy",
            p.start=c(1,1,1),
            p.min=function(d) c(0,0,0),
            p.max=function(d) c(Inf,Inf,Inf),
            cov.function=function(h,params) {
                params[3]*(1+(h/params[1])^2)^(-params[2])
            },vario=FALSE
            )

spatialPowerexp<-new("spatial",
           name="Powerexp",
           p.start=c(1,1,1),
           p.min=function(d) c(0,0,0),
           p.max=function(d) c(Inf,2,Inf),
           cov.function=function(h,params) {
               params[3]*exp(-(h/params[1])^params[2])
           },vario=FALSE
           )

spatialBessel<-new("spatial",
                   name="Bessel",
                   p.start=c(1,1,1),
                   p.min=function(d) c(0,(d-2)/2,0),
                   p.max=function(d) c(Inf,Inf,Inf),
                   cov.function=function(h,params) {
                       params[3]*ifelse(h==0,1,
                              ((2*params[1]/h)^params[2])*gamma(params[2]+1)*besselJ((h/params[1]),params[2]))
                   },vario=FALSE
                   )

spatialPower<-new("spatial",
                  name="Power",
                  p.start=c(1,1),
                  p.min=function(d) c(0,0),
                  p.max=function(d) c(Inf,2),
                  cov.function=function(h,params) {
                     (h/params[1])^params[2]
                  },vario=TRUE
                  )

var.from.points<-function(xy,family,params=NULL) {
    n.points<-dim(xy)[1]
    if(family@vario) {
        outer(1:n.points,1:n.points,
              function(i,j) {
                  family@cov.function(sqrt(sum(xy[i,]^2)),params)+
                      family@cov.function(sqrt(sum(xy[j,]^2)),params)-
                          family@cov.function(sqrt((xy[i,1]-xy[j,1])^2+(xy[i,2]-xy[j,2])^2),params)
              }
              )
    } else {
        family@cov.function(
            outer(1:n.points,1:n.points,
                  function(i,j)
                  sqrt((xy[i,1]-xy[j,1])^2+(xy[i,2]-xy[j,2])^2)
                  ),
            params)
    }
}

spatial.cor.matrix<-function(params,spatial) {
    var.from.points(spatial$sites,spatial$family,params)
}
