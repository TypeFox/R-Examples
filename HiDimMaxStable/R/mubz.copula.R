
mubz.copula.integrand<-Vectorize(function(g,b,z,
                                 params,
                                 copulas,
                                 margins,
                                 classes) {
    y<-prod(Vectorize(function(i) {
        pri.i<-copulas[[i]]@iPsi(margins[[i]]@p(g*z[classes==i],params[i*2]),
                             params[i*2-1])

        if(any(b[classes==i])) {
            copulas[[i]]@absdPsi(sum(pri.i),
                             params[i*2-1],sum(b[classes==i])) *
               prod( g * margins[[i]]@d(g*z[b&(classes==i)],params[i*2]) /
                    copulas[[i]]@absdPsi(pri.i[b[classes==i]],
                                         params[i*2-1],1) )
        } else {
            copulas[[i]]@psi(sum(pri.i),params[i*2-1])
        }
    })(1:max(classes)))
    ifelse(is.finite(y),y,0)
},vectorize.args=c("g"))

mubz.copula<-function(details=FALSE,...) {
    r<-integrate(mubz.copula.integrand,
                 0,Inf,
                 subdivisions=10000,stop.on.error=FALSE,...)
    ifelse(details,r,r$value)
}
