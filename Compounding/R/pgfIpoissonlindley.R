pgfIpoissonlindley <-
function(s,params) {
    xval<-length(s)
    theta<-params[1]
    for (i in 1:length(s)) {
        func<-function(x) pgfpoissonlindley(x,params)-s[i]
        xval[i]<-uniroot(func,lower=0,upper=1)$root
    }
    xval
}
