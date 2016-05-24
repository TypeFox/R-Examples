"concorgm" <-
function(x,px,y,py,r) { 
if (sum(px) != dim(x)[2] | sum(py) != dim(y)[2] ) stop("px or py IS NOT SUITABLE")
s<-svdbip(t(x)%*%y,px,py,r)
list(u=s$u,v=s$v,cov2=s$s2/dim(x)[1]^2)
}

