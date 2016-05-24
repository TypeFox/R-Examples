"s.estim" <-
function(l0,si,y,tlo,mxs) {
 r0 <- y - l0
 zs <- rysigm(r0,wgt=y,sigmai=2*si,np=1,itype=1,isigma=1,tol=tlo,maxis=mxs)
 c(zs$sigmaf,l0)}

