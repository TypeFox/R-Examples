`Amatdual2` <-
function (steps, pointsin, removelist, nbrs, weights, alpha) 
{

nn<-length(nbrs)
lpo<-length(pointsin)
lre<-length(removelist)
n<-lpo+lre

adual<-matrix(0,n-steps+1,n-steps+1)

tmp<-.C("amatdual",as.integer(steps),as.integer(pointsin),as.integer(removelist),as.integer(nbrs),as.double(weights),as.double(alpha),
		as.integer(lpo),as.integer(lre),as.integer(nn),Adual=as.double(t(adual)),PACKAGE="adlift")

adual<-matrix(tmp$Adual,byrow=TRUE,ncol=n-steps+1)

return(adual)
}

