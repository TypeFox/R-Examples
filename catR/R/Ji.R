Ji<-function(th,it,model=NULL, D=1){
pr<-Pi(th,it,model=model,D=D)
P <- pr$Pi
dP <- pr$dPi
d2P <- pr$d2Pi
d3P <- pr$d3Pi
if (is.null(model)){
Q <-1-P
Ji <- dP * d2P/(P * Q)
dJi <- (P * Q * (d2P^2 + dP * d3P) - dP^2 * d2P *(Q - P))/(P^2 * Q^2)
res <- list(Ji = Ji, dJi = dJi)
}
else{
prov<-dP*d2P/P
prov1<-(P*d2P^2+P*dP*d3P-dP^2*d2P)/P^2
res<-as.numeric(rowSums(prov,na.rm=TRUE))
resd<-as.numeric(rowSums(prov1,na.rm=TRUE))
res<-list(Ji=res,dJi=resd)
}
return(res)
}

