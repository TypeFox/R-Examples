Ii<-function(th,it,model=NULL,D=1){
pr<-Pi(th,it,model=model,D=D)
P<-pr$Pi
dP<-pr$dPi
d2P<-pr$d2Pi
d3P<-pr$d3Pi
if (is.null(model)){
Q<-1-P
Ii<-dP^2/(P*Q)
dIi<-dP*(2*P*Q*d2P-dP^2*(Q-P))/(P^2*Q^2)
d2Ii<-(2*P*Q*(d2P^2+dP*d3P)-2*dP^2*d2P*(Q-P))/(P^2*Q^2)-(3*P^2*Q*dP^2*d2P-P*dP^4*(2*Q-P))/(P^4*Q^2)+(3*P*Q^2*dP^2*d2P-Q*dP^4*(Q-2*P))/(P^2*Q^4)
}
else{
pr0<-dP^2/P
pr1<-2*dP*d2P/P-dP^3/P^2
pr2<-(2*d2P^2+2*dP*d3P)/P-2*dP^2*d2P/-3*dP*d2P/P^2+2*dP^4/P^3
Ii<-as.numeric(rowSums(pr0,na.rm=TRUE))
dIi<-as.numeric(rowSums(pr1,na.rm=TRUE))
d2Ii<-as.numeric(rowSums(pr2,na.rm=TRUE))
}
res<-list(Ii=Ii,dIi=dIi,d2Ii=d2Ii)
return(res)}

