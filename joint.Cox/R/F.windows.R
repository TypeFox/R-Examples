F.windows <-
function(time,widths,X,Z1,Z2,beta1,beta2,eta,theta,alpha,
                     g,h,xi1,xi3,Fplot=TRUE){
wmax=max(widths)
m=length(widths)
F_event_at_X=F_noevent=numeric(m)

for(i in 1:m){
  F.res=F.window(time,widths[i],X,Z1,Z2,beta1,beta2,eta,theta,
                 alpha,g,h,xi1,xi3,Fplot=FALSE)
  F_event_at_X[i]=F.res["F_event_at_X"]
  F_noevent[i]=F.res["F_noevent"]
}

if(Fplot==TRUE){
  num_grid=500
  x_grid=seq(X,time,length=num_grid)
  plot(x_grid,rep(0,num_grid),xlim=c(xi1,xi3),ylim=c(-0.05,1.05),type="l",lwd=3,
       xlab="t+w",ylab="Probability of death in ( t, t+w )",col="red")
  abline(h=0)
  abline(v=time,col="gray")
  points(X,0,lwd=4,col="red")
  points(time+widths,F_event_at_X,col="red",type="o",pch=16,lwd=2,cex=1.1)
  points(time+widths,F_noevent,col="blue",type="o",pch=17,lwd=2,cex=1.1)
  text(X,0.05,"X",cex=1)
  text(time,-0.05,"t",cex=1)
}

cbind(t=time,w=widths,X=X,F_event_at_X=F_event_at_X,F_noevent=F_noevent)
}
  
  
