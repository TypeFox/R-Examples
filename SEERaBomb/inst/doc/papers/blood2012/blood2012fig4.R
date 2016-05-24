aBombHome="~/data/abomb"
load(file.path(aBombHome,"Hema87.RData"));
names(d)
d=d[d$py>0,] #remove two recs with zero py
d=d[d$kerma==1,] # take only kerma < 4 Gy
d$py=10^4*d$py
d$calg=as.integer(cut(d$calg,c(0,2,4,6,8,10)))
head(d)
m=d[d$sex==1,]; f=d[d$sex==2,] 

flin<-function(x,df) {
  c1=x[1];k=x[2];L=x[3:7];   		
  with(df,{mn = exp(c1+k*age)*py + sv*exp(L[calg])*py;
           -sum(CML*log(mn) - mn)})	}

# the idea here is to let the data speak through L, a bit like a one way anova
X0=c(c1=-13,k=0.05,rep(-10,5))
sol=optim(X0,flin,df=m,method="L-BFGS-B",hessian=TRUE,control=list(maxit=400))  
waitm=exp(sol$par[3:7])
sol=optim(X0,flin,df=f,method="L-BFGS-B",hessian=TRUE,control=list(maxit=400))  
waitf=exp(sol$par[3:7])
# if(length(grep("linux",R.Version()$os))) windows <- function( ... ) X11( ... )
# windows(width=5,height=5)
# dev.new(width=5,height=5)

par(mfrow=c(1,1),mar=c(4.7,0,2.3,0.2),lwd=3,cex.lab=1.8,cex.axis=1.7,cex.main=1.7,oma=c(0,4.5,1,0.7))
years=c(1952,1958,1966,1976,1985)-1945
plot(years,waitm,cex=2,pch=1,xlab="Years Since Exposure",ylab="",yaxt="n",col="blue",ylim=c(0,8.2e-4))
points(years,waitf,pch=2,cex=2,col="red")
legend(25,8.6e-4,c("Males","Females"),pch=1:2,col=c("blue","red"),cex=1.8,bty="n")
mtext(expression(paste("Cases per ",10^4," Person-Year-Sv")),side=2,line=2.4,cex=1.8,outer=T)
axis(side=2,las=1, at=c(0,2e-4,4e-4,6e-4,8e-4),labels=c(0,2,4,6,8),outer=T)
title("IR-to-CML Latency",outer=T,line=-1)
par(mar=c(5.1,4.1,4.1,2.1),oma=c(0,0,0,0)) # reset to standards



