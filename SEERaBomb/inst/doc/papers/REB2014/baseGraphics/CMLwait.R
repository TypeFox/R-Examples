# Figure 2 of active MS  (differs from Blood Fig 4 only in computing mean T and M/F)
aBombHome="/data/abomb"
cols=c("city","sex","doseg","agexg","calg","kerma","PY","adjPY","num.entering",
		"age","agex","tsx","cal","sv","gam","neut","lymphoma","NHL","leukemia","AML","ALL","CML","ATL","MM")      
d<-read.table(file.path(aBombHome,"hema87.dat"), header=F,col.names=cols);
d=d[d$adjPY>0,] #remove two recs with zero py
d=d[d$kerma==1,] # take only kerma < 4 Gy
d$py=10^4*d$adjPY
m=d[d$sex==1,]; f=d[d$sex==2,] 
years=c(1951,1953,1956,1959,1963,1968,1973,1978,1983,1987)-1945
agem=0
flin<-function(x,df) {
	c0=x[1];k=x[2];L=x[3:12]; # let data speak through L, like a one way anova
	with(df,{mn = exp(c0+k*(age-agem))*py + sv*exp(L[calg])*py;
				-sum(CML*log(mn) - mn)})	}

X0=c(c0=ifelse(agem==0,-13,-10),k=0.05,rep(-10,10))
solm=optim(X0,flin,df=m,method="L-BFGS-B",hessian=TRUE,control=list(maxit=400))  
waitm=exp(solm$par[3:12])
solf=optim(X0,flin,df=f,method="L-BFGS-B",hessian=TRUE,control=list(maxit=400))  
waitf=exp(solf$par[3:12])
graphics.off()
windows(width=5,height=5)
par(mfrow=c(1,1),mar=c(4.7,0,2.3,0.2),lwd=3,cex.lab=1.8,cex.axis=1.7,cex.main=1.7,oma=c(0,4.5,1,0.7))
plot(years,waitm,cex=2,pch=1,xlab="Years Since Exposure",ylab="",yaxt="n",col="blue",ylim=c(0,8.2e-4))
points(years,waitf,pch=2,cex=2,col="red")
legend(25,8.6e-4,c("Males","Females"),pch=1:2,col=c("blue","red"),cex=1.8,bty="n")
mtext(expression(paste("Cases per ",10^4," Person-Year-Sv")),side=2,line=2.6,cex=1.8,outer=T)
axis(side=2,las=1, at=c(0,2e-4,4e-4,6e-4,8e-4),labels=c(0,2,4,6,8),outer=T)
title("IR-to-CML Latency",outer=T,line=-1)
print(sum(waitm)/sum(waitf)); print(MovF<-format(sum(waitm)/sum(waitf),digits=2))
pm=waitm/sum(waitm)
pf=waitf/sum(waitf)
print(taum<-sprintf("%s yrs",format(years%*%pm,digits=3)))
print(tauf<-sprintf("%s yrs",format(years%*%pf,digits=4)))
mtext(bquote(tau[m] == .(taum)),side=3,line=-6,cex=1.5,adj=.9,col="blue")
mtext(bquote(tau[f] == .(tauf)),side=3,line=-7.5,cex=1.5,adj=.9,col="red")
mtext(bquote(M/F == .(MovF)),side=3,line=-9.3,cex=1.3,adj=.9)
par(mar=c(5.1,4.1,4.1,2.1),oma=c(0,0,0,0)) # reset to standards