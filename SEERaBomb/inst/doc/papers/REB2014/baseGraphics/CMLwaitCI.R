# Figure S1: Plots the CI of the Fk estimates in Fig. 2
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
getCI<-function(sol) {
	if (det(sol$hessian)>0) {
		sig=sqrt(diag(solve(sol$hessian)))
		print("Hessian OK") 
	} else {
		sig=Inf
		print("Hessian Singular") 
	}
	upper=signif(sol$par+1.96*sig,5)
	lower=signif(sol$par-1.96*sig,5)
	point=signif(sol$par,5)
	CI=cbind(point,lower,upper)
	CI=exp(CI) 
	row.names(CI)=1:nrow(CI)
	CI=transform(CI,li=point-lower,ui=upper-point,len=upper-lower)
	CI
}
(mCI=getCI(solm))
(fCI=getCI(solf))

mui=mCI[3:12,"ui"]
mli=mCI[3:12,"li"]
mpt=mCI[3:12,"point"]

fui=fCI[3:12,"ui"]
fli=fCI[3:12,"li"]
fpt=fCI[3:12,"point"]

require(plotrix)
windows(width=5,height=5)
par(mfrow=c(1,1),mar=c(4.7,4.7,2.3,0.2),lwd=2,cex.lab=1.8,cex.axis=1.7,cex.main=1.7)
plotCI(x=years,y=mpt,uiw=mui,liw=mli,ylab="",yaxt="n",col="blue",xlab="Years Since Exposure",ylim=c(0,8.2e-4),pch="x",cex=2)
plotCI(x=years+.3,y=fpt,uiw=fui,liw=fli,col="red",add=T,pch="o",cex=2)
axis(side=2,las=1, at=c(0,2e-4,4e-4,6e-4,8e-4),labels=c(0,2,4,6,8))
mtext(expression(paste("Cases per ",10^4," Person-Year-Sv")),side=2,line=2.6,cex=1.8)
mtext("Males (x)",side=3,line=.6,cex=1.8,col="blue",adj=0.05)
mtext("Females (o)",side=3,line=.6,cex=1.8,col="red",adj=0.95)
