### Main function which finds the A and B constants using mean-mean, 
### mean-sigma, haebara and SL. It additionally implements H and SL using
### an asymetric cloglog ICC

irt.link<-function(parm,common,model,icc,D,...)
UseMethod("irt.link")


irt.link.default<-function(parm,common,model,icc,D,...){
if(icc=="cloglog" & model!="1PL") stop("The cloglog is not yet implemented for this model")

cl<-match.call()

colnames(parm) <- c("aJj", "bJj", "cJj", "aIj", "bIj", "cIj")
parm = parm[common,]

Amm = mean(parm[,4:6]$aIj)/mean(parm[,1:3]$aJj)
Bmm = mean(parm[,1:3]$bJj) - Amm*mean(parm[,4:6]$bIj)

Ams = sd(parm[,1:3]$bJj)/sd(parm[,4:6]$bIj)
Bms = mean(parm[,1:3]$bJj) - Ams*mean(parm[,4:6]$bIj)

if(model=="1PL"){
Amm=1
Ams=1
Bmm = mean(parm[,1:3]$bJj) - Amm*mean(parm[,4:6]$bIj)
Bms = mean(parm[,1:3]$bJj) - Ams*mean(parm[,4:6]$bIj)
Hae<-optim(Bmm,target,parm=parm,common=common,model=model,
icc=icc,D=D,meth="H",method="BFGS")$par
Haebara<-c(1,Hae)
StLo<-optim(Bmm,target,parm=parm,common=common,model=model,
icc=icc,D=D,meth="SL",method="BFGS")$par
StockLord<-c(1,StLo)
}
else{
Haebara<-optim(c(Amm,Bmm),target,parm=parm,common=common,model=model,
icc=icc,D=D,meth="H",method="BFGS")$par

StockLord<-optim(c(Amm,Bmm),target,parm=parm,common=common,model=model,
icc=icc,D=D,meth="SL",method="BFGS")$par
}
res<-list(call=cl,mm=c(Amm,Bmm),ms=c(Ams,Bms),Haebara=Haebara,StockLord=StockLord)
class(res)<-"irt.link"
		   return(res)
}


print.irt.link<-function(x,...)
{
	vals<-rbind(c(x$mm[1],x$mm[2]),
			c(x$ms[1],x$ms[2]),
			c(x$Haebara[1],x$Haebara[2]),
			c(x$StockLord[1],x$StockLord[2])) 
	dimnames(vals)<-list(c("Mean-Mean","Mean-Sigma","Haebara",
				"Stocking-Lord"),c("A","B"))

	cat("\nCall:\n")
	print(x$call)
	cat("\nIRT parameter-linking constants:\n")	
	cat("\n")
	print(vals)
	cat("\n")
}	

