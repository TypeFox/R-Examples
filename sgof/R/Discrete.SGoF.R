
Discrete.SGoF <- function(u,pCDFlist=NA,K=NA, alpha = 0.05, gamma = 0.05, method=NA, Discrete=TRUE, Sides=1,...) {



discrete.sgof <- function(u,pCDFlist=NA, K=NA, alpha = 0.05, gamma = 0.05, method=NA, Discrete=TRUE, Sides=1,...) {




robust.fdr<-function(p,sides=1,p2=1-p,discrete=F,use8=T)

{
	m<-length(p)
	ord<-order(p)
	pcumtab<-cumsum(table(p))/m
	F05<-mean(p<=0.5)	
	edf<-approx(as.numeric(names(pcumtab)),pcumtab,xout=p,rule=2)$y
	if (sides==2)
	{
		pi<-min(1,2*mean(p))
		loc.fdr<-pi*p/edf
	}
	else
	{
		p.new<-2*(p*(p<=p2)+p2*(p>p2))
		pi<-min(1,2*mean(p.new))
		if (discrete) 
		{
			if (use8) pi<-min(1,8*mean(p.new))
			else
			{
				lam<-max(p[p<=0.5])
				k<-1/(lam^2+(0.5-lam)^2)
				pi<-min(k*mean(p.new),1)
			
			}
		}
		loc.fdr<-pi*p/edf
		loc.fdr[p>0.5]<-(0.5*pi+edf[p>0.5]-F05)/edf[p>0.5]
	}
	
	
	return(list(p=p,loc.fdr=loc.fdr,pi=pi,ord=ord))
}




n=length(u)

if(is.na(K)==TRUE || (is.na(K)==FALSE&K<n)){stepf <- lapply(pCDFlist, function(x) stepfun(x, c(0, x)))}

if(is.na(K)==TRUE){prob.vec<-unlist(lapply(stepf,function(s){s(gamma)}))}else{if(is.na(K)==FALSE&K<n){prob.vec<-c(unlist(lapply(stepf,function(s){s(gamma)})),rep(gamma,K))}else{prob.vec=rep(gamma,n)}}

crit.value.disc<-numeric(n)

b=seq(0,n)



if(is.na(method)==TRUE&n<2000){crit.value.disc<- b[min(which((1 -ppoibin(b , prob.vec,method = "DFT-CF")) <= alpha))]}
if(is.na(method)==TRUE&n>=2000){crit.value.disc<- b[min(which((1 -ppoibin(b , prob.vec,method = "RNA")) <= alpha))]}
if(is.na(method)==FALSE&method=="DFT-CF"){crit.value.disc<- b[min(which((1 -ppoibin(b , prob.vec,method = "DFT-CF")) <= alpha))]}
if(is.na(method)==FALSE&method=="RNA"){crit.value.disc<- b[min(which((1 -ppoibin(b , prob.vec,method = "RNA")) <= alpha))]}


num.reject.0<-length(u[u<=gamma])  
rejections.disc<-max(num.reject.0-crit.value.disc,0)

Discrete.SGoF = min(rejections.disc, sum(as.integer(n *ecdf(u)(u)) <= rejections.disc))
        
su <- sort(u)
jj <- which(u == 1)
if (length(jj) != 0) pi0 <- 1 else pi0 <- min((-1/n) * sum(log(1 - u)),1)

if (Discrete.SGoF == 0) {
FDR_DS <- 0
}else {
FDR_DS <-round(sort(robust.fdr(su, sides=Sides, p2 = 1 - su, discrete=Discrete , use8=TRUE)$loc.fdr)[Discrete.SGoF],4)
}




return(c(list(Rejections = Discrete.SGoF, FDR = min(FDR_DS,1))))
}

n=length(u)

if(missing(u)){stop("data argument is required")}

if(missing(pCDFlist)& (is.na(K)==TRUE || (is.na(K)==FALSE & K<n))){stop("pCDFlist argument is required")}



if(is.na(method)==F&(method!="RNA"&method!="DFT-CF")){stop("The specified method in incorrect")}


if(length(is.na(pCDFlist))==1&K<length(u)){stop("The specified K is incorrect")}

u<-as.vector(u)
res<-discrete.sgof(u,pCDFlist,K, alpha, gamma,method,Discrete,Sides,...) 

 
res$pvalues<-u
res$alpha<-alpha
res$gamma<-gamma
res$K<-K
res$method<-method
res$Discrete<-Discrete
res$Sides<-Sides
res$call<-match.call()
class(res)<-"Discrete.SGoF"
return(res)
}


