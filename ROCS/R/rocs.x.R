rocs.x <-
function(x0, x1, s0=NULL, s1=NULL, n.perm=1000, do.plot=TRUE, FDR.cut=0.2) # receiver operational characteristic surface
# x0: raw values of class 0
# x1: raw values of class 1
# s0: confidence levels of class 0; between 0 and 1
# s1: confidence levels of class 1; between 0 and 1
{	
	if(is.null(s0))
	{
		s0<-rep(1, length(x0))
		s1<-rep(1, length(x1))
	}
	
	if(median(x1)<median(x0))
	{
		x1<- -x1
		x0<- -x0
	}
	
	x<-c(x0, x1)
	s<-c(s0, s1)
	y<-c(rep(0, length(x0)), rep(1,length(x1)))
	
	o<-order(x)
	s<-s[o]
	y<-y[o]
	x<-x[o]
	
	get.fptptdr<-function(x, s, y)
	{
		s[y==0]<-1-s[y==0]  # now s is the probability of each x being in the 1 class
		
		l<-length(s)		
		
		dd<-cumsum(s)  # false negative count
		aa<-cumsum(1-s)  # true negative count
		bb<-sum(1-s)-aa  # false positive count
		cc<-sum(s)-dd  # true positive count
		
		to.keep<-c(1, 1-(x[1:(l-1)]==x[2:l]))
		sel<-which(to.keep==1)	
		dd<-dd[sel]
		aa<-aa[sel]
		bb<-bb[sel]
		cc<-cc[sel]
		
		FP<-bb/sum(1-s)
		TP<-cc/sum(s)
		TDR<-cc/(cc+bb)
        
        sel<-which(!is.na(FP) & !is.na(TP) & !is.na(TDR))
        FP<-FP[sel]
        TP<-TP[sel]
        TDR<-TDR[sel]
		
		FP<-c(1, FP, 0)
		TP<-c(1, TP, 0)
		TDR<-c(min(TDR), TDR, 1)
		
		return(list(FP=FP, TP=TP, TDR=TDR))
	}
	
	
	get.vus<-function(x, s, y) # x are values ordered from smallest to largest, y are group labels
	{
		FP.TP.TDR<-get.fptptdr(x, s, y)
		FP<-FP.TP.TDR$FP
		TP<-FP.TP.TDR$TP
		TDR<-FP.TP.TDR$TDR
		l<-length(FP)
		
		abs.d.TP<- -diff(TP)
		mid.FP<- FP[1:(l-1)]+diff(FP)/2	
		mid.TDR<-TDR[1:(l-1)]+diff(TDR)/2
		vus<-sum(abs.d.TP*(1-mid.FP)*mid.TDR)	
		vus
	}
	
	
	vus<-get.vus(x, s, y)
	
	vus.perm<-matrix(0, ncol=2, nrow=n.perm)
	for(i in 1:n.perm)
	{
		new.o<-sample(1:length(y), length(y), replace=FALSE)
		new.y<-y[new.o]
		new.s<-s[new.o]
		vus.perm[i,]<-unlist(get.vus(x, new.s, new.y))
	}
	
#### generate plot
	
	fptptdr<-get.fptptdr(x, s, y)
	if(do.plot) 
    {
        rocs.fptp(fptptdr$FP, fptptdr$TP, fptptdr$TDR, FDR.cut=FDR.cut)
        fcauc.fptp(fptptdr$FP, fptptdr$TP, fptptdr$TDR, FDR.cut=FDR.cut)
    }
	
	fcauc<-fcauc.fptp(fptptdr$FP, fptptdr$TP, fptptdr$TDR, FDR.cut=FDR.cut, do.plot=FALSE)
#### return volume under surface
	list(vus=vus, fcauc=fcauc, vus.perm.pval=sum(vus.perm[,1] >= vus)/n.perm, n.perm=n.perm, fp=fptptdr$FP, tp=fptptdr$TP, tdr=fptptdr$TDR)
}
