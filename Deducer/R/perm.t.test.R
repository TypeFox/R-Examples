perm.t.test<-function(x,y,statistic=c("t","mean"),
			alternative=c("two.sided", "less", "greater"), midp=TRUE, B=10000){
	DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
	x<-na.omit(x)
	y<-na.omit(y)
	nx<-length(x)
	ny<-length(y)
	alternative<-match.arg(alternative)
	statistic<-match.arg(statistic)
	mult<- if(midp) .5 else 1
	var1<-var#function(s){.Internal(cov(s, NULL, TRUE, FALSE))}
	sum1<-.Primitive("sum")
	dat<-c(x,y)	
	if(statistic=="mean"){
		METHOD<-paste("Two-Sample permutation test using mean difference (B=", B,")")
		stat<-sum(y)
		nmin<-ny
		samp<-replicate(B,sum(sample(dat,nmin)))
		STAT<-mean(x)-mean(y)
		names(STAT)<-"Mean Difference"
	}else{
		METHOD<-paste("Two-Sample permutation test using Welsh's t (B=",B,")")
		stat<-(mean(x)- mean(y))/sqrt(var1(x)/nx+var1(y)/ny)
		samp<-replicate(B,sample(dat,nx+ny))
		sx<-apply(samp[1:nx,],2,function(s) c(sum1(s),var1(s)))
		sy<-apply(samp[(nx+1):(nx+ny),],2,function(s) c(sum1(s),var1(s)))
		samp<-(sx[1,]/nx-sy[1,]/ny)/sqrt(sx[2,]/nx+sy[2,]/ny)
		STAT<-stat
		names(STAT)<-"Welsh t-statistic"
	}
	lower<-sum(samp<stat)/(B+1)
	upper<-sum(samp>stat)/(B+1)
	equal<-sum(samp==stat)/(B+1)
	if(alternative=="two.sided")
		p.value<-2*min(lower,upper)+2*mult*equal
	else if(alternative=="less")
		p.value<-lower+mult*equal
	else
		p.value<-upper+mult*equal
	p.value<-min(p.value,1)
	RVAL<-list(statistic=STAT,p.value=p.value,method=METHOD,data.name=DNAME,
				alternative=alternative,B=B)
	class(RVAL)<-"htest"
	return(RVAL)
}
