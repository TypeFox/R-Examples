setClass("pssm",representation(call="call",
                               convergence="numeric",
                               loglike="function",
                               estimates="numeric",se.estimates="numeric",covariance.estimates="matrix",
                               estimates.progression="numeric",se.estimates.progression="numeric",
                               estimates.survival="numeric",se.estimates.survival="numeric",
                               hazard.progression="numeric",hazard.survival="numeric",
                               intervals="integer",rescale="numeric",
                               formula.progression="formula",formula.survival="formula",
                               progression.covariate.list="character",
                               survival.covariate.list="character",
                               message="character"))
setClass("pssm.summary",representation(call="character",
                                       convergence="character",
                                       coefficients="data.frame",
                                       confidence.bounds="data.frame"))
setMethod("print",signature(x="pssm"),
function(x, ...){
	lst1<-lst2<-NULL
	label1<-label2<-NULL
	cat("\nCall:\n", paste(deparse(x@call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
	if (x@convergence>0) cat("Did not Converge\n")
	cat("Coeficients:\n")
	m1=length(x@progression.covariate.list)
	if (m1>0){
	lst1=c(x@estimates.progression)
    label1=paste(rep("progression.",m1),x@progression.covariate.list,sep="")
}
m2=length(x@survival.covariate.list)
if(m2>0) {
	lst2=c(x@estimates.survival)
    label2=paste(rep("survival.",m2),x@survival.covariate.list,sep="")
}	
	
	
if (m1+m2==0){cat("No Coefficients\n")}
else {
		lst=c(lst1,lst2)
	names(lst)=c(label1,label2)
	print(format(lst,digits=4),print.gap=2,quote=FALSE)

}
}
)

setMethod("summary",
    signature(object = "pssm"),
   function(object, ...){
	x=object
	ests1<-ests2<-NULL
	label1<-label2<-NULL
	se1<-se2<-NULL
	call1=paste("Call: ",deparse(x@call), sep = "")
	if (x@convergence>0) converge="Convergence: No" else converge="Convergence: Yes"
	#cat("Coeficients:\n")
	m1=length(x@progression.covariate.list)
if (m1>0){
	ests1=c(x@estimates.progression)
	se1=c(x@se.estimates.progression)
    label1=paste(rep("progression.",m1),x@progression.covariate.list,sep="")
}
m2=length(x@survival.covariate.list)
if(m2>0) {
	ests2=c(x@estimates.survival)
	se2=c(x@se.estimates.survival)
    label2=paste(rep("survival.",m2),x@survival.covariate.list,sep="")
}
if (m1+m2==0){
  #cat("No Coefficients\n")
  coef=NULL
}
else {		
	lst=c(ests1,ests2)
	se=c(se1,se2)
	out1=data.frame(lst,se,lst/se,2*pnorm(-abs(lst/se)),
                  row.names=c(label1,label2))
	names(out1)<-c("coef","se(coef)","z","Pr(>|z|)")
	#print(format(out1,digits=4),print.gap=2,quote=FALSE)
	out2=data.frame(exp(lst),exp(-lst),exp(lst-1.96*se),exp(lst+1.96*se),
                  row.names=c(label1,label2))
	names(out2)<-c("exp(coef)","exp(-coef)","lower 0.95","upper 0.95")
    #print(format(out2,digits=4),print.gap=2,quote=FALSE)
  coef=out2
}
#print(...)

returns=new("pssm.summary",call=call1,convergence=converge,coefficients=out1,
            confidence.bounds=out2)
return(returns)
}
)
setMethod("print",signature(x="pssm.summary"),
          function(x){
            cat(x@call,"\n")
            cat(x@convergence,"\n")
            print(format(x@coefficients,digits=4),print.gap=2,quote=FALSE)
            print(format(x@confidence.bounds,digits=4),print.gap=2,quote=FALSE)
          })
          
setMethod("plot",
    signature(x = "pssm"),
    function (x,type,cov1=NULL,cov2=NULL) {
 		#x=object
		 m=x@intervals
		 mc=as.list(x@call)   
		sp=is.null(mc$progr)
        ss=is.null(mc$survv)
		both=!(sp||ss)
		if (nargs()==1){
		m=x@intervals
		absc=rep(0,m+1)
		mss=ifelse(both,m+1,1)
		fcns=matrix(0,mss,m+1)
        sg=0
		if(both||ss) haz<-exp(x@hazard.progression) else haz<-exp(x@hazard.survival)
		for(i in 1:m){
			absc[i+1]=i/x@rescale
			fcns[1,i+1]=fcns[1,i]+haz[i]
			if(both){
			for (j in (i:m)) fcns[1+i,j+1]=fcns[1+i,j]+exp(x@hazard.survival[sg+j-(i-1)])
			sg=sg+(m-i+1)}
}
        fc=matrix(0,mss,m*5+1)
		ab=seq(0,m/x@rescale,1/(5*x@rescale))
	    for (i in (1:mss)){
 			fc[i,]=approx(absc,fcns[i,],ab)$y
		}
		plot(ab,exp(-fc[1,]),type='l',main='Time to Survival/Progression',xlab='Time',
			ylab='Probability')
		if(both){for (i in 2:(m+1)) lines(ab,exp(-fc[i,]),lty=2)}
		
        }
    #MORE Then ONe arguement
		else{
			fcf=pssm.survivalcurv(x,cov1,cov2,timeToProgression=(type=='progression'),covariance=FALSE)
			ab=seq(0,m,1/5)/x@rescale
			fc=fcf(ab)
			if (type=='progression'&&(!(sp||ss)))fc=fc[fc$rep=='s2',]
			
			covs=c(x@progression.covariate.list,x@survival.covariate.list)
			l1=length(x@progression.covariate.list)
			l2=length(x@survival.covariate.list)
			ms=length(covs)
			us=data.frame(unique(fc[,3:(3+ms-1)]))
			lus=dim(us)[1]
			plot(1,1,xlim=c(0,max(ab)),ylim=c(0,1),type="n",main=type,xlab='Time',
			ylab='Probability')
		   lab=rep("                   ",lus)
		   for (ii in (1:lus)){
		    keep=rowSums(as.matrix(fc[,3:(3+ms-1)])==matrix(us[ii,],dim(fc)[1],ms,byrow=TRUE))==length(3:(3+ms-1))
			lines(ab,fc$estimate[keep],type='l',lty=ii,col=ii)
        	lab[ii]=paste(names(fc)[3:(3+ms-1)],"=",us[ii,],collapse=' ')		
		}
		    legend(x='bottomleft',legend=lab,lty=1:lus)
				
    }
}
	)

