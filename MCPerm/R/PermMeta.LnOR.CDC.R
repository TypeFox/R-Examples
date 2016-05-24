PermMeta.LnOR.CDC <-
function(PermMeta,plot_study="all",nrow=2,ncol=3,
    PermMeta.LnOR_pch=4,PermMeta.LnOR_col="black",
    LnOR_VAR_pch=18,LnOR_VAR_col="blue",
	 VAR_LnOR_pch=18,VAR_LnOR_col="red",
	 main="cumulative distribution curve for LnOR",
	 title=NULL,xlab="LnOR",ylab="cumulative probability",digits=3){
    if(class(PermMeta)!='PermMeta'){
	    stop("'PermMeta.LnOR.CDC' applied to an object of class 'PermMeta'.")
	 }
	 if(length(plot_study)==1 && plot_study=="all"){
	    study_num=PermMeta$study_num
		 plot_study=1:study_num
	 }else{
	    study_num=length(plot_study)
	 }
	 if(length(title)==0){
	    title=paste("study",plot_study,sep=" ")
	 }
	 plot_data=PermMeta$perm_LnOR[plot_study,,drop=FALSE]
	 LnOR=PermMeta$true_LnOR[plot_study]
	 LnOR_VAR=PermMeta$true_VARLnOR[plot_study]
	 sample=PermMeta$sample[plot_study]
	 VAR_LnOR=apply(plot_data,1,var)
	 repeatNum=PermMeta$repeatNum
	 opar=par(no.readonly=TRUE)
	 window=ceiling(study_num/(nrow*ncol))
	 j=1
	 for(i in 1:window){
	    par(mfrow=c(nrow,ncol))
		 if(length(main)!=0){
	       par(oma=c(0, 0, 2, 0))
	    }
	    while(ceiling(j/(nrow*ncol))==i){
		    xlim_value=range(plot_data[j,])
		    plot(ecdf(plot_data[j,]),ylim=c(0,1),
			    pch=PermMeta.LnOR_pch,cex=0.5,xlim=xlim_value,col=PermMeta.LnOR_col,
				 main=title[j],xlab=xlab,ylab=ylab)
			 
			 par(new=TRUE)
			 plot(pch=LnOR_VAR_pch,cex=0.5,
			    ecdf(rnorm(repeatNum,mean=0,sd=sqrt(LnOR_VAR[j]))),
			    col=LnOR_VAR_col,
			    xlim=xlim_value,ylim=c(0,1),main="",xlab='',ylab="")
	         
			 par(new=TRUE)
			 plot(pch=VAR_LnOR_pch,cex=0.5,
			    ecdf(rnorm(repeatNum,mean=0,sd=sqrt(VAR_LnOR[j]))),
				 col=VAR_LnOR_col,
			    xlim=xlim_value,ylim=c(0,1),main="",xlab='',ylab="")
			 
			 legend("topleft",
	          legend=c('perm_lnOR','pnorm_LnOR_VAR','pnorm_VAR_LnOR'),
			    pch=c(PermMeta.LnOR_pch,LnOR_VAR_pch,VAR_LnOR_pch),
		       col=c(PermMeta.LnOR_col,LnOR_VAR_col,VAR_LnOR_col),
				 border=c('white','white','white'),
				 text.col=c(PermMeta.LnOR_col,LnOR_VAR_col,VAR_LnOR_col),
				 cex=0.6)
			 legend("bottomright",
	             legend=c(paste("LnOR",round(LnOR[j],digits),sep="="),
					    paste("sample",sample[j],sep="="),
				       paste("LnOR_VAR",round(LnOR_VAR[j],digits),sep="="),
				       paste("VAR_LnOR",round(VAR_LnOR[j],digits),sep="=")),
					 pch=c(PermMeta.LnOR_pch,PermMeta.LnOR_pch,LnOR_VAR_pch,VAR_LnOR_pch),
				    col=c(PermMeta.LnOR_col,PermMeta.LnOR_col,LnOR_VAR_col,VAR_LnOR_col),
                text.col=c(PermMeta.LnOR_col,PermMeta.LnOR_col,LnOR_VAR_col,VAR_LnOR_col),
		          cex=0.6,ncol=1)
					 
          j=j+1		
          if(j>study_num){
			    break
			 }			 
		 }
		 mtext(text=main,side=3,outer=TRUE)	
		 if(j>study_num){
		     break
		 }
		 x11()
	    par(mfrow=c(nrow,ncol))	 
	 } 
	 par(opar)
	 result=list('plot_study'=plot_study,'LnOR'=LnOR,'sample'=sample,
	    'LnOR_VAR'=LnOR_VAR,"VAR_LnOR"=VAR_LnOR)
	 invisible(result)
}
