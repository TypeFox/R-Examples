PermMeta.LnOR.qqnorm <-
function(PermMeta,plot_study="all",nrow=2,ncol=2,
    main="qqnorm plot for LnOR",title=NULL,xlab="Theoretical Quantiles",ylab="Sample Quantiles",
	 scatter_col="black",line_col="red"){
    if(class(PermMeta)!='PermMeta'){
	    stop("Function 'PermMeta.LnOR.qqnorm' applied to an object of class 'PermMeta'.")
	 }
	 if(length(plot_study)==1 && plot_study=="all"){
	    study_num=PermMeta$study_num
		 plot_study=1:study_num
	 }else{
	    study_num=length(plot_study)
	 }
	 if(length(title)==0){
	    title=c()
		 for(i in 1:study_num){
		    title=c(title,paste("study",plot_study[i],sep=" "))
		 }
	 }
	 plot_data=PermMeta$perm_LnOR[plot_study,,drop=FALSE]
	 
	 opar=par(no.readonly=TRUE)
	 window=ceiling(study_num/(nrow*ncol))
	 j=1
	 for(i in 1:window){
	    par(mfrow=c(nrow,ncol))
	    if(length(main)!=0){
	       par(oma=c(0, 1, 3, 0))
	    }
	    while(ceiling(j/(nrow*ncol))==i){
		    qqnorm(plot_data[j,],main=title[j],xlab = xlab, ylab = ylab,col=scatter_col)
          qqline(plot_data[j,], distribution = qnorm,probs = c(0.25, 0.75), qtype = 7,col=line_col)  
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
}
