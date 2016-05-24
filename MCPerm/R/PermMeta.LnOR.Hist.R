PermMeta.LnOR.Hist <-
function(PermMeta,plot="LnOR",plot_study="all",nrow=2,ncol=2,
    main="Background distribution for LnOR",title=NULL,xlab="LnOR",
	 hist_border_col="black",arrows_col="red",digits=3){
    if(class(PermMeta)!='PermMeta'){
	    stop("Function 'PermMeta.LnOR.Hist' applied to an object of class 'PermMeta'.")
	 }
	 if (!is.element(plot, c("LnOR", "LnOR_VAR"))){
	     stop("'plot' should be 'LnOR', 'LnOR_VAR'.")
	 }
	 if(length(plot_study)==1 && plot_study=="all"){
	    study_num=PermMeta$study_num
		 plot_study=1:study_num
	 }else{
	    study_num=length(plot_study)
	 }
	 repeatNum=PermMeta$repeatNum
	 switch(plot,
	    LnOR={
		    plot_data=PermMeta$perm_LnOR[plot_study,,drop=FALSE]
			 true_data=PermMeta$true_LnOR[plot_study]
		 },
		 LnOR_VAR={
		    plot_data=PermMeta$perm_VARLnOR[plot_study,,drop=FALSE]
			 true_data=PermMeta$true_VARLnOR[plot_study]
		 }
	 )
	 
	 true_LnOR=PermMeta$true_LnOR[plot_study]
	 true_VARLnOR=PermMeta$true_VARLnOR[plot_study]
	 sample=PermMeta$sample[plot_study]
	 mean_value=apply(plot_data,1,mean)
	 var_value=apply(plot_data,1,var)
	 p_value=c()
	 for(i in 1:study_num){
	    temp=plot_data[i,]
		 temp=temp[temp>true_LnOR[i]]
		 p_value=c(p_value,((length(temp))/repeatNum))
	 }
	 if(length(title)==0){
	    title=c()
		 for(i in 1:study_num){
		    title=c(title,paste("study",plot_study[i],sep=" "))
		 }
	 }
	 
	 opar=par(no.readonly=TRUE)
	 window=ceiling(study_num/(nrow*ncol))
	 j=1
	 for(i in 1:window){
	    par(mfrow=c(nrow,ncol))
	    if(length(main)!=0){
	       par(oma=c(0, 1, 3, 0))
	    }
	    while(ceiling(j/(nrow*ncol))==i){
		    x_lim=range(plot_data[j,])
			 switch(plot,
			    LnOR={
				    x_lim=c(x_lim[1],x_lim[2]*1.7)
					 pch_value=c(2,20,20,20,20,20)
				 },
			    LnOR_VAR={
				    x_lim=c(x_lim[1],x_lim[2]*1.005)
					 pch_value=c(20,2,20,20,20,20)
				 }
			 )  
		    dens=hist(plot_data[j,],freq=FALSE,xlim=x_lim,main=title[j],xlab=xlab,
			    border=hist_border_col,font.main=1)$density
	       y1_value=max(dens)/10
	       points(true_data[j],y=-y1_value/5,pch=2,col=arrows_col,cex=0.6,lwd=1.3)
	       arrows(x0=true_data[j],y0=y1_value,x1=true_data[j],y1=0,
	          col=arrows_col,length=0.05,lwd=1.3)
	       legend("topright",
	          legend=c(paste('LnOR',round(true_LnOR[j],digits),sep="="),
				    paste('LnOR_VAR',round(true_VARLnOR[j],digits),sep="="),
					 paste('sample',sample[j],sep="="),
		          paste("mean",round(mean_value[j],digits),sep="="),
				    paste("var",round(var_value[j],digits),sep="="),
				    paste("p",round(p_value[j],digits),sep="=")), 
		       pch=pch_value,
		       col=c(arrows_col,arrows_col,arrows_col,'black','black','black'),
		       text.col=c(arrows_col,arrows_col,arrows_col,'black','black','black'),
		       cex=0.6)
		  
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
	 result=list('plot_study'=plot_study,'LnOR'=true_LnOR,'sample'=sample,
	    'LnOR_VAR'=true_VARLnOR,'VAR_LnOR'=var_value,'LnOR_mean'=mean_value,
		 'p'=p_value)
	 invisible(result)
}
