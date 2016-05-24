PermMeta.LnOR.boxplot <-
function(PermMeta,plot="LnOR",plot_study="all",order="no",
    main="LnOR,no order",true_value_pch=3,pch_col="red",pos=3,text_col="blue",digits=2){
    if (!is.element(order, c("LnOR","LnOR_VAR","VAR_LnOR","sample",'no'))){
	    stop("'order' should be 'LnOR','LnOR_VAR','VAR_LnOR','sample','no'.")
	 }
	 if (!is.element(plot, c("LnOR", "LnOR_VAR"))){
	     stop("'plot' should be 'LnOR', 'LnOR_VAR'.")
	 }
	 if(class(PermMeta)!='PermMeta'){
	    stop("Function 'PermMeta.LnOR.boxplot' applied to an object of class 'PermMeta'.")
	 }
	 if(length(plot_study)==1 && plot_study=="all"){
	    study_num=PermMeta$study_num
		 plot_study=1:study_num
	 }else{
	    study_num=length(plot_study)
	 }
	 
	 switch(plot,
	    LnOR={
		    plot_data=PermMeta$perm_LnOR[plot_study,,drop=FALSE]
			 y_lim=c(range(plot_data)[1],range(plot_data)[2]*1.5)
		 },
		 LnOR_VAR={
		    plot_data=PermMeta$perm_VARLnOR[plot_study,,drop=FALSE]
			 y_lim=c(range(plot_data)[1],range(plot_data)[2]*1.005)
		 }
	 )
    true_LnOR=PermMeta$true_LnOR[plot_study]
	 LnOR_VAR=PermMeta$true_VARLnOR[plot_study]
	 VAR_LnOR=apply(plot_data,1,var)
	 sample_value=PermMeta$sample[plot_study]
	 
	 rownames(plot_data)=paste("study",plot_study,sep=" ")
	 switch(order,
	    LnOR={
		    rank_index=rank(true_LnOR)
		 },
		 LnOR_VAR={
		    rank_index=rank(LnOR_VAR)
		 },
		 VAR_LnOR={
		    rank_index=rank(VAR_LnOR)
		 },
		 sample={
		    rank_index=rank(sample_value)
		 },
		 no={
		    rank_index=1:study_num
		 }
	 )
	 
	 order_index=c()
	 for(i in 1:study_num){
	    order_index=c(order_index,which(rank_index>=i & rank_index<i+1))
	 }
	 plot_data=plot_data[order_index,]
	 true_LnOR=true_LnOR[order_index]
	 LnOR_VAR=LnOR_VAR[order_index]
	 VAR_LnOR=VAR_LnOR[order_index]
	 sample_value=sample_value[order_index]
	 
	 ########################### 
	 #par(oma=c(0, 1, 3, 0))
	 boxplot(t(plot_data),las=2,ylim=y_lim)
	 legend('top',legend=main,cex=0.7,inset=-0.1,xpd=TRUE)
	 #mtext(text=main,side=3,outer=TRUE)
	 switch(plot,
	    LnOR={
		    points(1:study_num,true_LnOR,col=pch_col,pch=true_value_pch)
			 text(1:study_num,true_LnOR,labels=paste(round(true_LnOR,digits)," ",sep=""),
	          srt=0,pos=pos,cex=0.5,col=pch_col)	
		 },
		 LnOR_VAR={
		    points(1:study_num,LnOR_VAR,col=pch_col,pch=true_value_pch)
			 text(1:study_num,LnOR_VAR,labels=paste(round(LnOR_VAR,digits)," ",sep=""),
	          srt=0,pos=pos,cex=0.5,col=pch_col)
		 }
	 )
	 
	 # text(1:study_num,rep(-0.5,study_num),labels=paste(round(VAR_LnOR,digits)," ",sep=""),
	    # srt=30,pos=pos,cex=0.7,col=text_col)
	 # text(1:study_num,rep(-0.5,study_num),labels=paste(sample_value," ",sep=""),
	    # srt=30,pos=pos,cex=0.7,col=text_col)
	 
	 ## return value
	 result=list('plot'=plot,'plot_num'=study_num,'plot_order'=order,'order_index'=order_index,
	    'true_LnOR'=true_LnOR,'LnOR_VAR'=LnOR_VAR,'VAR_LnOR'=VAR_LnOR,'sample'=sample_value)
	 invisible(result)
}
