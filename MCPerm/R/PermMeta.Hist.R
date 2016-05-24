PermMeta.Hist <-
function(PermMeta,plot="Qp",fill_col=NULL,border_col='black',
    arrows_col='red',main="Hist plot for heterogeneity Q p_vlaue",
	 xlab="Q p_value",ylab="Density",digits=3){
	 if(class(PermMeta)!='PermMeta'){
	    stop("'PermMeta.Hist' applied to an object of class 'PermMeta'.")
	 }
	 if (!is.element(plot, c("Qp", "I2", "merged_LnOR", "merged_LnOR_VAR", "merged_LnOR_p"))){
	     stop("'plot' should be 'Qp', 'I2', 'merged_LnOR', 'merged_LnOR_VAR', 'merged_LnOR_p'.")
	 }
	 switch(plot,
	    Qp={plot_data=PermMeta$perm_Qp
		    xlim=range(c(plot_data,PermMeta$corrected_result[2,1]))},
		 I2={plot_data=PermMeta$perm_I2
		    xlim=range(c(plot_data,PermMeta$corrected_result[1,2]))},
		 merged_LnOR={plot_data=PermMeta$perm_merged_LnOR
		    xlim=range(c(plot_data,PermMeta$corrected_result[1,3]))},
		 merged_LnOR_VAR={plot_data=PermMeta$perm_merged_VARLnOR
		    xlim=range(c(plot_data,PermMeta$true_merged_LnOR_VAR))},
		 merged_LnOR_p={plot_data=PermMeta$perm_p
		    xlim=range(c(plot_data,PermMeta$corrected_result[2,3]))}
	 )
	 dens=hist(plot_data,freq=FALSE,main=main,col=fill_col,border=border_col,
	    xlim=xlim,xlab=xlab,ylab=ylab)$density
	 y_value=max(dens)/10

	 switch(plot,
	    Qp={
		    temp1=paste('Q_stat',round(PermMeta$corrected_result[1,1],digits),sep='=')
			 temp2=paste('Q_p',round(PermMeta$corrected_result[2,1],digits),sep='=')
			 temp3=paste('p.corrected',round(PermMeta$corrected_result[3,1],digits),sep='=')
			 temp=c(temp1,temp2,temp3)
			 
	       points(PermMeta$corrected_result[2,1],y=-y_value/5,
			    pch=2,col=arrows_col,cex=0.6,lwd=1.3)
	       arrows(x0=PermMeta$corrected_result[2,1],y0=y_value,
			    x1=PermMeta$corrected_result[2,1],y1=0,
	          col=arrows_col,length=0.05,lwd=1.3)
		    legend('top',legend=temp,pch=c(20,20,20),
			    text.col=c('black',arrows_col,'black'),
				 cex=0.6,inset=-0.06,xpd=TRUE,ncol=3)
		},
		 I2={
		    temp1=paste('I2_stat',round(PermMeta$corrected_result[1,2],digits),sep='=')
			 temp2=paste('p.corrected',round(PermMeta$corrected_result[3,2],digits),sep='=')
			 temp=c(temp1,temp2)
			 
	       points(PermMeta$corrected_result[1,2],y=-y_value/5,
			    pch=2,col=arrows_col,cex=0.6,lwd=1.3)
	       arrows(x0=PermMeta$corrected_result[1,2],y0=y_value,
			    x1=PermMeta$corrected_result[1,2],y1=0,
	          col=arrows_col,length=0.05,lwd=1.3)
		    legend('top',legend=temp,pch=c(20,20),
			    text.col=c(arrows_col,'black'),
				 cex=0.6,inset=-0.06,xpd=TRUE,ncol=2)
		},
		 merged_LnOR={
		    temp1=paste('merged_LnOR',round(PermMeta$corrected_result[1,3],digits),sep='=')
			 temp2=paste('merged_LnOR_VAR',round(PermMeta$true_merged_LnOR_VAR,3),sep='=')
			 temp3=paste('merged_LnOR_p',round(PermMeta$corrected_result[2,3],3),sep='=')
			 temp4=paste('p.corrected',round(PermMeta$corrected_result[3,3],3),sep='=')
			 temp=c(temp1,temp2,temp3,temp4)
			 
	       points(PermMeta$corrected_result[1,3],y=-y_value/5,
			    pch=2,col=arrows_col,cex=0.6,lwd=1.3)
	       arrows(x0=PermMeta$corrected_result[1,3],y0=y_value,
			    x1=PermMeta$corrected_result[1,3],y1=0,
	          col=arrows_col,length=0.05,lwd=1.3)
		    legend('top',legend=temp,pch=c(20,20,20,20),
			    text.col=c(arrows_col,'black','black','black'),
				 cex=0.6,inset=-0.06,xpd=TRUE,ncol=2)
		},
		 merged_LnOR_VAR={
		    temp1=paste('merged_LnOR',round(PermMeta$corrected_result[1,3],digits),sep='=')
			 temp2=paste('merged_LnOR_VAR',round(PermMeta$true_merged_LnOR_VAR,digits),sep='=')
			 temp3=paste('merged_LnOR_p',round(PermMeta$corrected_result[2,3],digits),sep='=')
			 temp4=paste('p.corrected',round(PermMeta$corrected_result[3,3],digits),sep='=')
			 temp=c(temp1,temp2,temp3,temp4)
			 
	       points(PermMeta$true_merged_LnOR_VAR,y=-y_value/5,
			    pch=2,col=arrows_col,cex=0.6,lwd=1.3)
	       arrows(x0=PermMeta$true_merged_LnOR_VAR,y0=y_value,
			    x1=PermMeta$true_merged_LnOR_VAR,y1=0,
	          col=arrows_col,length=0.05,lwd=1.3)
		    legend('top',legend=temp,pch=c(20,20,20,20),
			    text.col=c('black',arrows_col,'black','black'),
				 cex=0.6,inset=-0.06,xpd=TRUE,ncol=2)
		},
		
		 merged_LnOR_p={
		    temp1=paste('merged_LnOR',round(PermMeta$corrected_result[1,3],digits),sep='=')
			 temp2=paste('merged_LnOR_VAR',round(PermMeta$true_merged_LnOR_VAR,digits),sep='=')
			 temp3=paste('merged_LnOR_p',round(PermMeta$corrected_result[2,3],digits),sep='=')
			 temp4=paste('p.corrected',round(PermMeta$corrected_result[3,3],digits),sep='=')
			 temp=c(temp1,temp2,temp3,temp4)
			 
	       points(PermMeta$corrected_result[2,3],y=-y_value/5,
			    pch=2,col=arrows_col,cex=0.6,lwd=1.3)
	       arrows(x0=PermMeta$corrected_result[2,3],y0=y_value,
			    x1=PermMeta$corrected_result[2,3],y1=0,
	          col=arrows_col,length=0.05,lwd=1.3)
		    legend('top',legend=temp,pch=c(20,20,20,20),
			    text.col=c('black','black',arrows_col,'black'),
				 cex=0.6,inset=-0.06,xpd=TRUE,ncol=2)
		}
	 )
}
