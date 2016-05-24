VS.KS <-
function(Trad_data,MC_data,
	 scatter_alpha=0.01,line_alpha=0.001,
    scatter_col="black",line_col="red",
    xlab=NULL,ylab="KS test p_value",title="KS (Kolmogorov-Smirnov) test"){
	 
	 rowNum=nrow(Trad_data)
	 KS=c()
	 for(j in 1:rowNum){
		 temp=ks.test(Trad_data[j,],MC_data[j,])$p.value
		 KS=c(KS,temp)
	}

    index=KS[which(KS<scatter_alpha)]
	 
    plot(index,type="p",main=title,xlab=xlab,ylab=ylab,col=scatter_col,ylim=c(0,scatter_alpha))
	 abline(h=line_alpha,lwd=2,col=line_col)
	 result=list("KS_p"=KS)
	 invisible(result)
}
