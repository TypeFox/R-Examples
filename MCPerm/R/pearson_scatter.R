pearson_scatter <-
function(Trad_data,MC_data,scatter_col="gray28",line_col="black",
    title=NULL,xlab="TradPerm P-value",ylab="MCPerm P-value"){
	 temp=cor.test(Trad_data,MC_data,alternative="two.sided",method="pearson")
	 stat=round(temp$estimate[[1]],3)
	 r=paste("r = ",stat,sep="")
    pvalue=round(temp$p.value,4)
	 p=paste("p = ",pvalue,sep="")
	 plot(Trad_data,MC_data,type="p",main=title,col=scatter_col,xlab=xlab,ylab=ylab)
    legend("topleft",legend=c(r,p),cex=0.7)
    abline(a=0,b=1,lwd=2) 
}
