VS.QQ <-
function(Trad_data,MC_data,scatter_col="black",line_col="black",
    title="QQ plot",xlab="Quantile for TradPerm data)",
	 ylab="Quantile for MCPerm data"){
	 
	 qqplot(Trad_data,MC_data,main=title,col=scatter_col,
	       xlab=xlab,
	       ylab=ylab)
	 abline(a=0,b=1,lwd=1,col=line_col)	
}
