plotCov <-
function(cov)
{
	V	<- cov
	R	<- cov2cor(V)
	V[lower.tri(V)] 				<- 0
	R[upper.tri(R,diag=TRUE)] 		<- 0
	par(mar=c(7,7,4,4))
	image(V+R,main="Correlation (Upper); Covariance (Lower)",xaxt="n",yaxt="n")
	axis(2,at=seq(0,1,1/((dim(V)[1])-1)),labels=rownames(V),las=2,cex.axis=0.6)
	axis(1,at=seq(0,1,1/((dim(V)[1])-1)),labels=rownames(V),las=2,cex.axis=0.6)
	return(V+R)
}
