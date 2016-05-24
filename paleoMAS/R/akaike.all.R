akaike.all <-
function(x,y,interval=c(0.15,1,0.05))
{
{
	minima<-matrix(nrow=ncol(y),ncol=3)
	colnames(minima)<-c("alpha","degree","AIC")
	rownames(minima)<-colnames(y)
	for(i in 1:ncol(y)){
		akaike.l(x,y[,i],interval=interval,plot=FALSE,
			parameters=FALSE)$minimum->minima[i,]
		}
}
return(minima)
}

