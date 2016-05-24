codcut<-function(data,cutpoint,var.cod)
 {

	coup3<-NULL
	
	if (missing(cutpoint))
	{
		X1<-NULL
	}else{
		X1<-matrix(nrow=length(var.cod),ncol=nrow(cutpoint))
	#A partir de chaque point de coupure defini dans quantile on d�finie une variable binaire
		for (i in 1:nrow(cutpoint))
		{
			X1[,i]<- cut(var.cod,breaks=c(min(var.cod),cutpoint[i,1:(ncol(cutpoint))],max(var.cod)),include.lowest = TRUE)
		}
		
		for (i in 1:nrow(cutpoint))
		{
			coup3[i]<-unique(max(X1[,i]))-1
		}
		X1<-as.matrix(X1)-1
	}
	res<-list(X1,coup3)
 }