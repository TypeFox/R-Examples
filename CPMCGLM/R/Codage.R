Codage<-function(data,quantile,var.cod)
 {
	if (class(quantile)=="NULL")
	{ 
		X1<-NULL
		coup2<-NULL
	}else{
		X1<-matrix(nrow=length(var.cod),ncol=nrow(quantile))
	#A partir de chaque point de coupure defini dans quantile on d�finie une variable binaire
		for (i in 1:nrow(quantile)){
			X1[,i]<- cut(var.cod,breaks=c(min(var.cod),quantile(var.cod,quantile[i,1:(ncol(quantile))]),max(var.cod)),include.lowest = TRUE)
		}
		coup2<-NULL
		for (i in 1:nrow(quantile)){
			coup2[i]<-unique(max(X1[,i]))-1
		}
		X1<-as.matrix(X1)-1
	}
	res<-list(X1,coup2)
 }