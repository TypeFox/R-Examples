fossil.values <-
function(distan,age,quant=0.05,detrend=FALSE,crossv)
{
{
	values.par<-matrix(nrow=length(age),ncol=2)
	values.par[,1]<-age
	colnames(values.par)<-c("Age","Estimated")
	for(i in 2:(length(age)+1)){
		mean(distan[1,(which(distan[i,]<quantile(distan[i,
			],quant)))])->values.par[i-1,2]
		}	
	if(detrend==TRUE){
		coef<-crossv$coef
		values.par<-cbind(values.par,c(1:length(age)))
		colnames(values.par)<-c("Age","Estimated",
			"Detrended estimated")
		values.par[,3]<-values.par[,2]+coef[1,
			1]+coef[1,2]*values.par[,2]
		values.par[,3]<-values.par[,3]-crossv$transl
		}
}
return(values.par)
}

