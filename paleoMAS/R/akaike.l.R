akaike.l <-
function(x,y,interval=c(0.15,1,0.05),plot=TRUE,parameters=FALSE)
{
{
	alphas1<-seq(interval[1],interval[2],interval[3])
	alphas<-rep(alphas1,2)
	degrees<-rep(c(1,2),each=length(alphas1))
	param<-cbind(alphas,degrees)
	aic.loess<-matrix(nrow=length(alphas),ncol=3)
	colnames(aic.loess)<-c("alpha","degree","AIC")
	aic.loess[,1]<-alphas
	aic.loess[,2]<-degrees
	l.funct<-as.list(c(1:nrow(aic.loess)))
	l.p<-matrix(nrow=nrow(aic.loess),ncol=7)
	colnames(l.p)<-c("alpha","degree","n","s2","d1",
		"d2","tr(H)")
	l.p[,1:2]<-param
	for(i in 1:nrow(aic.loess)){
		loess(y~x,span=param[i,1],degree=param[i,2])->l.funct[[i]]
		l.funct[[i]]$n->l.p[i,3]
		sum((l.funct[[i]]$residuals)^2)/l.funct[[i]]$n->l.p[i,4]
		l.funct[[i]]$one.delta->l.p[i,5]
		l.funct[[i]]$two.delta->l.p[i,6]
		l.funct[[i]]$trace.hat->l.p[i,7]
		l.p[i,3]*log(l.p[i,4])+l.p[i,3]*((l.p[i,5]/l.p[i,
			6])*(l.p[i,3]+l.p[i,7]))/(((l.p[i,5]^2)/l.p[i,
			6])-2)->aic.loess[i,3]
		}
	if(plot==TRUE){
	layout(c(1,2))
	plot(aic.loess[1:nrow(aic.loess)/2,c(1,3)],type="l",
		main="degree 1")
	plot(aic.loess[((nrow(aic.loess)/2)+1):nrow(aic.loess),
		c(1,3)],type="l",main="degree 2")
		}
	minimum<-aic.loess[which(aic.loess[,
		3]==min(aic.loess[,3])),]
}
if(parameters==TRUE){
	results<-list(l.p,aic.loess,minimum)
	names(results)<-c("l.p","aic.loess","minimum")
	return(results)
	}
else{
	results<-list(aic.loess,minimum)
	names(results)<-c("aic.loess","minimum")
	return(results)
	}	
}

