vcrossv.da <-
function(x,f,fold,nsimulat,funct)
{
library(MASS)
{
	ns<-nrow(x)/nsimulat
	leave.out<-matrix(ncol=fold,nrow=floor(ns/fold))
	leave.out[,]<-sample(ns,length(leave.out))	
	leave.out1<-as.list(c(1:nrow(leave.out)))
	for(i in 1:nrow(leave.out)){
		leave.out1[[i]]<-matrix(ncol=2,nrow=ncol(leave.out))
		for(j in 1:ncol(leave.out)){
			leave.out[i,
				j]*nsimulat-nsimulat+1->leave.out1[[i]][j,1]
			leave.out[i,j]*nsimulat->leave.out1[[i]][j,2]
			}
		}
	leave.out3<-matrix(nrow=nrow(leave.out),ncol=fold*nsimulat)
	leave.out2<-as.list(c(1:nrow(leave.out)))
	for(i in 1:nrow(leave.out)){
		leave.out2[[i]]<-matrix(ncol=nsimulat,
			nrow=ncol(leave.out))
		for(j in 1:ncol(leave.out)){
			leave.out1[[i]][j,1]:leave.out1[[i]][j,
				2]->leave.out2[[i]][j,]
			}
		leave.out3[i,]<-as.vector(t(leave.out2[[i]]))
		}
	df.f<-as.list(c(1:nrow(leave.out3)))
	predicted<-as.list(c(1:nrow(leave.out3)))
	predicted1<-matrix(ncol=2,nrow=nsimulat*fold*nrow(leave.out))
	a<-matrix(ncol=2,nrow=nrow(leave.out))
	a[,1]<-seq(1,nsimulat*fold*nrow(leave.out),by=fold*nsimulat)
	a[,2]<-seq(fold*nsimulat,nsimulat*fold*nrow(leave.out),
		by=fold*nsimulat)
	for(i in 1:nrow(leave.out3)){
		funct(x[-leave.out3[i,],],f[-leave.out3[i,]])->df.f[[i]]
		predict(df.f[[i]],
			x[leave.out3[i,],])$posterior->predicted[[i]]
		predicted1[a[i,1]:a[i,2],]<-predicted[[i]]
		}
	rownames(predicted1)<-rownames(x[as.vector(t(leave.out3)),])
	colnames(predicted1)<-colnames(predicted[[1]])
	b<-matrix(ncol=2,nrow=nrow(predicted1)/nsimulat)
	b[,1]<-seq(1,nrow(predicted1),by=nsimulat)
	b[,2]<-seq(nsimulat,nrow(predicted1),by=nsimulat)
	posterior<-b
	rownames(posterior)<-rownames(predicted1)[b[,1]]
	colnames(posterior)<-colnames(predicted1)
	for(i in 1:nrow(b)){
		mean(predicted1[b[i,1]:b[i,2],1])->posterior[i,1]
		mean(predicted1[b[i,1]:b[i,2],2])->posterior[i,2]
		}
	classification<-ifelse(posterior[,]>0.5,1,0)
	fact<-levels(f)
	prior<-matrix(ncol=2,nrow=(nrow(x)/nsimulat))
	taxan<-seq(1,length(f),by=nsimulat)
	colnames(prior)<-fact
	rownames(prior)<-rownames(x)[taxan]
	fact1<-f[taxan]
	prior[,1]<-ifelse(fact1==fact[1],1,0)
	prior[,2]<-ifelse(fact1==fact[2],1,0)
	taxan1<-intersect(rownames(classification),rownames(prior))
	error1<-abs(prior[taxan1,]-classification[taxan1,])
	accuracy<-(1-sum(error1[,1])/nrow(error1))*100
}
results<-list(posterior,classification,accuracy)
names(results)<-c("posterior","classification","accuracy")
return(results)
}

