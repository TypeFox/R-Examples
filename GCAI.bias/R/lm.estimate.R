lm.estimate <-
function(obj.train,fit.cut.train=5){

	dat.train<-cbind(obj.train$counts,obj.train$index)
	dat.train<-dat.train[which(dat.train[,1]>0&dat.train[,2]>0&rowSums(dat.train[,1:2])>fit.cut.train),]
	dat.train[,1]<-log2(dat.train[,1]/dat.train[,2])
	dat.train<-dat.train[,-2]
	colnames(dat.train)[1]<-"lr"

	colnames(dat.train)[-1]<-paste("X",colnames(dat.train)[-1],sep="")

	vars <- colnames(dat.train)

	f.lm <- as.formula(paste(vars[1], "~", paste(vars[-1], collapse="+"),collapse=""))

	out.lm<-glm(f.lm, data=dat.train)

	coe.lm<-matrix(out.lm$coe,ncol=1)

	coe.lm[which(is.na(coe.lm[,1])),]<-0

	return(coe.lm)
}
