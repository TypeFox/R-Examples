boscat.f <-
function(data,nv,col,mb,crit="f2")
{
	## nv=nodevec,mb=minbucket
	##anova(lm(data),lm(newdata))$F Value[2]
	##(y2rsq-y1rsq)*(n1-n2)/(1-y2rsq)
	critn<--1
	if(crit=="F-value") {critn<-2}
	if(crit=="R2change") {critn<-1}
	if(crit=="f2") {critn<-0}
	x<-col[nv==1]
	y<-resid.f(data)
	y<-y[nv==1]
	meanx<-data.frame(sort(tapply(y,x,mean)))
	ordernames<-row.names(meanx)
	trx<-as.numeric(factor(x,levels=ordernames))
	trx2<-as.numeric(factor(col,levels=ordernames))
	trx[is.na(trx)]<-0
	trx2[is.na(trx2)]<-0
	m<-dim(data)[1]
	n<-dim(data)[2]
	mm<-length(trx)
	res<-numeric(2)
  	res<-.Fortran("rs_boscat",res=as.double(res),data=as.matrix(data),
	     as.integer(nv),as.integer(m),as.integer(n),
	     as.double(trx),as.double(trx2),as.integer(mm),
	     as.integer(mb),as.integer(critn),package="stima")$res
	trx2[trx2==0]<-NA
	object<-list(res=res,ordernames=ordernames,trx=trx2)
}
