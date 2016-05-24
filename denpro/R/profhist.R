profhist<-function(dendat,binlkm,cvol=TRUE,ccen=TRUE,cfre=FALSE)
{
#esim. dendat<-matrix(rnorm(20),10) on 10*2 matriisi

epsi<-0
hi<-histo(dendat,binlkm,epsi)
recs<-hi$recs
hisfrekv<-hi$values

pr<-profgene(values=hisfrekv,recs=recs,frekv=hisfrekv,cvol=cvol,ccen=ccen,
cfre=cfre)

return(list(parent=pr$parent,level=pr$level,invalue=pr$invalue,
volume=pr$volum,center=pr$center,nodefrek=pr$nodefrek,recs=recs,
hisfrekv=t(hisfrekv),lsets=pr$lsets))
}










