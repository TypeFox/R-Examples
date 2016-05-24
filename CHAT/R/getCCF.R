getCCF <- function(VCFdir,sAGPfile,output="new.vcf",nt=FALSE,nc=10,tc=11,AD=3,filter=TRUE,TCGA=TRUE){
	new.dd=get(load(sAGPfile))
	newfile<-c()
	sID<-unique(rownames(new.dd))
    if(nt){
        tmpList<-mclapply(sID,function(x)getSampleCCF(x,new.dd,VCFdir,nc=nc,tc=tc,filter=filter,AD=AD,TCGA=TCGA))
    }else tmpList<-lapply(sID,function(x)getSampleCCF(x,new.dd,VCFdir,nc=nc,tc=tc,filter=filter,AD=AD,TCGA=TCGA))
	for(tmp in tmpList)newfile<-rbind(newfile,tmp)
	write.table(newfile,file=output,row.names=F,col.names=F,quote=F,sep='\t')
}