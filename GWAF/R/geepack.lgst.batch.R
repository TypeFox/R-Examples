#dichotomous traits only
#phenfile: phenotype file name in quotation marks,must provide
#genfile" genotype file name in quotation marks, must provide
#outfile: output file name in quotation marks,must provide
#library: path of the library with GEE packge  
#model can be "a", "d", "r","g"
#pedfile: famid id fa mo sex 
#famid is cluster id
geepack.lgst.batch=function(genfile,phenfile,pedfile,outfile,phen,covars=NULL,model="a",col.names=T,sep.ped=",",sep.phe=",",sep.gen=","){
  print(paste("phenotype data = ", phenfile))
  print(paste("genotype data = ", genfile))
  print(paste("pedigree data = ", pedfile))
  print(paste("Result of GEE analyses =",outfile))
  if(is.null(covars)){
    print("Covariates = NONE")
  }else{
    print(paste("Covariate(s) =",covars,collapse=", "))
  }
  print("Running GEE")

  read.in.data <- function(phenfile,genfile,pedfile,sep.ped=sep.ped,sep.phe=sep.phe,sep.gen=sep.gen) {
  print("Reading in Data")
  ped.dat <- read.table(genfile,header=TRUE,na.strings="",sep=sep.gen)
  snp.names <- names(ped.dat)[-1]
  pedigree <- read.table(pedfile,header=TRUE,sep=sep.ped)
  gntp.all <- merge(pedigree,ped.dat,by="id")

#read in phenotype data
  phen.dat=read.table(phenfile,header=TRUE,sep=sep.phe)
  phen.name=colnames(phen.dat)[-1]
  n.snp=length(names(gntp.all))

  if(length(grep("^sex$",colnames(phen.dat)))==0) {
  phensnp.dat<-merge(gntp.all,phen.dat,by=c("id"))
  } else {
## sex is one of the columns in the phenotype file
  phensnp.dat<-merge(gntp.all,phen.dat,by=c("id","sex"))
  }
  print("Done reading in data")
  return(list(data=phensnp.dat,snps=snp.names,phen.name=phen.name))
}

  phensnp.dat <- read.in.data(phenfile,genfile,pedfile,sep.ped=sep.ped,sep.phe=sep.phe,sep.gen=sep.gen)
  snplist<-phensnp.dat$snps
  test.dat<-phensnp.dat$data
  if (!is.null(covars) & sum(snplist %in% covars)>=1) {
     names(test.dat)[which(names(test.dat) %in% paste(snplist[snplist %in% covars],".x",sep=""))] <- snplist[snplist %in% covars]
     covars[covars %in% snplist] <- paste(covars[covars %in% snplist],".y",sep="")
  }
  test.dat<-test.dat[order(test.dat$famid),]
  if (!is.null(covars)) {
     covars.dat <- na.omit(test.dat[,covars])
     single.cov <- F
     if (length(covars)==1) single.cov <- var(covars.dat)==0 else {
        single.cov <- any(apply(covars.dat,2,var)==0)
        if (single.cov) stop(paste("Single category in covariates!"))
        for (i in covars){
            cov1 <- covars.dat[,i]
            if (!is.numeric(cov1)) cov1 <-as.numeric(as.factor(cov1))
            for (j in covars[covars!=i]){
                cov2 <- covars.dat[,j]
                if (!is.numeric(cov2)) cov2 <-as.numeric(as.factor(cov2))
                if (abs(cor(cov1,cov2))>0.99999999) stop(paste("Highly correlated covariates ",i," and ",j,"!!",sep=""))
            }
        }
     }
  }   
  
  #library(geepack)

  final1<-c()
  if (is.null(covars)){
     print(paste("No Covariates, Running:",phen))
     if (length(snplist)<2) { 
        temp.out <- c(phen,snplist,geepack.lgst(snp=test.dat[,snplist],phen=phen,test.dat=test.dat,model=model))
     } else {
        temp.out <-as.data.frame(apply(test.dat[,phensnp.dat$snps],2,geepack.lgst,phen=phen,test.dat=test.dat,model=model))
  	 temp.out <-cbind(rep(phen,ncol(temp.out)),colnames(temp.out),t(temp.out))
       }
     final1 <- rbind(final1,temp.out)
  } else {
     print(paste("Covariates, Running:",phen))
     if (length(snplist)<2) { 
        temp.out <- c(phen,snplist,geepack.lgst(snp=test.dat[,snplist],phen=phen,test.dat=test.dat,covar=covars,model=model))
     } else {
        temp.out <-as.data.frame(apply(test.dat[,phensnp.dat$snps],2,geepack.lgst,phen=phen,test.dat=test.dat,covar=covars,model=model))
  	 temp.out <-cbind(rep(phen,ncol(temp.out)),colnames(temp.out),t(temp.out))
       }
     final1 <- rbind(final1,temp.out)
  }


  if (model %in% c("a","d","r")) {
  	
	colnames(final1)=c("phen","snp","n0","n1","n2","nd0","nd1","nd2","miss.0","miss.1","miss.diff.p","beta",
				"se","chisq","df","model","remark","pval")
  }else
	colnames(final1)=c("phen","snp","n0","n1","n2","nd0","nd1","nd2","miss.0","miss.1","miss.diff.p",
			"beta10","beta20","beta21",
			"se10","se20","se21","chisq","df","model","remark","pval")

 
  write.table(as.matrix(final1),outfile,col.names=col.names, row.names=F,quote=F,sep=",",na="",append=T)
}
