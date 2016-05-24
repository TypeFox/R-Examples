#dichotomous traits only
#phenfile: phenotype file name in quotation marks,must provide
#genfile" genotype file name in quotation marks, must provide
#outfile: output file name in quotation marks,must provide
#library: path of the library with GEE packge  
#pedfile: famid id fa mo sex
#famid is cluster id
geepack.lgst.int.batch=function(genfile,phenfile,pedfile,outfile,phen,covars,cov.int,sub="N",col.names=T,sep.ped=",",sep.phe=",",sep.gen=","){
  print(paste("phenotype data = ", phenfile))
  print(paste("genotype data = ", genfile))
  print(paste("pedigree data = ", pedfile))
  print(paste("Result of GEE analyses =",outfile))

  if (missing(covars) | missing(cov.int)| length(cov.int)!=1 | sum(cov.int %in% covars)!=1) stop('no covariates or no covariate for interaction or other covariate issue')

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
  if (sum(phensnp.dat$phen.name %in% covars)==length(covars)) phenlist<-phensnp.dat$phen.name[!phensnp.dat$phen.name %in% covars] else  
     stop('some covariates are not available') 
 
  test.dat<-phensnp.dat$data
  if (length(table(test.dat[,cov.int]))!=2 & !is.na(sub) & sub=="Y") stop('No subset analysis for non-binary interaction covariate!') ####061209  

  test.dat<-test.dat[order(test.dat$famid),]
  if (sum(is.na(covars))==0 & sum(snplist %in% covars)>=1) { 
     names(test.dat)[which(snplist %in% covars)+6] <- snplist[snplist %in% covars] 
     covars[covars %in% snplist] <- paste(covars[covars %in% snplist],".y",sep="") 
  } 

  cov.int.snp <- NA 
  if (sum(is.na(cov.int))==0 & sum(snplist %in% cov.int)==1) {
      cov.int.snp <- snplist[snplist %in% cov.int]
      cov.int <- paste(cov.int,".y",sep="")
  } 

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

  #library(geepack)
  final1<-c()
  print(paste("Covariates, Running:",phen))
  if (length(snplist)<2) { 
     temp.out <- c(phen,snplist,geepack.lgst.int(snp=test.dat[,snplist],phen=phen,test.dat=test.dat,covar=covars,cov.int=cov.int,sub=sub))
     final1 <- rbind(final1,temp.out)
  } else {
     for (i in 1:length(snplist)) {
	  temp.out <-c(phen,snplist[i],geepack.lgst.int(snp=test.dat[,snplist[i]],phen=phen,test.dat=test.dat,covar=covars,cov.int=cov.int,sub=sub))
         final1 <- rbind(final1,temp.out)  
     }
  }


if (ncol(final1)==16) {    
 colnames(final1)<-c("phen","snp","covar_int","n","AF","nd","AFd","cov_beta_snp_beta_int","model","beta_snp","se_snp","pval_snp","beta_int","se_int","pval_int","remark") 
    final1 <- final1[,c(1:14,16,15)] } else {    
    colnames(final1)<-c("phen","snp","covar_int","n","AF","nd","AFd","model","beta_snp","se_snp","pval_snp","beta_snp_cov0",
			"se_snp_cov0","pval_snp_cov0","beta_snp_cov1","se_snp_cov1","pval_snp_cov1","beta_int","se_int","pval_int","remark")  
    final1 <- final1[,c(1:19,21,20)]
}

if (sum(is.na(cov.int.snp))==0 & length(cov.int.snp)==1) { 
     final1[,"covar_int"] <- cov.int.snp
} 
 
write.table(as.matrix(final1),outfile,col.names=T, row.names=F,quote=F,sep=",",na="",append=T)
}

