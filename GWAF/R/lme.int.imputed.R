lmepack.int.batch.imputed <- function(phenfile,genfile,pedfile,phen,kinmat,covars,cov.int,sub="N",outfile,col.names=T,sep.ped=",",sep.phe=",",sep.gen=","){

###########################################################
  #library(coxme)
  kmat <- NULL; rm(kmat)
  #check the existence of kinship matrix
  trykin<-try(load(kinmat))
  if (inherits(trykin,"try-error"))
        stop(paste('kinship matrix does not exist at ',kinmat))

  #do not run when 1) no covariates or no cov.int; 2) the length of cov.int is not 1; 3) cov.int is not in covariates 
  if (sum((!is.na(cov.int))*(sum(is.na(covars))!=1))!=1 | length(cov.int)!=1 | sum(cov.int %in% covars)!=1) stop('one interaction covariate is allowed and it has to be in covariates')
  
  cor.snp <- function(y,x){  #########
   if(!is.numeric(y))y<-as.numeric(as.factor(y)) ########
   return(sd(y)==0 || abs(cor(y,x,use="complete"))>0.99999999 )} #########

  if (sum(is.na(sub))==1) sub <- "N"  ####061009
  assign("phen", phen, pos=-1,inherits=T)
  assign("cov.int", cov.int, pos=-1,inherits=T)

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

#####################main programs##########################
 
  phensnp.dat <- read.in.data(phenfile,genfile,pedfile,sep.ped=sep.ped,sep.phe=sep.phe,sep.gen=sep.gen)
  snplist<-phensnp.dat$snps
  if (sum(phensnp.dat$phen.name %in% covars)==length(covars)) phenlist<-phensnp.dat$phen.name[!phensnp.dat$phen.name %in% covars] else  #######
     stop('some covariates are not available') #######

  test.dat <- phensnp.dat$data
  if (length(table(test.dat[,cov.int]))!=2 & sub=="Y") stop('No subset analysis for non-binary interaction covariate!') ####061009
  result <- NULL

  if (sum(is.na(covars))==0 & sum(snplist %in% covars)>=1) {
     names(test.dat)[which(snplist %in% covars)+6] <- snplist[snplist %in% covars]
     covars[covars %in% snplist] <- paste(covars[covars %in% snplist],".y",sep="")
  }

  cov.int.snp <- NA #100708
  if (sum(is.na(cov.int))==0 & sum(snplist %in% cov.int)==1) {
      cov.int.snp <- snplist[snplist %in% cov.int]
      cov.int <- paste(cov.int,".y",sep="")
  } #100708

  idlab <- "id"

      for (i in snplist) {
          test2.dat <- na.omit(test.dat[,c(i,phen,idlab,covars)]) #######
          assign("test2.dat", test2.dat, pos=-1,inherits=T)
          x.covar<-as.matrix(test2.dat[,covars])     
          assign("x.covar", x.covar, pos=-1,inherits=T)
          id <- test2.dat[,idlab]  ####### 5)
          assign("id", id, pos=-1,inherits=T)
          bin.flag <- length(table(test2.dat[,cov.int])) 
          snp <- test2.dat[,i]
          assign("snp", snp, pos=-1,inherits=T)
          if (bin.flag==2) {
             bin <- sort(as.numeric(names(table(test2.dat[,cov.int])))) 
             assign("bin", bin, pos=-1,inherits=T)
          } 

          if (length(covars)>1) colinear <- apply(x.covar,2,cor.snp,x=test2.dat[,i]) else colinear <- cor.snp(x.covar,test2.dat[,i]) #######
          x.int<-test2.dat[,cov.int]*snp;   #######071409
          assign("x.int", x.int, pos=-1,inherits=T)
          colinear <- c(colinear,cor.snp(snp,x.int),cor.snp(test2.dat[,cov.int],x.int))  

          imaf <- mean(test2.dat[,i])/2
          count<-table(snp)
          count1 <- c(0,0)
          count1[1] <- length(test2.dat[,i])
          count1[2] <- imaf
          testiMAF <- 2*count1[1]*imaf*(1-imaf)
			 		
          if (bin.flag<2) result<-rbind(result,c(phen,i,cov.int,count1[1:2],rep(NA,8))) else { 
             if (length(count)==1 | sum(colinear,na.rm=T)>0 | length(unique(snp))==1 || testiMAF <= 1 | length(unique(x.int))==1){
	          if (bin.flag>2) result<-rbind(result,c(phen,i,cov.int,count1[1:2],rep(NA,8))) else {
                    if (sub=="Y") result<-rbind(result,c(phen,i,cov.int,count1[1:2],rep(NA,13))) else result<-rbind(result,c(phen,i,cov.int,count1[1:2],rep(NA,8)))} 			    
		   } else {
                    if (bin.flag>2) {lme.out<-try(lmekin(test2.dat[,phen]~snp+x.int+x.covar+(1|id),varlist=kmat,na.action=na.omit)) 
                       if (!class(lme.out)%in%"try-error"){
                          tmp<-c(lme.out$var[2,3],"additive",lme.out$coef$fixed[2],sqrt(lme.out$var[2,2]),pchisq(lme.out$coef$fixed[2]^2/lme.out$var[2,2],1,lower.tail=F),
                               lme.out$coef$fixed[3],sqrt(lme.out$var[3,3]),pchisq(lme.out$coef$fixed[3]^2/lme.out$var[3,3],1,lower.tail=F))     
                       } else tmp <- rep(NA,8)
                    } else {
                         if (sub=="Y") { #061009
                            x.covar1 <- as.matrix(test2.dat[,covars[covars!=cov.int]])         
                            assign("x.covar1", x.covar1, pos=-1,inherits=T)
                            #tab1 <- table(snp[x.covar[,cov.int]==bin[1]]) ######05132009 imputed.X
                            #tab2 <- table(snp[x.covar[,cov.int]==bin[2]])
                            #if (length(tab1)==1 | length(tab2)==1) result<-rbind(result,c(phen,i,cov.int,count1[1:2],rep(NA,13))) else { ######05132009 imputed.X                                            
                               lme.out1<-try(lmekin(test2.dat[,phen]~snp+x.covar1+(1|id),varlist=kmat,na.action=na.omit,subset=test2.dat[,cov.int]==bin[1]))
                               lme.out2<-try(lmekin(test2.dat[,phen]~snp+x.covar1+(1|id),varlist=kmat,na.action=na.omit,subset=test2.dat[,cov.int]==bin[2]))
                               lme.out<-try(lmekin(test2.dat[,phen]~snp+x.int+x.covar+(1|id),varlist=kmat,na.action=na.omit)) 
                               if (!class(lme.out)%in%"try-error" & !class(lme.out1)%in%"try-error" & !class(lme.out2)%in%"try-error"){
                                  tmp<-c("additive",lme.out$coef$fixed[2],sqrt(lme.out$var[2,2]),pchisq(lme.out$coef$fixed[2]^2/lme.out$var[2,2],1,lower.tail=F),
                                       lme.out1$coef$fixed[2],sqrt(lme.out1$var[2,2]),pchisq(lme.out1$coef$fixed[2]^2/lme.out1$var[2,2],1,lower.tail=F),
                                       lme.out2$coef$fixed[2],sqrt(lme.out2$var[2,2]),pchisq(lme.out2$coef$fixed[2]^2/lme.out2$var[2,2],1,lower.tail=F),
                                       lme.out$coef$fixed[3],sqrt(lme.out$var[3,3]),pchisq(lme.out$coef$fixed[3]^2/lme.out$var[3,3],1,lower.tail=F)) #090808
                               } else tmp <- rep(NA,13)
                            #}                         
                            } else {lme.out<-try(lmekin(test2.dat[,phen]~snp+x.int+x.covar+(1|id),varlist=kmat,na.action=na.omit)) #061009
                            if (!class(lme.out)%in%"try-error"){
                            tmp<-c(lme.out$var[2,3],"additive",lme.out$coef$fixed[2],sqrt(lme.out$var[2,2]),pchisq(lme.out$coef$fixed[2]^2/lme.out$var[2,2],1,lower.tail=F),
                               lme.out$coef$fixed[3],sqrt(lme.out$var[3,3]),pchisq(lme.out$coef$fixed[3]^2/lme.out$var[3,3],1,lower.tail=F))  
                            } else tmp <- rep(NA,8)}                          
   		            }
              result<-rbind(result,c(phen,i,cov.int,count1[1:2],tmp))
		}
	   }
      }

if (length(table(test.dat[,cov.int]))>2) colnames(result)<-c("phen","snp","covar_int","n","AF","cov_beta_snp_beta_int","model","beta_snp","se_snp","pval_snp","beta_int","se_int","pval_int")  else {
   if (sub=="Y") colnames(result)<-c("phen","snp","covar_int","n","AF","model","beta_snp","se_snp","pval_snp","beta_snp_cov0","se_snp_cov0","pval_snp_cov0","beta_snp_cov1","se_snp_cov1","pval_snp_cov1","beta_int","se_int","pval_int") else  
       colnames(result)<-c("phen","snp","covar_int","n","AF","cov_beta_snp_beta_int","model","beta_snp","se_snp","pval_snp","beta_int","se_int","pval_int")}
  
if (sum(is.na(cov.int.snp))==0 & length(cov.int.snp)==1) { 
     result[,"covar_int"] <- cov.int.snp
} 

write.table(result, outfile, quote=F,row.names=F, col.names=T,sep=",",na="",append=T)

}