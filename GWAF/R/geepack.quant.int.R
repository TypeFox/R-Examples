geepack.quant.int.batch <- function(phenfile,genfile,pedfile,phen,covars,cov.int,sub="N",outfile,col.names=T,sep.ped=",",sep.phe=",",sep.gen=","){

######################check input files, parameters#########
###########################################################
  
  #do not run when 1) no covariates or no cov.int; 2) the length of cov.int is not 1; 3) cov.int is not in covariates 
  if (sum((!is.na(cov.int))*(sum(is.na(covars))!=1))!=1 | length(cov.int)!=1 | sum(cov.int %in% covars)!=1) stop('one interaction covariate is allowed and it has to be in covariates')

  cor.snp <- function(y,x){  #########
  if(!is.numeric(y))y<-as.numeric(as.factor(y)) ########
  return(sd(y)==0 || abs(cor(y,x,use="complete"))>0.99999999 )} #########

  if (sum(is.na(sub))==1) sub <- "N"  ####061509
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
  test.dat <- test.dat[order(test.dat$famid),]
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

  idlab="famid"

      for (i in snplist) {
          test2.dat <- na.omit(test.dat[,c(i,phen,idlab,covars)]) #####
          x.covar<-as.matrix(test2.dat[,covars])  
          famid <- test2.dat[,idlab] 
          bin.flag <- length(table(test2.dat[,cov.int])) 
          if (bin.flag==2) bin <- sort(as.numeric(names(table(test2.dat[,cov.int])))) 
          maf <- mean(test2.dat[,i],na.rm=T)/2 
          n <- nrow(test2.dat) 
          snp <- test2.dat[,i]
          count<-table(snp)         
          gntps<-names(count)
          count1<-rep(0,3) 
          count1[as.numeric(gntps)+1]<-count
 
          if (length(covars)>1 & length(count)>1) colinear <- apply(x.covar,2,cor.snp,x=test2.dat[,i]) else #######
             if (length(covars)==1 & length(count)>1) colinear <- cor.snp(x.covar,test2.dat[,i]) else #######
                if (length(count)==1) colinear <- F
           
          x.int<-test2.dat[,cov.int]*snp;   #######071409
          colinear <- c(colinear,cor.snp(snp,x.int),cor.snp(test2.dat[,cov.int],x.int))  #######071409

          if (bin.flag<2) { ###061509 sub
             if (sub=="Y") result<-rbind(result,c(phen,i,cov.int,n,maf,rep(NA,13))) else result<-rbind(result,c(phen,i,cov.int,n,maf,rep(NA,8)))
          } else { ####if statemenet for bin.flag #090808        
          #if (maf<0.01 | sort(count1)[1]+sort(count1)[2]<10 | sum(colinear,na.rm=T)>0 | length(table(x.int))==1) { 
          if (sort(count1)[1]+sort(count1)[2]<1 | sum(colinear,na.rm=T)>0 | length(table(x.int))==1) {
             if (bin.flag>2) result<-rbind(result,c(phen,i,cov.int,n,maf,rep(NA,8))) else {
                 if (sub=="Y") result<-rbind(result,c(phen,i,cov.int,n,maf,rep(NA,13))) else result<-rbind(result,c(phen,i,cov.int,n,maf,rep(NA,8)))
             }
          } else { ####bin.flag>=2 and good SNP
          mod.lab <- "additive"; 
          if (bin.flag>2) {
                   gee.test.v <- try(geese(test2.dat[,phen]~snp+x.int+x.covar,id=famid,corstr="independence"))
                   if (!"try-error" %in% class(gee.test.v) & gee.test.v$error==0) {           
                      gee.test <- summary(gee.test.v)        
                      tmp <- c(gee.test.v$vbeta[names(gee.test.v$beta)=="snp",names(gee.test.v$beta)=="x.int"],"additive",gee.test$mean["snp",1:2],pchisq(gee.test$mean["snp",3],1,lower.tail=F),
                              gee.test$mean["x.int",1:2],pchisq(gee.test$mean["x.int",3],1,lower.tail=F))
                   } else tmp <- c(rep(NA,8))  
          } else { ###bin.flag==2          
          if (sub=="Y") {
             #if (sum(table(test2.dat[,cov.int],snp)[1,-1])==0 | sum(table(test2.dat[,cov.int],snp)[2,-1])==0) tmp<-rep(NA,13) else {
                      if (sum(!covars %in% cov.int)>0) { ####replace !covar %in% cov.int 
                         xcovar.bin<-as.matrix(test2.dat[,covars[!covars %in% cov.int]])
                         gee.test1 <- try(summary(geese(test2.dat[,phen]~snp+xcovar.bin,id=famid,corstr="independence",subset=test2.dat[,cov.int]==bin[1])))
                         gee.test2 <- try(summary(geese(test2.dat[,phen]~snp+xcovar.bin,id=famid,corstr="independence",subset=test2.dat[,cov.int]==bin[2])))
                         gee.test <- try(summary(geese(test2.dat[,phen]~snp+x.int+x.covar,id=famid,corstr="independence")))
                      } else {
                         gee.test1 <- try(summary(geese(test2.dat[,phen]~snp,id=famid,corstr="independence",subset=test2.dat[,cov.int]==bin[1])))
                         gee.test2 <- try(summary(geese(test2.dat[,phen]~snp,id=famid,corstr="independence",subset=test2.dat[,cov.int]==bin[2])))
                         gee.test <- try(summary(geese(test2.dat[,phen]~snp+x.int+x.covar,id=famid,corstr="independence")))
                      }
                      if (sum(c(class(gee.test1),class(gee.test2),class(gee.test))=="try-error")==0 & sum(c(gee.test$error,gee.test1$error,gee.test2$error))==0) { ###061209
                         tmp <- c("additive",gee.test$mean["snp",1:2],pchisq(gee.test$mean["snp",3],1,lower.tail=F),gee.test1$mean["snp",1:2],
                                 pchisq(gee.test1$mean["snp",3],1,lower.tail=F),gee.test2$mean["snp",1:2],pchisq(gee.test2$mean["snp",3],1,lower.tail=F),
                                 gee.test$mean["x.int",1:2],pchisq(gee.test$mean["x.int",3],1,lower.tail=F))
                  } else {
                         if ("try-error" %in% c(class(gee.test1),class(gee.test2)) & !"try-error" %in% class(gee.test) & gee.test$error==0 & sum(c(gee.test1$error,gee.test2$error))!=0) { 
                            tmp <- c("additive",gee.test$mean["snp",1:2],pchisq(gee.test$mean["snp",3],1,lower.tail=F),rep(NA,6),
                                    gee.test$mean["x.int",1:2],pchisq(gee.test$mean["x.int",3],1,lower.tail=F))
                        } else tmp <- c(rep(NA,13))
                    }
             #         }      
          } else { 
                      gee.test.v <- try(geese(test2.dat[,phen]~snp+x.int+x.covar,id=famid,corstr="independence"))
                      if (!"try-error" %in% class(gee.test.v) & gee.test.v$error==0) {
                         gee.test <- summary(gee.test.v)        
                         tmp <- c(gee.test.v$vbeta[names(gee.test.v$beta)=="snp",names(gee.test.v$beta)=="x.int"],"additive",gee.test$mean["snp",1:2],
                                 pchisq(gee.test$mean["snp",3],1,lower.tail=F),gee.test$mean["x.int",1:2],pchisq(gee.test$mean["x.int",3],1,lower.tail=F))
                      } else tmp <- c(rep(NA,8))
                      }          
                 }
            result<-rbind(result, c(phen,i,cov.int,n,maf,unlist(tmp))) 
          } #######bin.flag>=2 and good SNP
          } ###if statemenet for bin.flag #090808
      }    

if (length(table(test.dat[,cov.int]))>2) colnames(result)<-c("phen","snp","covar_int","n","AF","cov_beta_snp_beta_int","model","beta_snp","se_snp","pval_snp","beta_int","se_int","pval_int") else {
   if (sub=="Y") colnames(result)<-c("phen","snp","covar_int","n","AF","model","beta_snp","se_snp","pval_snp","beta_snp_cov0","se_snp_cov0","pval_snp_cov0","beta_snp_cov1","se_snp_cov1","pval_snp_cov1","beta_int","se_int","pval_int") 
else colnames(result)<-c("phen","snp","covar_int","n","AF","cov_beta_snp_beta_int","model","beta_snp","se_snp","pval_snp","beta_int","se_int","pval_int")
}


if (sum(is.na(cov.int.snp))==0 & length(cov.int.snp)==1) { #100708
     result[,"covar_int"] <- cov.int.snp
} #100708

write.table(result, outfile, quote=F,row.names=F, col.names=T,sep=",",na="",append=T)

}