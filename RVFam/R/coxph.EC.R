coxph.EC <- function(snp,phen,test.dat,covar,chr,time){


  maf.fun <- function(x){
    if (all(is.na(x))) maf <- NA else {
       x <- na.omit(x)
       maf <- (sum(x==1) + sum(x==2) * 2) / (2*length(x))
       maf <- ifelse(maf<0.5,maf,1-maf)
    }   
    maf
  }

  maf.fun.x <- function(x,sex){
    if (all(is.na(x))) maf <- NA else {
       x <- na.omit(cbind(x,sex))
       maf <- (sum(x[x[,2]==2,1])+sum(x[x[,2]==1,1])/2)/(2*sum(x[,2]==2)+sum(x[,2]==1))
       maf <- ifelse(maf<0.5,maf,1-maf)
    }  
    maf
  }

  cor.snp <- function(y,x){
    return(sd(y)==0 || abs(cor(y,x,use="complete"))>0.99999999 )
  }

  getIlike <- function(modobj){   
    modobj$history$`frailty(famid)`$c.loglik
  }


  na.ind <-  if (missing(covar)) !is.na(snp) & !is.na(test.dat[,phen]) else 
    if (length(covar)>1) !is.na(snp) & !is.na(test.dat[,phen]) & apply(!is.na(test.dat[,covar]),1,all) else
      !is.na(snp) & !is.na(test.dat[,phen]) & !is.na(test.dat[,covar])
  
  snp1=snp[na.ind]; sex=test.dat[na.ind,"sex.in.ped"]
  phen1=test.dat[na.ind,phen]; survival_time=test.dat[na.ind,time]
  
  if (!missing(covar)) x.covar <- as.matrix(test.dat[na.ind,covar]) 
  
  ###########################################
  count <- table(snp1)
  gntps <- names(count)
  count1 <- rep(0,3) 
  count1[as.numeric(gntps)+1] <- count
  if (count1[1]>count1[3]) MAC <- count1[2]+count1[3]*2 else MAC <- count1[2]+2*count1[1]
  ntotal <- sum(count)
  nmiss <- length(snp)-ntotal
  if (chr=="X") maf_ntotal=maf.fun.x(snp1,sex) else maf_ntotal=maf.fun(snp1)
  #########################################  

  if (length(unique(phen1))<=1 ) {
      coxph.out= list(ntotal,nmiss,maf_ntotal, NA, NA, NA,"single category in phen",NA,MAC,count1[1],count1[2],count1[3])
     return(coxph.out)
  }
  
  ###############################################

  if (length(count)==1 ) {
     coxph.out <- list(ntotal,nmiss,maf_ntotal, NA, NA, NA,"single category in SNP",NA,MAC,count1[1],count1[2],count1[3])
     return(coxph.out)
  } 
 
  ###############################################
  if (length(unique(snp1))>1 && !missing(covar)){ 
    if (length(covar)>1) colinear <- apply(x.covar,2,cor.snp,x=snp1) else colinear<-cor.snp(x.covar,snp1)
    if(sum(colinear)>0){ 
      coxph.out= list(ntotal,nmiss,maf_ntotal, NA, NA, NA,"collinarity",NA,MAC,count1[1],count1[2],count1[3])
      return(coxph.out)
    }  
  }

##############################################################################################
#
# Perform COXPH
#
 ################################
  test1<-test.dat[na.ind,]
  if (!missing(covar)){
     coxph.test <- try(coxph(Surv(survival_time,phen1)~snp1+x.covar+frailty(famid),data=test1))
     coxph.test0 <- try(coxph(Surv(survival_time,phen1)~x.covar+frailty(famid),data=test1))
  } else {
     coxph.test <- try(coxph(Surv(survival_time,phen1)~snp1+frailty(famid),data=test1))
     coxph.test0 <- try(coxph(Surv(survival_time,phen1)~frailty(famid),data=test1))
  } 
  if (!"try-error" %in% class(coxph.test0) && !"try-error" %in% class(coxph.test) && (coxph.test$var[1,1]!=0)){ 
     LRT = 2*(getIlike(coxph.test)-getIlike(coxph.test0))
     if (LRT<0) LRT <- 0
     Z <- sign(coxph.test$coef[1])*sqrt(LRT)
     coxph.out <- list(ntotal,nmiss,maf_ntotal,coxph.test$coef[1],sqrt(coxph.test$var[1,1]),Z,NA,pchisq(LRT,df=1,lower.tail=F),MAC,count1[1],count1[2],count1[3])
  } else coxph.out= list(ntotal,nmiss,maf_ntotal, NA, NA, NA,"try-error",NA,MAC,count1[1],count1[2],count1[3])	
  return(coxph.out)
}

