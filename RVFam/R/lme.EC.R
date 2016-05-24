lme.EC <- function(snp,phen,test.dat,covar,kmat=kmat,chr){

  maf.fun <- function(x){
    if (all(is.na(x))) maf <- NA else {
       x <- na.omit(x)
       maf <- (sum(x==1)+sum(x==2)*2)/(2*length(x))
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

  na.ind <-  if (missing(covar)) !is.na(snp) & !is.na(test.dat[,phen]) else 
              if (length(covar)>1) !is.na(snp) & !is.na(test.dat[,phen]) & apply(!is.na(test.dat[,covar]),1,all) else
                 !is.na(snp) & !is.na(test.dat[,phen]) & !is.na(test.dat[,covar])

  snp1=snp[na.ind]; sex=test.dat[na.ind,"sex.in.ped"]
  phen1=test.dat[na.ind,phen]; id=test.dat[na.ind,"id"]
  if (!missing(covar)) x.covar <- as.matrix(test.dat[na.ind,covar]) 

  ###########################################
  count <- table(snp1)
  gntps <- names(count)
  count1 <- rep(0,3) 
  count1[as.numeric(gntps)+1] <- count
  if (count1[1]>count1[3]) MAC <- count1[2]+count1[3]*2 else MAC <- count1[2]+2*count1[1]
  ntotal <- sum(count)
  nmiss <- length(snp)-ntotal
  if (chr=="X") maf_ntotal <- maf.fun.x(snp1,sex) else maf_ntotal <- maf.fun(snp1)

  ###########################################  
  if (length(unique(phen1))<=1 ) {
      lme.out= list(ntotal,nmiss,maf_ntotal, NA, NA, NA,"single category in phen",NA,MAC,count1[1],count1[2],count1[3])
     return(lme.out)
  }
  
  ###########################################
  if (length(count)==1 ) {
     lme.out <- list(ntotal,nmiss,maf_ntotal, NA, NA, NA,"single category in SNP",NA,MAC,count1[1],count1[2],count1[3])
     return(lme.out)
  } 
 
  ###########################################
  if (length(unique(snp1))>1 && !missing(covar)){ 
    if (length(covar)>1) colinear <- apply(x.covar,2,cor.snp,x=snp1) else colinear <- cor.snp(x.covar,snp1)
    if (sum(colinear)>0){ 
       lme.out= list(ntotal,nmiss,maf_ntotal, NA, NA, NA,"collinarity",NA,MAC,count1[1],count1[2],count1[3])
       return(lme.out)
    }  
  }

#############################################
# Perform LME
#############################################
  if (!missing(covar)) lme.test <- try(lmekin(phen1 ~ snp1+x.covar+(1|id),varlist=kmat)) else lme.test <- try(lmekin(phen1 ~ snp1+(1|id),varlist=kmat))
  if (!"try-error" %in% class(lme.test)){ 
     Z <- lme.test$coef$fixed[2]/sqrt(lme.test$var[2,2])
     lme.out <- list(ntotal,nmiss,maf_ntotal,lme.test$coef$fixed[2],sqrt(lme.test$var[2,2]),Z,NA,pchisq(Z^2,1,lower.tail=F),MAC,count1[1],count1[2],count1[3])
  } else lme.out= list(ntotal,nmiss,maf_ntotal, NA, NA, NA,"try-error",NA,MAC,count1[1],count1[2],count1[3])	
  return(lme.out)
}

