glmm.EC <- function(snp,phen,test.dat,covar,chr){

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
  phen1=test.dat[na.ind,phen];  famid=test.dat[na.ind,"famid"]
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

  #########################################  
  ###########################################  
  if (length(unique(phen1))<=1 | length(unique(phen1))>2) {
     glmm.out <- list(ntotal,nmiss,maf_ntotal, NA, NA, NA,"single category in phen",NA,MAC,count1[1],count1[2],count1[3])
     return(glmm.out)
  }
  
  ###########################################
  if (length(count)==1 ) {
     glmm.out <- list(ntotal,nmiss,maf_ntotal, NA, NA, NA,"single category in SNP",NA,MAC,count1[1],count1[2],count1[3])
     return(glmm.out)
  } 
 
  ###########################################
  if (length(unique(snp1))>1 && !missing(covar)){ 
     if (length(covar)>1) colinear <- apply(x.covar,2,cor.snp,x=snp1) else colinear <- cor.snp(x.covar,snp1)
     if (sum(colinear)>0){ 
        glmm.out= list(ntotal,nmiss,maf_ntotal, NA, NA, NA,"collinarity",NA,MAC,count1[1],count1[2],count1[3])
        return(glmm.out)
     }  
  }

#############################################################################################
#
# Perform GLMM
#
#############################################################################################

  if (!missing(covar)){
     gee.test.0 <- try(glmer(phen1 ~ x.covar+(1|famid), family="binomial",control=glmerControl(optimizer="bobyqa")))
     gee.test <- try(glmer(phen1 ~ snp1+x.covar+(1|famid), family="binomial",control=glmerControl(optimizer="bobyqa")))
  } else {
     gee.test.0 <- try(glmer(phen1 ~ (1|famid), family="binomial",control=glmerControl(optimizer="bobyqa")))
     gee.test <- try(glmer(phen1 ~ snp1+(1|famid), family="binomial",control=glmerControl(optimizer="bobyqa")))
  } 

  if (!"try-error" %in% class(gee.test.0) && !"try-error" %in% class(gee.test)){ 
     conv <- is.null(unlist(gee.test@optinfo$conv$lme4)) | (!is.null(unlist(gee.test@optinfo$conv$lme4)) && gee.test@optinfo$conv$lme4$code==0)
     if (is.null(unlist(gee.test@optinfo$warnings)) & conv & gee.test@optinfo$conv$opt==0) warning <- NA else warning <- "converging issue"
     lrt <- -2*(logLik(gee.test.0)-logLik(gee.test))[1]
     gee.test.s <- try(summary(gee.test))
     if (!"try-error" %in% class(gee.test.s)){ 
        if (lrt>0) Z <- sign(coef(gee.test.s)["snp1",1])*sqrt(lrt) else Z <- 0
        p <- pchisq(lrt,1,lower.tail=F)
        glmm.out <- list(ntotal,nmiss,maf_ntotal,coef(gee.test.s)["snp1",1],coef(gee.test.s)["snp1",2],Z,warning,p,MAC,count1[1],count1[2],count1[3])
     } else glmm.out= list(ntotal,nmiss,maf_ntotal, NA, NA, NA,"try-error",NA,MAC,count1[1],count1[2],count1[3])
  } else glmm.out= list(ntotal,nmiss,maf_ntotal, NA, NA, NA,"try-error",NA,MAC,count1[1],count1[2],count1[3])
  return(glmm.out)
}

