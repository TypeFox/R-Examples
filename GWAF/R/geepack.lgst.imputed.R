geepack.lgst.imputed=function(snp,phen,test.dat,covar=NULL){
  ##get rid of  observations missing either genotype or phenotyp
   na.ind <-  if(is.null(covar)) !is.na(snp) & !is.na(test.dat[,phen]) else if
                        (length(covar)>1) 
			!is.na(snp) & !is.na(test.dat[,phen]) &  
			apply(!is.na(test.dat[,covar]),1,all) else
                        !is.na(snp) & !is.na(test.dat[,phen]) &
                        !is.na(test.dat[,covar])
  snp1=snp[na.ind]
  phen1=test.dat[na.ind,phen]
  famid=test.dat[na.ind,"famid"]
  if (!is.null(covar)) x.covar=test.dat[na.ind,covar] 
  n=length(snp1)
  nd=sum(phen1==max(phen1))
  ###########################################
  #produce genotyp count
  count=table(round(snp1))
  gntps<-names(count)
  count1<-rep(0,3) # if count of certain gntp=0, match counts with gntps in ouput
  count1[as.numeric(gntps)+1]<-count

  #produce genotype count in affected (higher level)
  cnt.tbl<-table(phen1,round(snp1))
  gntps<-unlist(dimnames(cnt.tbl)[2])
  count.d<-rep(0,3)
  count.d[as.numeric(gntps)+1]<-cnt.tbl[2,]
  #########################################  
  #check categories in y, if more than 2, stop the run, if <=1, output NA
  length.unique<-length(unique(na.omit(phen1)))

  if (length.unique >2) 
	stop (paste("More than two categories in the phenotype: ",phen, 
                    ". Program stopped because of this error",sep="")) 

  if (length.unique<=1 ) {
      	gee.out= matrix(c(n,rep(NA,5),"single category in phen",NA),ncol=1)
  	return(gee.out)
  }
  
  x.snp <- snp1
  imaf <- mean(x.snp)/2
  imafd <- mean(x.snp[phen1==max(phen1)])/2
  ###############################################
 ##non-informative SNP; 
  if (length(count)==1 ) {
     #print("Fail 1.5")
     gee.out= matrix(c(n,nd,imaf,imafd,rep(NA,4)),ncol=1)
     return(gee.out)
  }

 ##two genotypes, minimum count <10,skip 
  if (length(count)==2 && min(count)<10) {
     #print("Fail 2")
   ##do not run gee if genotype count <10 for 2 gntp SNP, 
     gee.out= matrix(c(n,nd,imaf,imafd,rep(NA,4)),ncol=1)
     return(gee.out)
  } 

 ##three genotypes, but the sum of the two lower counts < 10, skip
  if (length(count)==3 && (sort(count)[1]+sort(count)[2])<10) { 
     #print("Fail 3")
     gee.out= matrix(c(n,nd,imaf,imafd,rep(NA,4)),ncol=1)
     return(gee.out)
  } 

## if covariates, check the collinearity, skip of collinearity
#function to check collinarity
#############################################
######NEW cor.snp handles factor variables by making them linear (not perfect solution):
  cor.snp <- function(y,x){
   if(!is.numeric(y))y<-as.numeric(as.factor(y))
   return(sd(y)==0 || abs(cor(y,x,use="complete"))>0.99999999 )}
###############################################
  if (length(count)>1 && !is.null(covar)){ 
     if (length(covar)>1) colinear <- apply(x.covar,2,cor.snp,x=x.snp) else colinear<-cor.snp(x.covar,x.snp)
    ##check colinearity
     if (sum(colinear)>0){ ## only 1 gntp | colinearity between covar and SNP
      	 #print("Fail 1")
        gee.out= matrix(c(n,nd,imaf,imafd,rep(NA,2),"collinarity",NA),ncol=1) 
        return(gee.out)
     }
  }
  # Perform GLM if number of clusters with 2 or more individuals <minfs2
  #or if minimum of the expect counts <5
  #or any zero gntp counts in affected or unaffected
  famsiz=table(table(famid))
  famsiz2=sum(famsiz[as.numeric(names(famsiz))>=2])#count # >2 sibs
  minfs2=10
      
  tmp.table <- apply(cnt.tbl>0,1,sum)      
  not.enough <- sum(tmp.table<length(table(round(x.snp)))) > 0 #zero count of a gntp in affected/unaffected

############################################################ deal with binary covariates 
  cell0 <- F
  if (!is.null(covar)){
     if (length(covar)==1) cat.covar <- length(unique(x.covar)) else cat.covar <- apply(x.covar,2,function(x)length(unique(x))) 
     if (any(cat.covar==1)) {
        gee.out= matrix(c(n,nd,imaf,imafd,rep(NA,2),"non-informative covariate",NA),ncol=1) 
        return(gee.out)
     } else {
          if (any(cat.covar==2)) {
             bin.cov <- covar[cat.covar==2]
             if (length(bin.cov)==1) {
                if (length(covar)==1) cell0 <- sum(table(x.covar,round(x.snp))[1,-1])==0 | sum(table(x.covar,round(x.snp))[2,-1])==0 else 
                   cell0 <- sum(table(x.covar[,bin.cov],round(x.snp))[1,-1])==0 | sum(table(x.covar[,bin.cov],round(x.snp))[2,-1])==0 
             } else cell0 <- any(apply(x.covar[,bin.cov],2,function(x,snp=round(x.snp))sum(table(x,snp)[1,-1])==0 | sum(table(x,snp)[2,-1])==0,snp=round(x.snp)))
          }     
       }
  }
############################################################# 

  if (famsiz2<=minfs2 | cell0 || not.enough ){ 
     if (min(chisq.test(table(phen1,round(x.snp)))$expected)<5) warning="logistic reg & exp count<5" else warning="logistic reg"
     if (!is.null(covar)){
#####################################
# Change how this is specified so that it works with factor covars:
        x.covar<-as.data.frame(test.dat[na.ind,covar])
        dat1<-cbind(phen1,x.covar,x.snp)
        form<-as.formula(paste("phen1~",paste(colnames(dat1)[2:length(colnames(dat1))],sep="+",collapse="+")))
        gee.test<-try(glm(form,data=dat1, family="binomial", na.action=na.omit)) 
#####################################
     } else {
          gee.test <- try(glm(phen1 ~ x.snp, family="binomial", na.action=na.omit))
	}
     if (!"try-error" %in% class(gee.test)){  
	 gee.test.s<-summary(gee.test)
        gee.out <- matrix(c(n,nd,imaf,imafd,gee.test.s$coef["x.snp",1:2],warning,gee.test.s$coef["x.snp",4]),ncol=1)
	} else gee.out= matrix(c(n,nd,imaf,imafd,rep(NA,2),"try-error",NA),ncol=1) 
     return(gee.out)
  }

##############################################################################################
#
# Perform GEE
#
 ################################
  if (!is.null(covar)){
#####################################
# Change how this is specified so that it works with factor covars:
     x.covar<-as.data.frame(test.dat[na.ind,covar])
     dat1<-cbind(phen1,x.covar,x.snp)
     form<-as.formula(paste("phen1~",paste(colnames(dat1)[2:length(colnames(dat1))],sep="+",collapse="+")))               
		gee.test<-try(summary(geese(form,data=dat1,id=famid,family="binomial",na.action=na.omit,corstr="independence")))
#####################################
  } else {
		gee.test <- try(summary(geese(phen1 ~ x.snp, id=famid, family="binomial", na.action=na.omit,corstr="independence")))
    } ## end if/else(!is.null(covar))
  if (!"try-error" %in% class(gee.test)){                                              
     if (gee.test$error!=0){
	 if (min(chisq.test(table(phen1,round(x.snp)))$expected)<5) warning="not converged and exp count<5" else warning="not converged"
     } else {
          if (min(chisq.test(table(phen1,round(x.snp)))$expected)<5) warning="exp count<5" else warning=NA
	}
     x.snp.pos <- substr(rownames(gee.test$mean),start=1,stop=5)=="x.snp"
     chisq=try(gee.test$mean[x.snp.pos,"wald"])
 
     ##summarize output
     if (class(chisq)!="try-error"){
        gee.out <- matrix(c(n,nd,imaf,imafd,gee.test$mean["x.snp","estimate"],gee.test$mean["x.snp","san.se"],warning,pchisq(chisq,1,lower.tail=F)),ncol=1)  
     } else {##try-error chisq
        gee.out= matrix(c(n,nd,imaf,imafd,rep(NA,2),"try-error",NA),ncol=1)
	}
  } else {##try-error gee.test
      	gee.out= matrix(c(n,nd,imaf,imafd,rep(NA,2),"try-error",NA),ncol=1)
    }
  return(gee.out) 
}


