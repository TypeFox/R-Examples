glmm.lgst=function(snp,phen,test.dat,covar=NULL,model="a"){
  #print(str(test.dat))
  ##get ride of  observations missing either genotype or phenotyp
   na.ind <-  if(missing(covar)) !is.na(snp) & !is.na(test.dat[,phen]) else if
                        (length(covar)>1) 
			!is.na(snp) & !is.na(test.dat[,phen]) &  
			apply(!is.na(test.dat[,covar]),1,all) else
                        !is.na(snp) & !is.na(test.dat[,phen]) &
                        !is.na(test.dat[,covar])
  snp1=snp[na.ind]
  phen1=test.dat[na.ind,phen]
  famid=test.dat[na.ind,"famid"]
  if(!missing(covar)) x.covar=test.dat[na.ind,covar] ###082608

  ###########################################
  #produce genotyp count
  count=table(snp1)
  gntps<-names(count)
  count1<-rep(0,3) # if count of certain gntp=0, match counts with gntps in ouput
  count1[as.numeric(gntps)+1]<-count
  #print(count1)

  #########################################  
  #check categories in y, if more than 2, stop the run
  length.unique<-length(unique(na.omit(phen1)))

  if (length.unique >2) 
	stop (paste("More than two categories in the phenotype: ",phen, 
                    ". Program stopped because of this error",sep="")) 

  if (length.unique<=1 ) {
   	if (model %in% c("a","d","r","fa")) {      #16 columns #######fa
        	gee.out= matrix(c(count1,rep(NA,11),"single category in phen",NA),ncol=1)
   	}else   #20 columns
        	gee.out= matrix(c(count1,rep(NA,15),"single categoty in phen",NA),ncol=1)

  	return(gee.out)
  }
  
  #produce genotype count in affected (higher level)
  cnt.tbl<-table(phen1,snp1)
  gntps<-dimnames(cnt.tbl)$snp1
  count.d<-rep(0,3)
  count.d[as.numeric(gntps)+1]<-cnt.tbl[2,]

  ###############################################
 ##non-informative SNP; or only one category or less in y,skip#####
 if (length(count)==1 ) {
    print("Fail 1.5")
   if (model %in% c("a","d","r","fa")) {      #######fa
        gee.out= matrix(c(count1,count.d,rep(NA,10)),ncol=1)
   }else
        gee.out= matrix(c(count1,count.d,rep(NA,14)),ncol=1)
  return(gee.out)
 } else if (length(count)==2) {
     if (model %in% c("a","fa")) {x.snp = snp1     #######fa
     } else if (model=="g") {x.snp = ifelse(snp1!=2,snp1,1); 
     }else if (model=="d") {x.snp = ifelse(snp1!=2,snp1,1); 
                          if (length(table(x.snp))==1) {
                             gee.out= matrix(c(count1,count.d,rep(NA,10)),ncol=1)
                             return(gee.out)
                          }
     }else if (model=="r" & !count1[3]<10) {x.snp=ifelse(snp1==2,1,0)
     }else if (model=="r" & count1[3]<10) {gee.out= matrix(c(count1,count.d,rep(NA,10)),ncol=1)
                                           return(gee.out)
   }     
 } else if (length(count)==3) {
     if (model %in% c("a","fa")) {x.snp = snp1    #######fa
     } else if (model=="g") {x.snp = factor(snp1)
     } else if (model=="d") {x.snp = ifelse(snp1!=2,snp1,1)
     } else if (model=="r") {x.snp=ifelse(snp1==2,1,0)
   } 

 }

 ## if covariates, check the collinearity, skip of collinearity
  #function to check collinarity
#  cor.snp <- function(y,x) sd(y)==0 || abs(cor(y,x,use="complete"))>0.99999999 
#############################################
######NEW cor.snp handles factor variables by making them linear (not perfect solution):
  cor.snp <- function(y,x){
   if(!is.numeric(y))y<-as.numeric(as.factor(y))
   return(sd(y)==0 || abs(cor(y,x,use="complete"))>0.99999999 )}
###############################################
  if(length(table(x.snp))>1 && !missing(covar) && model!="g" | (length(table(x.snp))==2 && !missing(covar) && model=="g")){ #090908
  #  print(covar)
#############################################
#  if covars are non-numeric, as.matrix() will bomb.  Don't need to make it a matrix to do the testing:
#    x.covar=as.matrix(test.dat[na.ind,covar])
#############################################
#    x.covar=test.dat[na.ind,covar] #######082608
#############################################
# Fix to work if there is only one covariate:  
###############################################
#    colinear <- apply(x.covar,2,cor.snp,x=x.snp)
    if(length(covar)>1)colinear <- apply(x.covar,2,cor.snp,x=x.snp) else colinear<-cor.snp(x.covar,x.snp)
    ##check colinearity
    if( sum(colinear)>0){ ## only 1 gntp | colinearity between covar and SNP
      	print("Fail 1")
      	if (model %in% c("a","d","r","fa")) {         #######fa
        	gee.out= matrix(c(count1,count.d,rep(NA,8),"collinarity",NA),ncol=1) 
      	}else
        	gee.out= matrix(c(count1,count.d,rep(NA,12),"collinarity",NA),ncol=1)
    return(gee.out)
    }
  }

#print("OK Running")
  
### test differential missingness between cases and controls
  snp.na <- rep(1,length(snp))
  snp.na[is.na(snp)] <- 0
  if (length(table(snp.na))==2) miss.mat <- matrix(table(test.dat[,phen],snp.na),2,2) else miss.mat <- cbind(c(0,0),table(test.dat[,phen],snp.na))
  miss.diff.p <- fisher.test(miss.mat)$p.value
  if (miss.diff.p>1) miss.diff.p <- 1
  miss.0 <- miss.mat[1,1]/sum(miss.mat[1,])
  miss.1 <- miss.mat[2,1]/sum(miss.mat[2,])

##############################################################################################
#
# Perform GLMM
#
 ################################change 4 ends here
  if (model !="g" || (model=="g" & length(table(x.snp))==2)) {   #########################change 5 starts and ends here
	if (!missing(covar)){
#####################################
          x.covar <- as.matrix(test.dat[na.ind,covar])
 	   gee.test <- try(summary(lmer(phen1 ~ x.snp+x.covar+(1|famid), family="binomial")))
#####################################
   	} else {
 	   gee.test <- try(summary(lmer(phen1 ~ x.snp+(1|famid), family="binomial")))
  	} 
 	if (!"try-error" %in% class(gee.test)){ 
          if (!model %in% c("a")){       
             if (min(chisq.test(table(phen1,x.snp))$expected)<5) warning="exp count<5" else warning=NA
	   }  else {
             x.snp.d = ifelse(snp1!=2,snp1,1)
	      if (min(chisq.test(table(phen1,x.snp.d))$expected)<5) warning="exp count<5" else warning=NA
	      }
          chisq=as.numeric((coef(gee.test)["x.snp",1]/coef(gee.test)["x.snp",2])^2)
	   if (model !="g" ) {      
             if (model=="fa") model <- "a"  
	      gee.out <- matrix(c(count1,count.d,miss.0,miss.1,miss.diff.p,coef(gee.test)["x.snp",1:2],chisq,"1",model,warning,pchisq(chisq,1,lower.tail=F)),ncol=1)
   	   } else gee.out <- matrix(c(count1,count.d,miss.0,miss.1,miss.diff.p,coef(gee.test)["x.snp",1],rep(NA,2),coef(gee.test)["x.snp",2],rep(NA,2),chisq,"1","d",warning,pchisq(chisq,1,lower.tail=F)),ncol=1)  		  
	} else {
	    if (model %in% c("a","d","r","fa")) gee.out= matrix(c(count1,count.d,rep(NA,8),"try-error",NA),ncol=1) else gee.out= matrix(c(count1,count.d,rep(NA,12),"try-error",NA),ncol=1)		
	} 
  } else {#model="g" & 3 levels
   	if (!missing(covar)){
     	   x.covar=as.matrix(test.dat[na.ind,covar])
 	   gee.test <- try(summary(lmer(phen1 ~ x.snp+x.covar+(1|famid), family="binomial")))
   	} else {
 	   gee.test <- try(summary(lmer(phen1 ~ x.snp+(1|famid), family="binomial")))
   	}
       if (!"try-error" %in% class(gee.test)){ ######changed from if (class(gee.test) != "try-error")  
          if (min(chisq.test(table(phen1,x.snp))$expected)<5) warning="exp count<5" else warning=NA
 	   chisq=as.numeric(t(matrix(coef(gee.test)[c("x.snp1","x.snp2"),1],ncol=1))%*%solve(vcov(gee.test)[2:3,2:3])%*%matrix(coef(gee.test)[c("x.snp1","x.snp2"),1],ncol=1))  				
	   gee.out=matrix(c(count1,count.d,miss.0,miss.1,miss.diff.p,coef(gee.test)["x.snp1",1],coef(gee.test)["x.snp2",1],coef(gee.test)["x.snp2",1]-coef(gee.test)["x.snp1",1],
		   coef(gee.test)["x.snp1",2],coef(gee.test)["x.snp2",2],sqrt(vcov(gee.test)[2,2]+vcov(gee.test)[3,3]-2*vcov(gee.test)[2,3]),chisq,"2",model,warning,pchisq(chisq,2,lower.tail=F)),ncol=1)
	} else gee.out= matrix(c(count1,count.d,rep(NA,12),"try-error",NA),ncol=1)
  } 
  return(gee.out)
}

