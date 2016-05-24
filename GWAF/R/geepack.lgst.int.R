geepack.lgst.int=function(snp,phen,test.dat,covar,cov.int,sub="N"){ ####061409 sub
  ##get rid of  observations missing either genotype or phenotyp
   na.ind <- if (length(covar)>1) !is.na(snp) & !is.na(test.dat[,phen]) &  apply(!is.na(test.dat[,covar]),1,all) else
                        !is.na(snp) & !is.na(test.dat[,phen]) & !is.na(test.dat[,covar])
  snp1=snp[na.ind]
  maf=mean(snp1)/2 
  phen1=test.dat[na.ind,phen]
  famid=test.dat[na.ind,"famid"]
  x.covar=test.dat[na.ind,covar] 
  n=sum(na.ind) 
  model <- "additive"
  ###########################################
  #produce genotyp count
  count=table(snp1)
  gntps<-names(count)
  count1<-rep(0,3) # if count of certain gntp=0, match counts with gntps in ouput
  count1[as.numeric(gntps)+1]<-count

  #########################################  
  ###identify the # of categories in cov.int, stop if it is 1
  bin.flag <- length(table(test.dat[na.ind,cov.int])) 
  if (bin.flag==2) bin <- sort(as.numeric(names(table(test.dat[na.ind,cov.int])))) 
  if (bin.flag==1 ) { ####061409 sub
     if (sub=="Y") gee.out= matrix(c(cov.int,n,maf,rep(NA,15),"single categoty in cov.int"),ncol=1) else
        gee.out= matrix(c(cov.int,n,maf,rep(NA,10),"single categoty in cov.int"),ncol=1) 
     return(gee.out) 
  } 
       
  #check categories in y, if more than 2, stop the run
  length.unique<-length(unique(na.omit(phen1)))

  if (length.unique >2) 
	stop (paste("More than two categories in the phenotype: ",phen, 
                    ". Program stopped because of this error",sep="")) 

  if (length.unique<=1 ) {
   	if (bin.flag>2) {      
        	gee.out= matrix(c(cov.int,n,maf,rep(NA,10),"single category in phen"),ncol=1)  
   	} else { ####061409 sub
       if (sub=="Y") gee.out= matrix(c(cov.int,n,maf,rep(NA,15),"single categoty in phen"),ncol=1) else
          gee.out= matrix(c(cov.int,n,maf,rep(NA,10),"single categoty in phen"),ncol=1)
       }
  	return(gee.out)
  }
  
  #produce genotype count in affected (higher level)
  cnt.tbl<-table(phen1,snp1)
  gntps<-dimnames(cnt.tbl)$snp1
  count.d<-rep(0,3)
  count.d[as.numeric(gntps)+1]<-cnt.tbl[2,]
  nd=sum(count.d) 
  mafd <- (count.d %*% 0:2)/sum(count.d)/2

  ###############################################
 ##non-informative SNP; or only one category or less in y,skip
 if (length(count)==1 ) {
    #print("Fail 1.5")
    if (bin.flag>2) {      
      	gee.out= matrix(c(cov.int,n,maf,nd,mafd,rep(NA,8),"non-informative SNP"),ncol=1) 
   	} else { ####061409 sub
       if (sub=="Y") gee.out= matrix(c(cov.int,n,maf,nd,mafd,rep(NA,13),"non-informative SNP"),ncol=1) else
          gee.out= matrix(c(cov.int,n,maf,nd,mafd,rep(NA,8),"non-informative SNP"),ncol=1)
       }
    return(gee.out)
 }

 ##two genotypes, minimum count <10,skip 
 if (length(count)==2 && min(count)<10) {
    #print("Fail 2")
   ##do not run gee if genotype count <10 for 2 gntp SNP, 
    if (bin.flag>2) {      
      	gee.out= matrix(c(cov.int,n,maf,nd,mafd,rep(NA,8),"low genotype count"),ncol=1) 
   	} else { ####061409 sub
       if (sub=="Y") gee.out= matrix(c(cov.int,n,maf,nd,mafd,rep(NA,13),"low genotype count"),ncol=1) else 
          gee.out= matrix(c(cov.int,n,maf,nd,mafd,rep(NA,8),"low genotype count"),ncol=1)
       }
    return(gee.out)
 } 
   
 ##three genotypes, but the sum of the two lower counts < 10, skip
 if (length(count)==3 && (sort(count)[1]+sort(count)[2])<10) { 
   #print("Fail 3")
    if (bin.flag>2) {      
      	gee.out= matrix(c(cov.int,n,maf,nd,mafd,rep(NA,8),"low genotype count"),ncol=1)
   	} else { ####061409 sub
       if (sub=="Y")	gee.out= matrix(c(cov.int,n,maf,nd,mafd,rep(NA,13),"low genotype count"),ncol=1) else
          gee.out= matrix(c(cov.int,n,maf,nd,mafd,rep(NA,8),"low genotype count"),ncol=1)
       }
  return(gee.out)
 }  

 x.snp = snp1    
 x.int <- x.snp*test.dat[na.ind,cov.int]  #####071509
 if (length(unique(x.int))==1 ) { ##### 061209 
     if (bin.flag>2) {      
        gee.out= matrix(c(cov.int,n,maf,nd,mafd,rep(NA,8),"single category in interaction"),ncol=1)  
     } else { ####061409 sub
     if (sub=="Y") gee.out= matrix(c(cov.int,n,maf,nd,mafd,rep(NA,13),"single categoty in interaction"),ncol=1) else
        gee.out= matrix(c(cov.int,n,maf,nd,mafd,rep(NA,8),"single category in interaction"),ncol=1)  
     }
     return(gee.out)
 }  ######

 ##stop if sample MAF<0.01                      
 if (maf<0.01) {
    if (bin.flag>2) {      
      	gee.out= matrix(c(cov.int,n,maf,nd,mafd,rep(NA,8),"sample MAF<0.01"),ncol=1) 
   	} else { ####061409 sub
       if (sub=="Y")	gee.out= matrix(c(cov.int,n,maf,nd,mafd,rep(NA,13),"sample MAF<0.01"),ncol=1) else
          gee.out= matrix(c(cov.int,n,maf,nd,mafd,rep(NA,8),"sample MAF<0.01"),ncol=1)
       }
    return(gee.out)
 }

##check the collinearity between snp and covariates, stop if collinearity exists
#############################################
######NEW cor.snp handles factor variables by making them linear (not perfect solution):
  cor.snp <- function(y,x){
   if(!is.numeric(y))y<-as.numeric(as.factor(y))
   return(sd(y)==0 || abs(cor(y,x,use="complete"))>0.99999999 )}
###############################################
  if(length(table(x.snp))>1){ 
    if (length(covar)>1) {
       colinear <- apply(x.covar,2,cor.snp,x=x.snp) 
       colinear <- c(colinear,cor.snp(x.int,x.snp),cor.snp(x.covar[,cov.int],x.int),(length(table(x.snp))==2 & cor(x.snp,x.int)>0.95))  #####071509
    } else {
       colinear<-cor.snp(x.covar,x.snp)
       colinear <- c(colinear,cor.snp(x.int,x.snp),cor.snp(x.covar,x.int),(length(table(x.snp))==2 & cor(x.snp,x.int)>0.95))
    }
    if (sum(colinear)>0){
      	print("Fail 1")
       if (bin.flag>2) gee.out= matrix(c(cov.int,n,maf,nd,mafd,rep(NA,8),"collinearity"),ncol=1) else { ####061409 sub
           if (sub=="Y") gee.out=matrix(c(cov.int,n,maf,nd,mafd,rep(NA,13),"collinearity"),ncol=1) else
              gee.out=matrix(c(cov.int,n,maf,nd,mafd,rep(NA,8),"collinearity"),ncol=1)
           }
       return(gee.out)
    }
  }

#print("OK Running")
  
  # Perform GLM if number of clusters with 2 or more individuals <minfs2
  #or if minimum of the expect counts <5
  #or any zero gntp counts in affected or unaffected
  #or any 0 gntp counts in table(x.int,x.snp) when x.snp and cov.int both have 2 levels 
  #or sum of # of genotypes 1 and 2 is 0 in either category of bin.cov in table(bin.cov,x.snp), where bin.cov is a binary covariate
  #or x.snp and x.int have the same count in either unaffected or affected sample
  famsiz=table(table(famid))
  famsiz2=sum(famsiz[as.numeric(names(famsiz))>=2])#count # >2 sibs
  minfs2=10
      
  not.enough <- prod(table(phen1,x.snp))==0 #zero count of a gntp in affected/unaffected
  tab.phn <- identical(table(phen1,x.snp)[1,],table(phen1,x.int)[1,]) | identical(table(phen1,x.snp)[2,],table(phen1,x.int)[2,])

############################################################ deal with binary covariates 
  cell0 <- not.enough.subs <- F
  if (length(covar)==1) cat.covar <- length(unique(x.covar)) else cat.covar <- apply(x.covar,2,function(x)length(unique(x))) #####102208
  if (any(cat.covar==1)) {
     if (bin.flag>2) gee.out= matrix(c(cov.int,n,maf,nd,mafd,rep(NA,8),"single category in covariate(s)"),ncol=1) else { ####061409 sub
        if (sub=="Y") gee.out= matrix(c(cov.int,n,maf,nd,mafd,rep(NA,13),"single category in covariate(s)"),ncol=1) else
           gee.out= matrix(c(cov.int,n,maf,nd,mafd,rep(NA,8),"single category in covariate(s)"),ncol=1)
        }
     return(gee.out)
  }  else {
  if (any(cat.covar==2)) {
     bin.cov <- covar[cat.covar==2]
     if (length(bin.cov)==1) {
        if (length(covar)==1) cell0 <- sum(table(x.covar,x.snp)[1,-1])==0 | sum(table(x.covar,x.snp)[2,-1])==0 else #####102208
           cell0 <- sum(table(x.covar[,bin.cov],x.snp)[1,-1])==0 | sum(table(x.covar[,bin.cov],x.snp)[2,-1])==0 
     } else {           
           cell0 <- any(apply(x.covar[,bin.cov],2,function(x,snp=x.snp)sum(table(x,snp)[1,-1])==0 | sum(table(x,snp)[2,-1])==0,snp=x.snp))           
           if (bin.flag==2){
           bin.no.int <- bin.cov[bin.cov!=cov.int]
           for (i in bin.no.int){
               tab.bin.cov00 <- sum(table(phen1[x.covar[,i]==unique(x.covar[,i])[1] & x.covar[,cov.int]==bin[1]],
                  x.snp[x.covar[,i]==unique(x.covar[,i])[1] & x.covar[,cov.int]==bin[1]])[,-1])==0
               tab.bin.cov01 <- sum(table(phen1[x.covar[,i]==unique(x.covar[,i])[1] & x.covar[,cov.int]==bin[2]],
                  x.snp[x.covar[,i]==unique(x.covar[,i])[1] & x.covar[,cov.int]==bin[2]])[,-1])==0
               tab.bin.cov10 <- sum(table(phen1[x.covar[,i]==unique(x.covar[,i])[2] & x.covar[,cov.int]==bin[1]],
                  x.snp[x.covar[,i]==unique(x.covar[,i])[2] & x.covar[,cov.int]==bin[1]])[,-1])==0
               tab.bin.cov11 <- sum(table(phen1[x.covar[,i]==unique(x.covar[,i])[2] & x.covar[,cov.int]==bin[2]],
                  x.snp[x.covar[,i]==unique(x.covar[,i])[2] & x.covar[,cov.int]==bin[2]])[,-1])==0
               cell0 <- tab.bin.cov00 | tab.bin.cov01 | tab.bin.cov10 | tab.bin.cov11 | cell0
                                }
                           }
            }
  }
     }

############################################################# 
#############################################102408
 if (bin.flag==2) {
    x.cov.int <- test.dat[na.ind,cov.int]
    count.sub1 <- count.sub2 <- rep(0,3)
    count.sub1[as.numeric(names(table(x.snp[x.cov.int==bin[1]])))+1] <- table(x.snp[x.cov.int==bin[1]])
    count.sub2[as.numeric(names(table(x.snp[x.cov.int==bin[2]])))+1] <- table(x.snp[x.cov.int==bin[2]])
    if (sum(sort(count.sub1)[1:2])<10 | sum(sort(count.sub2)[1:2])<10) { ####061409 sub
       if (sub=="Y") gee.out= matrix(c(cov.int,n,maf,nd,mafd,rep(NA,13),"low genotype count in subset(s)"),ncol=1) else
          gee.out= matrix(c(cov.int,n,maf,nd,mafd,rep(NA,8),"low genotype count in subset(s)"),ncol=1)       
       return(gee.out)
    }
    not.enough.subs <- prod(table(phen1[x.cov.int==bin[1]],x.snp[x.cov.int==bin[1]]))*prod(table(phen1[x.cov.int==bin[2]],x.snp[x.cov.int==bin[2]]))==0
 }
#############################################

  if (famsiz2<=minfs2 | cell0 | (prod(table(x.int,x.snp))+bin.flag==2 & length(table(x.int,x.snp))==4) | tab.phn | not.enough.subs ||not.enough){ 
     if (min(chisq.test(table(phen1,x.snp))$expected)<5) warning="logistic reg & exp count<5" else warning="logistic reg"               
     xcovar<-as.matrix(test.dat[na.ind,covar])
     if (bin.flag>2) {
        gee.test<-try(glm(phen1~x.snp+x.int+xcovar, family="binomial")) 
        if (!"try-error" %in% class(gee.test)) {
           gee.test.s<-summary(gee.test)
           gee.out <- c(cov.int,n,maf,nd,mafd,gee.test.s$cov.scaled["x.snp","x.int"],model,gee.test.s$coef["x.snp",1:2],pchisq(gee.test.s$coef["x.snp",3]^2,1,lower.tail=F),
			  gee.test.s$coef["x.int",1:2],pchisq(gee.test.s$coef["x.int",3]^2,1,lower.tail=F),warning) #### add gee.test.s$cov.scaled["x.snp","x.int"] 091409
        } else gee.out <- c(cov.int,n,maf,nd,mafd,rep(NA,8),"try-error")
     } else { ####061409 sub
     if (sub=="Y"){
        tab1 <- table(test.dat[na.ind & test.dat[,cov.int]==bin[1],cov.int],snp[na.ind & test.dat[,cov.int]==bin[1]])
        tab2 <- table(test.dat[na.ind & test.dat[,cov.int]==bin[2],cov.int],snp[na.ind & test.dat[,cov.int]==bin[2]])
        if (length(tab1)==1 | length(tab2)==1) gee.out <- c(cov.int,n,maf,nd,mafd,rep(NA,13),"non-informative SNP in subset") else {
           if (sum(!covar %in% cov.int)>0) { ###########replace !covar %in% cov.int
              xcovar.bin<-as.matrix(test.dat[na.ind,covar[!covar %in% cov.int]])         
              gee.test1<-try(glm(phen1~x.snp+xcovar.bin, family="binomial", subset=test.dat[na.ind,cov.int]==bin[1])) 
              gee.test2<-try(glm(phen1~x.snp+xcovar.bin, family="binomial", subset=test.dat[na.ind,cov.int]==bin[2]))
              gee.test<-try(glm(phen1~x.snp+x.int+xcovar, family="binomial"))
           } else {
              gee.test1<-try(glm(phen1~x.snp, family="binomial", subset=test.dat[na.ind,cov.int]==bin[1])) 
              gee.test2<-try(glm(phen1~x.snp, family="binomial", subset=test.dat[na.ind,cov.int]==bin[2]))
              gee.test<-try(glm(phen1~x.snp+x.int+xcovar, family="binomial"))
           }                 
           if (!"try-error" %in% c(class(gee.test1),class(gee.test2),class(gee.test))) {
              gee.test.s1<-summary(gee.test1);gee.test.s2<-summary(gee.test2);gee.test.s<-summary(gee.test);
              gee.out <- c(cov.int,n,maf,nd,mafd,model,gee.test.s$coef["x.snp",1:2],pchisq(gee.test.s$coef["x.snp",3]^2,1,lower.tail=F),
                          gee.test.s1$coef["x.snp",1:2],pchisq(gee.test.s1$coef["x.snp",3]^2,1,lower.tail=F),
			     gee.test.s2$coef["x.snp",1:2],pchisq(gee.test.s2$coef["x.snp",3]^2,1,lower.tail=F),
                          gee.test.s$coef["x.int",1:2],pchisq(gee.test.s$coef["x.int",3]^2,1,lower.tail=F),warning)
           } else {
           if ("try-error" %in% c(class(gee.test1),class(gee.test2)) & !"try-error" %in% class(gee.test)) { ####061409 output interaction even subset analyses fail
              gee.test.s<-summary(gee.test);
              gee.out <- c(cov.int,n,maf,nd,mafd,model,gee.test.s$coef["x.snp",1:2],pchisq(gee.test.s$coef["x.snp",3]^2,1,lower.tail=F),                          
                          rep(NA,6),gee.test.s$coef["x.int",1:2],pchisq(gee.test.s$coef["x.int",3]^2,1,lower.tail=F),warning)
           } else gee.out <- c(cov.int,n,maf,nd,mafd,rep(NA,13),"try-error")
           }
        }
     } else {####061409 sub="N" or NA
        gee.test<-try(glm(phen1~x.snp+x.int+xcovar, family="binomial")) 
        if (!"try-error" %in% class(gee.test)) {
           gee.test.s<-summary(gee.test)
           gee.out <- c(cov.int,n,maf,nd,mafd,gee.test.s$cov.scaled["x.snp","x.int"],model,gee.test.s$coef["x.snp",1:2],pchisq(gee.test.s$coef["x.snp",3]^2,1,lower.tail=F),
			  gee.test.s$coef["x.int",1:2],pchisq(gee.test.s$coef["x.int",3]^2,1,lower.tail=F),warning)
        } else gee.out <- c(cov.int,n,maf,nd,mafd,rep(NA,8),"try-error")
     }
     }
     return(gee.out)
  } else {
##############################################################################################
#
# Perform GEE
#
##############################################################################################
  xcovar<-as.matrix(test.dat[na.ind,covar])
  if (bin.flag>2) {
     gee.test.v <- try(geese(phen1~x.snp+x.int+xcovar,id=famid,family="binomial",corstr="independence"))
     if (!"try-error" %in% class(gee.test.v)) {
        if (gee.test.v$error!=0) {
           if (min(chisq.test(table(phen1,x.snp))$expected)<5) warning="not converged and exp count<5" else warning="not converged"
        } else {
        if (min(chisq.test(table(phen1,x.snp))$expected)<5) warning="exp count<5" else warning=NA
        }
        gee.test <- summary(gee.test.v)
        gee.out <- c(cov.int,n,maf,nd,mafd,gee.test.v$vbeta[names(gee.test.v$beta)=="x.snp",names(gee.test.v$beta)=="x.int"],model,
                        unlist(gee.test$mean["x.snp",1:2]),pchisq(gee.test$mean["x.snp",3],1,lower.tail=F),unlist(gee.test$mean["x.int",1:2]),
                        pchisq(gee.test$mean["x.int",3],1,lower.tail=F),warning) #### add gee.test.s$cov.scaled["x.snp","x.int"] 091409
     } else gee.out <- c(cov.int,n,maf,nd,mafd,rep(NA,8),"try-error")
     if (gee.test.v$error!=0) gee.out= c(cov.int,n,maf,nd,mafd,rep(NA,8),"not converged")
  } else {
  if (sub=="Y") { ####061409 sub
     tab1 <- table(test.dat[na.ind & test.dat[,cov.int]==bin[1],cov.int],snp[na.ind & test.dat[,cov.int]==bin[1]])
     tab2 <- table(test.dat[na.ind & test.dat[,cov.int]==bin[2],cov.int],snp[na.ind & test.dat[,cov.int]==bin[2]])
     if (length(tab1)==1 | length(tab2)==1) gee.out <- c(cov.int,n,maf,nd,mafd,rep(NA,13),"non-informative SNP in subset") else { 
        if (sum(!covar %in% cov.int)>0) { ####replace !covar %in% cov.int
           xcovar.bin<-as.matrix(test.dat[na.ind,covar[!covar %in% cov.int]])
           gee.test1 <- try(summary(geese(phen1~x.snp+xcovar.bin,id=famid,family="binomial",corstr="independence",subset=test.dat[na.ind,cov.int]==bin[1])))
           gee.test2 <- try(summary(geese(phen1~x.snp+xcovar.bin,id=famid,family="binomial",corstr="independence",subset=test.dat[na.ind,cov.int]==bin[2])))
           gee.test <- try(summary(geese(phen1~x.snp+x.int+xcovar,id=famid,family="binomial",corstr="independence")))
        } else {
          gee.test1 <- try(summary(geese(phen1~x.snp,id=famid,family="binomial",corstr="independence",subset=test.dat[na.ind,cov.int]==bin[1])))
          gee.test2 <- try(summary(geese(phen1~x.snp,id=famid,family="binomial",corstr="independence",subset=test.dat[na.ind,cov.int]==bin[2])))
          gee.test <- try(summary(geese(phen1~x.snp+x.int+xcovar,id=famid,family="binomial",corstr="independence")))
        }

        if (!"try-error" %in% c(class(gee.test1),class(gee.test2),class(gee.test))) { ###061209
           if (gee.test1$error==0 & gee.test2$error==0 & gee.test$error==0) {
               warning=NA                         
               gee.out <- c(cov.int,n,maf,nd,mafd,model,unlist(gee.test$mean["x.snp",1:2]),pchisq(gee.test$mean["x.snp",3],1,lower.tail=F),
                           unlist(gee.test1$mean["x.snp",1:2]),pchisq(gee.test1$mean["x.snp",3],1,lower.tail=F),unlist(gee.test2$mean["x.snp",1:2]),
                           pchisq(gee.test2$mean["x.snp",3],1,lower.tail=F),unlist(gee.test$mean["x.int",1:2]),pchisq(gee.test$mean["x.int",3],1,lower.tail=F),warning)
           } else 
           if (gee.test1$error+gee.test2$error!=0 & gee.test$error==0) {
              warning="not converged only in at least one subset analysis"
              gee.out <- c(cov.int,n,maf,nd,mafd,model,unlist(gee.test$mean["x.snp",1:2]),pchisq(gee.test$mean["x.snp",3],1,lower.tail=F),rep(NA,6),
                          unlist(gee.test$mean["x.int",1:2]),pchisq(gee.test$mean["x.int",3],1,lower.tail=F),warning)
           } else {
              warning="not converged" 
              gee.out <- c(cov.int,n,maf,nd,mafd,rep(NA,13),warning)               
           }
        } else {
           if ("try-error" %in% c(class(gee.test1),class(gee.test2)) & !"try-error" %in% class(gee.test)) { ####061209 output interaction even subset analyses fail
              if (gee.test$error!=0) {
                 warning="not converged" 
                 gee.out <- c(cov.int,n,maf,nd,mafd,rep(NA,13),warning)
              } else {
                 warning=NA                         
                 gee.out <- c(cov.int,n,maf,nd,mafd,model,unlist(gee.test$mean["x.snp",1:2]),pchisq(gee.test$mean["x.snp",3],1,lower.tail=F),rep(NA,6),
                             unlist(gee.test$mean["x.int",1:2]),pchisq(gee.test$mean["x.int",3],1,lower.tail=F),warning)
              }
           } else gee.out <- c(cov.int,n,maf,nd,mafd,rep(NA,13),"try-error")
        }
     }
  } else { ####061409 sub="N" or NA
    gee.test.v <- try(geese(phen1~x.snp+x.int+xcovar,id=famid,family="binomial",corstr="independence"))
    if (!"try-error" %in% class(gee.test.v)) {
       if (gee.test.v$error!=0) warning="not converged" else warning=NA              
       gee.test <- summary(gee.test.v)
       gee.out <- c(cov.int,n,maf,nd,mafd,gee.test.v$vbeta[names(gee.test.v$beta)=="x.snp",names(gee.test.v$beta)=="x.int"],model,
                   unlist(gee.test$mean["x.snp",1:2]),pchisq(gee.test$mean["x.snp",3],1,lower.tail=F),unlist(gee.test$mean["x.int",1:2]),
                   pchisq(gee.test$mean["x.int",3],1,lower.tail=F),warning)
    } else gee.out <- c(cov.int,n,maf,nd,mafd,rep(NA,8),"try-error")
    if (gee.test.v$error!=0) gee.out= c(cov.int,n,maf,nd,mafd,rep(NA,8),"not converged")
  }
  }
      return(gee.out)
  }
}     