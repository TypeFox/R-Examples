geepack.lgst.int.imputed=function(snp,phen,test.dat,covar,cov.int,sub="N"){ ####061209 sub
  ##get rid of  observations missing either genotype or phenotyp
   na.ind <- if (length(covar)>1) !is.na(snp) & !is.na(test.dat[,phen]) &  apply(!is.na(test.dat[,covar]),1,all) else
                        !is.na(snp) & !is.na(test.dat[,phen]) & !is.na(test.dat[,covar])

  snp1=snp[na.ind]
  phen1=test.dat[na.ind,phen]
  famid=test.dat[na.ind,"famid"]
  x.covar=test.dat[na.ind,covar]  
  model <- "additive"
  bin.flag <- length(table(test.dat[na.ind,cov.int])) 
  if (bin.flag==2) bin <- sort(as.numeric(names(table(test.dat[na.ind,cov.int])))) 
  if (bin.flag==1 ) { 
     if (sub=="Y") gee.out= matrix(c(cov.int,rep(NA,17),"single categoty in cov.int"),ncol=1) else 
        gee.out= matrix(c(cov.int,rep(NA,12),"single categoty in cov.int"),ncol=1)      
     return(gee.out) 
  } 
  
  count1<-rep(0,3) # if count of certain gntp=0, match counts with gntps in ouput
    
  bigN = length(snp1) 
  imaf = mean(snp1)/2 
  varimaf = var(snp1)/2 
  
  bigNd=table(phen1)[2]
  imafd=(tapply(snp1,phen1,mean,na.rm=T)[2])/2 
  varimafd=(tapply(snp1,phen1,var,na.rm=T)[2])/2 
  imafc=(tapply(snp1,phen1,mean,na.rm=T)[1])/2  

  count1[1]=bigN
  count1[2]=imaf
  count1[3]=varimaf
  
  count.d=count1
     
  count.d[1]=bigNd
  count.d[2]=imafd
  count.d[3]=varimafd
      
  #########################################  
      
  testiMAF=2*bigN*imaf*(1-imaf)
  testiMAFd=2*bigNd*imafd*(1-imafd)
  testiMAFc=2*(bigN-bigNd)*imafc*(1-imafc)
      
  if (testiMAF+bigN*imaf^2 <=10) {
      #print("Fail 1.5 - expected (heterozygote+minor allele homozygote) <=10")
      if (bin.flag>2) {      
         gee.out= matrix(c(cov.int,bigN,imaf,bigNd,imafd,rep(NA,9)),ncol=1)  
      } else { ####061209 sub
      if (sub=="Y") gee.out= matrix(c(cov.int,bigN,imaf,bigNd,imafd,rep(NA,14)),ncol=1) else
         gee.out= matrix(c(cov.int,bigN,imaf,bigNd,imafd,rep(NA,9)),ncol=1)  
      }
      return(gee.out)
  }

  if (testiMAFd < 1) { #####20090125
      #print("Fail 1.6 - expected heterozygote<1 in cases")
      if (bin.flag>2) {      
       	gee.out= matrix(c(cov.int,bigN,imaf,bigNd,imafd,rep(NA,9)),ncol=1)  
      } else { ####061209 sub
      if (sub=="Y") gee.out= matrix(c(cov.int,bigN,imaf,bigNd,imafd,rep(NA,14)),ncol=1) else
         gee.out= matrix(c(cov.int,bigN,imaf,bigNd,imafd,rep(NA,9)),ncol=1)  
      }
      return(gee.out)
  }
  
  length.unique<-length(unique(na.omit(phen1)))
  
  ###########################################
  #produce genotyp count
  
  if (length.unique >2) 
	stop (paste("More than two categories in the phenotype: ",phen, 
                    ". Program stopped because of this error",sep="")) 

  if (length.unique<=1 ) {
      if (bin.flag>2) {      
       	gee.out= matrix(c(cov.int,bigN,imaf,bigNd,imafd,rep(NA,8),"single category in phen"),ncol=1)  
      } else { ####061209 sub
      if (sub=="Y") gee.out= matrix(c(cov.int,bigN,imaf,bigNd,imafd,rep(NA,13),"single categoty in phen"),ncol=1) else 
         gee.out= matrix(c(cov.int,bigN,imaf,bigNd,imafd,rep(NA,8),"single category in phen"),ncol=1)  
      }
      return(gee.out)  
  }
  
  ###############################################
  ##non-informative SNP; or only one category or less in y,skip
  x.snp = snp1
  x.int <- x.snp*test.dat[na.ind,cov.int]  ####071509

  if (length(unique(x.snp))==1 ) { ######05132009 imputed.X
      if (bin.flag>2) {      
       	gee.out= matrix(c(cov.int,bigN,imaf,bigNd,imafd,rep(NA,8),"single category in SNP"),ncol=1)  
      } else { ####061209 sub
      if (sub=="Y") gee.out= matrix(c(cov.int,bigN,imaf,bigNd,imafd,rep(NA,13),"single categoty in SNP"),ncol=1) else
         gee.out= matrix(c(cov.int,bigN,imaf,bigNd,imafd,rep(NA,8),"single category in SNP"),ncol=1)  
      }
    	return(gee.out)
  }  ######05132009 imputed.X

  if (length(unique(x.int))==1 ) { ##### 061209 
      if (bin.flag>2) {      
       	gee.out= matrix(c(cov.int,bigN,imaf,bigNd,imafd,rep(NA,8),"single category in interaction"),ncol=1)  
      } else {
      if (sub=="Y") gee.out= matrix(c(cov.int,bigN,imaf,bigNd,imafd,rep(NA,13),"single categoty in interaction"),ncol=1) else
         gee.out= matrix(c(cov.int,bigN,imaf,bigNd,imafd,rep(NA,8),"single category in interaction"),ncol=1)  
      }
    	return(gee.out)
  }  ######

 if (imaf<0.01) {
    if (bin.flag>2) {      
      	gee.out= matrix(c(cov.int,bigN,imaf,bigNd,imafd,rep(NA,8),"sample iMAF<0.01"),ncol=1)
    } else { ####061209 sub
    if (sub=="Y") gee.out= matrix(c(cov.int,bigN,imaf,bigNd,imafd,rep(NA,13),"sample iMAF<0.01"),ncol=1) else
       gee.out= matrix(c(cov.int,bigN,imaf,bigNd,imafd,rep(NA,8),"sample iMAF<0.01"),ncol=1)
    }
    return(gee.out)
 }
                            
##check the collinearity, skip of collinearity
#############################################
######NEW cor.snp handles factor variables by making them linear (not perfect solution):
  cor.snp <- function(y,x){
   if(!is.numeric(y))y<-as.numeric(as.factor(y))
   return(sd(y)==0 || abs(cor(y,x,use="complete"))>0.99999999 )}
###############################################
  if(length(table(x.snp))>1){
    if (length(covar)>1)colinear <- apply(x.covar,2,cor.snp,x=x.snp) else colinear<-cor.snp(x.covar,x.snp)
    colinear <- c(colinear,cor.snp(x.int,x.snp),cor.snp(x.covar[,cov.int],x.int),(length(table(x.snp))==2 & cor(x.snp,x.int)>0.95))  #####071509
    if (sum(colinear)>0){
      	#print("Fail 1")
       if (bin.flag>2) gee.out= matrix(c(cov.int,bigN,imaf,bigNd,imafd,rep(NA,8),"colinearity"),ncol=1) else { ####061209 sub
          if (sub=="Y") gee.out= matrix(c(cov.int,bigN,imaf,bigNd,imafd,rep(NA,13),"colinearity"),ncol=1) else
             gee.out= matrix(c(cov.int,bigN,imaf,bigNd,imafd,rep(NA,8),"colinearity"),ncol=1) 
       }
       return(gee.out)
    }
  }

#############################################102208
 if (bin.flag==2) {
    x.cov.int <- test.dat[na.ind,cov.int]
    imafd.sub1 <- tapply(x.snp[x.cov.int==bin[1]],phen1[x.cov.int==bin[1]],mean)[2]
    imafd.sub2 <- tapply(x.snp[x.cov.int==bin[2]],phen1[x.cov.int==bin[2]],mean)[2]
    sub.flag <- sum(x.cov.int==bin[1] & phen1==1)*imafd.sub1*(2-imafd.sub1)<5 | sum(x.cov.int==bin[2] & phen1==1)*imafd.sub2*(2-imafd.sub2)<5
    if (sub.flag) { ####061209 sub
       if (sub=="Y") gee.out= matrix(c(cov.int,bigN,imaf,bigNd,imafd,rep(NA,13),"low expected count in affected subset(s)"),ncol=1) else 
          gee.out= matrix(c(cov.int,bigN,imaf,bigNd,imafd,rep(NA,8),"low expected count in affected subset(s)"),ncol=1)    
    return(gee.out)
    }
 }
#############################################
#print("OK Running")
  
### test differential missingness between cases and controls

  # Perform GLM if number of clusters with 2 or more individuals <minfs2
  #or if minimum of the expect counts <5
  #or any zero gntp counts in affected or unaffected
  
  famsiz=table(table(famid))
  famsiz2=sum(famsiz[as.numeric(names(famsiz))>=2])#count # >2 sibs
  minfs2=10

  ### uses testiMAF to determine if gee or logistic regression should be used 
  
  if (famsiz2<=minfs2 || (testiMAFc+(bigN-bigNd)*imafc^2 <5) || (testiMAFd+bigNd*imafd^2<5) ){
     warning="logistic reg"	
     xcovar<-as.matrix(test.dat[na.ind,covar])
     if (bin.flag>2) {
        gee.test<-try(glm(phen1~x.snp+x.int+xcovar, family="binomial")) 
        if (!"try-error" %in% class(gee.test)) {
           gee.test.s<-summary(gee.test)
           gee.out <- c(cov.int,bigN,imaf,bigNd,imafd,gee.test.s$cov.scaled["x.snp","x.int"],"additive",gee.test.s$coef["x.snp",1:2],pchisq(gee.test.s$coef["x.snp",3]^2,1,lower.tail=F),
			  gee.test.s$coef["x.int",1:2],pchisq(gee.test.s$coef["x.int",3]^2,1,lower.tail=F),warning) #### add gee.test.s$cov.scaled["x.snp","x.int"] 091409
        } else gee.out <- c(cov.int,bigN,imaf,bigNd,imafd,rep(NA,8),"try-error")
     } else {
     if (sub=="Y") {  ####061209 sub
        tab1 <- table(xcovar[xcovar[,cov.int]==bin[1],cov.int],x.snp[xcovar[,cov.int]==bin[1]])
        tab2 <- table(xcovar[xcovar[,cov.int]==bin[2],cov.int],x.snp[xcovar[,cov.int]==bin[2]])
        if (length(tab1)==1 | length(tab2)==1) gee.out <- c(cov.int,bigN,imaf,bigNd,imafd,rep(NA,13),"non-informative SNP in subset") else {
           if (sum(!covar %in% cov.int)>0) { ####replace !covar %in% cov.int
              xcovar.bin<-as.matrix(test.dat[na.ind,covar[!covar %in% cov.int]])         
              gee.test1<-try(glm(phen1~x.snp+xcovar.bin, family="binomial", subset=test.dat[na.ind,cov.int]==bin[1])) 
              gee.test2<-try(glm(phen1~x.snp+xcovar.bin, family="binomial", subset=test.dat[na.ind,cov.int]==bin[2]))
              gee.test<-try(glm(phen1~x.snp+x.int+xcovar, family="binomial"))
           } else {
              gee.test1<-try(glm(phen1~x.snp, family="binomial", subset=test.dat[na.ind,cov.int]==bin[1])) 
              gee.test2<-try(glm(phen1~x.snp, family="binomial", subset=test.dat[na.ind,cov.int]==bin[2]))
              gee.test<-try(glm(phen1~x.snp+x.int+xcovar, family="binomial"))
           }                 
           if (sum(c(class(gee.test1),class(gee.test2),class(gee.test))=="try-error")==0) { ####061209 
              gee.test.s1<-summary(gee.test1);gee.test.s2<-summary(gee.test2);gee.test.s<-summary(gee.test);
              gee.out <- c(cov.int,bigN,imaf,bigNd,imafd,"additive",gee.test.s$coef["x.snp",1:2],pchisq(gee.test.s$coef["x.snp",3]^2,1,lower.tail=F),
                          gee.test.s1$coef["x.snp",1:2],pchisq(gee.test.s1$coef["x.snp",3]^2,1,lower.tail=F),
			     gee.test.s2$coef["x.snp",1:2],pchisq(gee.test.s2$coef["x.snp",3]^2,1,lower.tail=F),
                          gee.test.s$coef["x.int",1:2],pchisq(gee.test.s$coef["x.int",3]^2,1,lower.tail=F),warning)
           } else {
           if ("try-error" %in% c(class(gee.test1),class(gee.test2)) & !"try-error" %in% class(gee.test)) { ####061209 output interaction even subset analyses fail
              gee.test.s<-summary(gee.test);
              gee.out <- c(cov.int,bigN,imaf,bigNd,imafd,"additive",gee.test.s$coef["x.snp",1:2],pchisq(gee.test.s$coef["x.snp",3]^2,1,lower.tail=F),                          
                          rep(NA,6),gee.test.s$coef["x.int",1:2],pchisq(gee.test.s$coef["x.int",3]^2,1,lower.tail=F),warning)
           } else gee.out <- c(cov.int,bigN,imaf,bigNd,imafd,rep(NA,13),"try-error")
           }
        }
     } else { ####061209 sub="N" or NA
        gee.test<-try(glm(phen1~x.snp+x.int+xcovar, family="binomial")) 
        if (!"try-error" %in% class(gee.test)) {
           gee.test.s<-summary(gee.test)
           gee.out <- c(cov.int,bigN,imaf,bigNd,imafd,gee.test.s$cov.scaled["x.snp","x.int"],"additive",gee.test.s$coef["x.snp",1:2],pchisq(gee.test.s$coef["x.snp",3]^2,1,lower.tail=F),
			  gee.test.s$coef["x.int",1:2],pchisq(gee.test.s$coef["x.int",3]^2,1,lower.tail=F),warning)
        } else gee.out <- c(cov.int,bigN,imaf,bigNd,imafd,rep(NA,8),"try-error")
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
            if (gee.test.v$error!=0) warning="not converged" else warning=NA              
            gee.test <- summary(gee.test.v)
            gee.out <- c(cov.int,bigN,imaf,bigNd,imafd,gee.test.v$vbeta[names(gee.test.v$beta)=="x.snp",names(gee.test.v$beta)=="x.int"],"additive",
                        unlist(gee.test$mean["x.snp",1:2]),pchisq(gee.test$mean["x.snp",3],1,lower.tail=F),
                        unlist(gee.test$mean["x.int",1:2]),pchisq(gee.test$mean["x.int",3],1,lower.tail=F),warning)
         } else gee.out <- c(cov.int,bigN,imaf,bigNd,imafd,rep(NA,8),"try-error")
         if (gee.test.v$error!=0) gee.out= c(cov.int,bigN,imaf,bigNd,imafd,rep(NA,8),"not converged")
      } else {
      if (sub=="Y") { ####061209 sub
         tab1 <- table(xcovar[xcovar[,cov.int]==bin[1],cov.int],x.snp[xcovar[,cov.int]==bin[1]]) ######05132009 imputed.X
         tab2 <- table(xcovar[xcovar[,cov.int]==bin[2],cov.int],x.snp[xcovar[,cov.int]==bin[2]])
         if (length(tab1)==1 | length(tab2)==1) gee.out <- c(cov.int,bigN,imaf,bigNd,imafd,rep(NA,13),"non-informative SNP in subset") else { ######05132009 imputed.X
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
           if (sum(c(class(gee.test1),class(gee.test2),class(gee.test))=="try-error")==0) { ###061209
               if (gee.test1$error==0 & gee.test2$error==0 & gee.test$error==0) {
                  warning=NA                         
                  gee.out <- c(cov.int,bigN,imaf,bigNd,imafd,"additive",unlist(gee.test$mean["x.snp",1:2]),pchisq(gee.test$mean["x.snp",3],1,lower.tail=F),
                              unlist(gee.test1$mean["x.snp",1:2]),pchisq(gee.test1$mean["x.snp",3],1,lower.tail=F),unlist(gee.test2$mean["x.snp",1:2]),
                              pchisq(gee.test2$mean["x.snp",3],1,lower.tail=F),unlist(gee.test$mean["x.int",1:2]),pchisq(gee.test$mean["x.int",3],1,lower.tail=F),warning)
               } else 
               if (gee.test1$error+gee.test2$error!=0 & gee.test$error==0) {
                  warning="not converged only in at least one subset analysis"
                  gee.out <- c(cov.int,bigN,imaf,bigNd,imafd,"additive",unlist(gee.test$mean["x.snp",1:2]),pchisq(gee.test$mean["x.snp",3],1,lower.tail=F),rep(NA,6),
                             unlist(gee.test$mean["x.int",1:2]),pchisq(gee.test$mean["x.int",3],1,lower.tail=F),warning)
               } else {
                  warning="not converged" 
                  gee.out <- c(cov.int,bigN,imaf,bigNd,imafd,rep(NA,13),warning)               
               }
            } else {
            if ("try-error" %in% c(class(gee.test1),class(gee.test2)) & !"try-error" %in% class(gee.test)) { ####061209 output interaction even subset analyses fail
               if (gee.test$error!=0) {
                  warning="not converged" 
                  gee.out <- c(cov.int,bigN,imaf,bigNd,imafd,rep(NA,13),warning)
               } else {
                  warning=NA                         
                  gee.out <- c(cov.int,bigN,imaf,bigNd,imafd,"additive",unlist(gee.test$mean["x.snp",1:2]),pchisq(gee.test$mean["x.snp",3],1,lower.tail=F),rep(NA,6),
                             unlist(gee.test$mean["x.int",1:2]),pchisq(gee.test$mean["x.int",3],1,lower.tail=F),warning)
               }
            } else gee.out <- c(cov.int,bigN,imaf,bigNd,imafd,rep(NA,13),"try-error")
            }
         } 
      } else { ####061209 sub="N" or NA
         gee.test.v <- try(geese(phen1~x.snp+x.int+xcovar,id=famid,family="binomial",corstr="independence"))
         if (!"try-error" %in% class(gee.test.v)) {
            if (gee.test.v$error!=0) warning="not converged" else warning=NA              
            gee.test <- summary(gee.test.v)
            gee.out <- c(cov.int,bigN,imaf,bigNd,imafd,gee.test.v$vbeta[names(gee.test.v$beta)=="x.snp",names(gee.test.v$beta)=="x.int"],"additive",
                        unlist(gee.test$mean["x.snp",1:2]),pchisq(gee.test$mean["x.snp",3],1,lower.tail=F),unlist(gee.test$mean["x.int",1:2]),
                        pchisq(gee.test$mean["x.int",3],1,lower.tail=F),warning)
         } else gee.out <- c(cov.int,bigN,imaf,bigNd,imafd,rep(NA,8),"try-error")
         if (gee.test.v$error!=0) gee.out= c(cov.int,bigN,imaf,bigNd,imafd,rep(NA,8),"not converged")
      }
      }
      return(gee.out)
  }  
}
 