geepack.lgst=function(snp,phen,test.dat,covar=NULL,model="a"){
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
  if(!is.null(covar)) x.covar=test.dat[na.ind,covar] 

  ###########################################
  #produce genotyp count
  count=table(snp1)
  gntps<-names(count)
  count1<-rep(0,3) # if count of certain gntp=0, match counts with gntps in ouput
  count1[as.numeric(gntps)+1]<-count

  #########################################  
  #check categories in y, if more than 2, stop the run
  length.unique<-length(unique(na.omit(phen1)))

  if (length.unique >2) 
	stop (paste("More than two categories in the phenotype: ",phen, 
                    ". Program stopped because of this error",sep="")) 

  if (length.unique<=1 ) {
   	if (model %in% c("a","d","r")) {      
        	gee.out= matrix(c(count1,rep(NA,11),"single category in phen",NA),ncol=1)
   	}else  
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
    #print("Fail 1.5")
   if (model %in% c("a","d","r")) {      
        gee.out= matrix(c(count1,count.d,rep(NA,10)),ncol=1)
   }else
        gee.out= matrix(c(count1,count.d,rep(NA,14)),ncol=1)

  return(gee.out)
 }
 ##two genotypes, minimum count <10,skip 
 if (length(count)==2 && min(count)<10) {
    #print("Fail 2")
   ##do not run gee if genotype count <10 for 2 gntp SNP, 
    if (model %in% c("a","d","r")) {      
        gee.out= matrix(c(count1,count.d,rep(NA,10)),ncol=1)
    }else
        gee.out= matrix(c(count1,count.d,rep(NA,14)),ncol=1)
  return(gee.out)
 }else if (length(count)==2 && !min(count)<10) {   
   if (model=="a") {x.snp = snp1     
   }else if (model=="g") {x.snp = ifelse(snp1!=2,snp1,1); 
   }else if (model=="d") {x.snp = ifelse(snp1!=2,snp1,1)
                          if (length(table(x.snp))==1) {
                             gee.out= matrix(c(count1,count.d,rep(NA,10)),ncol=1)
                             return(gee.out)
                          }
   }else if (model=="r" & !count1[3]<10) {x.snp=ifelse(snp1==2,1,0)
   }else if (model=="r" & count1[3]<10) {gee.out= matrix(c(count1,count.d,rep(NA,10)),ncol=1)
                                         return(gee.out)
   }   
}

 ##three genotypes, but the sum of the two lower counts < 10, skip
 if (length(count)==3 && (sort(count)[1]+sort(count)[2])<10) { 
   #print("Fail 3")
   if (model %in% c("a","d","r")) {      
        gee.out= matrix(c(count1,count.d,rep(NA,10)),ncol=1)
   }else
        gee.out= matrix(c(count1,count.d,rep(NA,14)),ncol=1)
  return(gee.out)
 }else if (length(count)==3 && min(count)<10) {   
           if (model %in% c("d","g")) {x.snp = ifelse(snp1!=2,snp1,1); 
		#only change model when model="a" because model="g" has more fields in the results
           }else if (model=="r") {
                    x.snp=ifelse(snp1==2,1,0)
                    if (min(table(x.snp))<10) { gee.out= matrix(c(count1,count.d,rep(NA,10)),ncol=1)
                                 return(gee.out)
                                              }
           }else if (model=="a") x.snp = snp1       
       }else if (length(count)==3 && !min(count)<10) {
                 if (model=="a") {x.snp = snp1   
                 }else if (model=="g") {x.snp = factor(snp1)
                 }else if (model=="d") {x.snp = ifelse(snp1!=2,snp1,1)
                 }else if (model=="r") {x.snp=ifelse(snp1==2,1,0)
                 }   
 }                      


 ## if covariates, check the collinearity, skip of collinearity
  #function to check collinarity

#############################################
######NEW cor.snp handles factor variables by making them linear (not perfect solution):
  cor.snp <- function(y,x){
   if(!is.numeric(y))y<-as.numeric(as.factor(y))
   return(sd(y)==0 || abs(cor(y,x,use="complete"))>0.99999999 )}
###############################################
  if(length(table(x.snp))>1 && !is.null(covar) && model!="g" | (length(table(x.snp))==2 && !is.null(covar) && model=="g")){
 
#############################################
# Fix to work if there is only one covariate:  
###############################################
    if(length(covar)>1)colinear <- apply(x.covar,2,cor.snp,x=x.snp) else colinear<-cor.snp(x.covar,x.snp)
    ##check colinearity
    if( sum(colinear)>0){ ## only 1 gntp | colinearity between covar and SNP
      	#print("Fail 1")
      	if (model %in% c("a","d","r")) {         
        	gee.out= matrix(c(count1,count.d,rep(NA,8),"collinarity",NA),ncol=1) 
      	}else
        	gee.out= matrix(c(count1,count.d,rep(NA,12),"collinarity",NA),ncol=1)
    return(gee.out)
    }
  }

  
### test differential missingness between cases and controls
  snp.na <- rep(1,length(snp))
  snp.na[is.na(snp)] <- 0
  if (length(table(snp.na))==2) miss.mat <- matrix(table(test.dat[,phen],snp.na),2,2) else miss.mat <- cbind(c(0,0),table(test.dat[,phen],snp.na))
  miss.diff.p <- fisher.test(miss.mat)$p.value
  if (miss.diff.p>1) miss.diff.p <- 1
  miss.0 <- miss.mat[1,1]/sum(miss.mat[1,])
  miss.1 <- miss.mat[2,1]/sum(miss.mat[2,])

  # Perform GLM if number of clusters with 2 or more individuals <minfs2
  #or if minimum of the expect counts <5
  #or any zero gntp counts in affected or unaffected
  famsiz=table(table(famid))
  famsiz2=sum(famsiz[as.numeric(names(famsiz))>=2])#count # >2 sibs
  minfs2=10
      
  tmp.table <- apply(table(test.dat[na.ind,phen],x.snp)>0,1,sum)      
  not.enough <- sum(tmp.table<length(table(x.snp))) > 0 #zero count of a gntp in affected/unaffected

############################################################ deal with binary covariates ######082608
  cell0 <- F
  if (!is.null(covar)){
     if (length(covar)==1) cat.covar <- length(unique(x.covar)) else cat.covar <- apply(x.covar,2,function(x)length(unique(x))) #####102208
     if (any(cat.covar==1)) {
        if (model %in% c("a","d","r")) gee.out= matrix(c(count1,count.d,rep(NA,8),"covariate",NA),ncol=1) else gee.out= matrix(c(count1,count.d,rep(NA,12),"covariate",NA),ncol=1)
        return(gee.out)
     }  else {
     if (any(cat.covar==2)) {
        if (model=="g" && length(table(x.snp))==3) cell0 <- F else {
           bin.cov <- covar[cat.covar==2]
           if (length(bin.cov)==1) {
              if (length(covar)==1) cell0 <- sum(table(x.covar,x.snp)[1,-1])==0 | sum(table(x.covar,x.snp)[2,-1])==0 else #####102208
                  cell0 <- sum(table(x.covar[,bin.cov],x.snp)[1,-1])==0 | sum(table(x.covar[,bin.cov],x.snp)[2,-1])==0 
           } else {
                cell0 <- any(apply(x.covar[,bin.cov],2,function(x,snp=x.snp)sum(table(x,snp)[1,-1])==0 | sum(table(x,snp)[2,-1])==0,snp=x.snp))
                for (i in bin.cov){
                   for (j in bin.cov[bin.cov!=i]){
               tab.bin.cov00 <- sum(table(phen1[x.covar[,i]==sort(unique(x.covar[,i]))[1] & x.covar[,j]==sort(unique(x.covar[,j]))[1]],
                  x.snp[x.covar[,i]==sort(unique(x.covar[,i]))[1] & x.covar[,j]==sort(unique(x.covar[,j]))[1]])[,-1])==0
               tab.bin.cov01 <- sum(table(phen1[x.covar[,i]==sort(unique(x.covar[,i]))[1] & x.covar[,j]==sort(unique(x.covar[,j]))[2]],
                  x.snp[x.covar[,i]==sort(unique(x.covar[,i]))[1] & x.covar[,j]==sort(unique(x.covar[,j]))[2]])[,-1])==0
               tab.bin.cov10 <- sum(table(phen1[x.covar[,i]==sort(unique(x.covar[,i]))[2] & x.covar[,j]==sort(unique(x.covar[,j]))[1]],
                  x.snp[x.covar[,i]==sort(unique(x.covar[,i]))[2] & x.covar[,j]==sort(unique(x.covar[,j]))[1]])[,-1])==0
               if (sum(x.snp[x.covar[,i]==sort(unique(x.covar[,i]))[2] & x.covar[,j]==sort(unique(x.covar[,j]))[2]])!=0)   #####042010 dealing exclusive indicators
               tab.bin.cov11 <- sum(table(phen1[x.covar[,i]==sort(unique(x.covar[,i]))[2] & x.covar[,j]==sort(unique(x.covar[,j]))[2]],
                  x.snp[x.covar[,i]==sort(unique(x.covar[,i]))[2] & x.covar[,j]==sort(unique(x.covar[,j]))[2]])[,-1])==0 else tab.bin.cov11 <- F
               #print(c(i,j,tab.bin.cov00,tab.bin.cov10,tab.bin.cov01,tab.bin.cov11))
               cell0 <- tab.bin.cov00 | tab.bin.cov01 | tab.bin.cov10 | tab.bin.cov11 | cell0
                                                 }
                                  }
                  }
           }     
        }
     }
  }
############################################################# 

  if (famsiz2<=minfs2 | cell0 ||not.enough ){ 
	if(!model %in% c("a")){  
		if(min(chisq.test(table(phen1,x.snp))$expected)<5)
			warning="logistic reg & exp count<5" else warning="logistic reg"

	}
	else{
		x.snp.d = ifelse(snp1!=2,snp1,1) #warning for additive model is the same as for dominant model
		if (min(chisq.test(table(phen1,x.snp.d))$expected)<5)
		   warning="logistic reg & exp count<5" else warning="logistic reg"
	}
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
              if (!"x.snp" %in% rownames(gee.test.s$coef)) {
                 if (model %in% c("a","d","r")) gee.out <- matrix(c(count1,count.d,rep(NA,8),"collinarity",NA),ncol=1) else 
                    gee.out= matrix(c(count1,count.d,rep(NA,12),"collinarity",NA),ncol=1)
              } else {
  		  if (model!="g" | (model=="g" && length(table(x.snp))==2)) {
     		  ## additive,domninant,recessive summary
   	 	 	if (model %in% c("a","d","r")) {                                      
				gee.out <- matrix(c(count1,count.d,miss.0,miss.1,miss.diff.p,gee.test.s$coef["x.snp",1:2],
				gee.test.s$coef["x.snp",3]^2,"1",model,
				warning,gee.test.s$coef["x.snp",4]),ncol=1)
  	  		}else{
                        	x.snp.pos <- substr(rownames(gee.test.s$coef),start=1,stop=5)=="x.snp"
				gee.out <- matrix(c(count1,count.d,miss.0,miss.1,miss.diff.p,gee.test.s$coef[x.snp.pos,1],
					rep(NA,2),gee.test.s$coef[x.snp.pos,2],rep(NA,2),
					gee.test.s$coef[x.snp.pos,3]^2,"1","d",
					warning,gee.test.s$coef[x.snp.pos,4]),ncol=1)
			}
		  } else {   #(model=="g") && 3 levels     			
			chisq=try(t(gee.test.s$coef[c("x.snp1","x.snp2"),1])%*%
				solve(gee.test.s$cov.scaled[c("x.snp1","x.snp2"),c("x.snp1","x.snp2")])%*%
					gee.test.s$coef[c("x.snp1","x.snp2"),1])
			if (class(chisq)!="try-error"){
				gee.out <- matrix(c(count1,count.d,miss.0,miss.1,miss.diff.p,gee.test.s$coef["x.snp1",1],
						gee.test.s$coef["x.snp2",1],
						gee.test.s$coef["x.snp1",1]-gee.test.s$coef["x.snp2",1],
						gee.test.s$coef["x.snp1",2],gee.test.s$coef["x.snp2",2],
						sqrt(gee.test.s$cov.scaled["x.snp1","x.snp1"]+gee.test.s$cov.scaled["x.snp2","x.snp2"]-2*gee.test.s$cov.scaled["x.snp1","x.snp2"]),
					chisq,"2",model,warning,pchisq(chisq,1,lower.tail=F)),ncol=1)
			}else gee.out=matrix(c(count1,count.d,rep(NA,11),"try-error",NA,NA),ncol=1)
                } 
  		}
	} else {
		if (model %in% c("a","d","r")) {   
        		gee.out= matrix(c(count1,count.d,rep(NA,8),"try-error",NA),ncol=1) 
      		}else
        		gee.out= matrix(c(count1,count.d,rep(NA,12),"try-error",NA),ncol=1)

	}
	return(gee.out)
  }

##############################################################################################
#
# Perform GEE
#
 ################################change 4 ends here
  if(model !="g" || (model=="g" & length(table(x.snp))==2)) {   #########################change 5 starts and ends here
	if(!is.null(covar)){
#####################################
# Change how this is specified so that it works with factor covars:
               x.covar<-as.data.frame(test.dat[na.ind,covar])
               dat1<-cbind(phen1,x.covar,x.snp)
               form<-as.formula(paste("phen1~",paste(colnames(dat1)[2:length(colnames(dat1))],sep="+",collapse="+")))
               
		gee.test.v<-try(geese(form,data=dat1,id=famid,family="binomial",na.action=na.omit,corstr="independence"))

#####################################
   	} else {
     		gee.test.v <- try(geese(phen1 ~ x.snp, id=famid, family="binomial", na.action=na.omit,corstr="independence"))
  	} ## end if/else(!is.null(covar))
 	if(!"try-error" %in% class(gee.test.v)){ ######changed from if (class(gee.test) != "try-error") because class(gee.test) could be "gee" "glm", which gives warning
              gee.test <- summary(gee.test.v)                                              
 		if (gee.test$error!=0){
			if (!model %in% c("a")){     
				if(min(chisq.test(table(phen1,x.snp))$expected)<5)
              			warning="not converged and exp count<5" else warning="not converged"
			}
			else {
				x.snp.d = ifelse(snp1!=2,snp1,1) #warning for additive model is the same as for dominant model
                		if (min(chisq.test(table(phen1,x.snp.d))$expected)<5)
				warning="not converged and exp count<5" else warning="not converged"
			}
		}else{
			if (!model %in% c("a")){       
                        	if(min(chisq.test(table(phen1,x.snp))$expected)<5)
	                        warning="exp count<5" else warning=NA
			}
			else {
				x.snp.d = ifelse(snp1!=2,snp1,1)
				if (min(chisq.test(table(phen1,x.snp.d))$expected)<5)
				warning="exp count<5" else warning=NA
			}

		}
              x.snp.pos <- substr(rownames(gee.test$mean),start=1,stop=5)=="x.snp"
              chisq=try(gee.test$mean[x.snp.pos,"wald"])
  
  		##summarize output
		if(class(chisq)!="try-error"){
  			if (model !="g" ) {                            				
				gee.out <- matrix(c(count1,count.d,miss.0,miss.1,miss.diff.p,gee.test$mean["x.snp",1:3],"1",model,
					warning,pchisq(chisq,1,lower.tail=F)),ncol=1)
   
  			}else{
                		gee.out <- matrix(c(count1,count.d,miss.0,miss.1,miss.diff.p,
					gee.test$mean["x.snp",1],rep(NA,2),         #####################change 6 starts here
                		gee.test$mean["x.snp",2],rep(NA,2),chisq,"1","d",warning,   ################change 6 ends here
                    			pchisq(chisq,1,lower.tail=F)),ncol=1)
  			}
		}else {##try-error chisq
			if (model %in% c("a","d","r")) {   
                		gee.out= matrix(c(count1,count.d,rep(NA,8),"try-error",NA),ncol=1)
        		}else  
                		gee.out= matrix(c(count1,count.d,rep(NA,12),"try-error",NA),ncol=1)
		}
	}else {##try-error gee.test
		if (model %in% c("a","d","r")) {   
                	gee.out= matrix(c(count1,count.d,rep(NA,8),"try-error",NA),ncol=1)
        	}else  
                	gee.out= matrix(c(count1,count.d,rep(NA,12),"try-error",NA),ncol=1)
	}
       if (gee.test$error!=0) {
          if (model %in% c("a","d","r","fa")) gee.out= matrix(c(count1,count.d,rep(NA,8),"not converged",NA),ncol=1) else
             gee.out= matrix(c(count1,count.d,rep(NA,12),"not converged",NA),ncol=1)
       }

	return(gee.out)
  } else {#model="g" & 3 levels
   	if(!is.null(covar)){
     		x.covar=as.matrix(test.dat[na.ind,covar])
       	gee.test.v <- try(geese(phen1 ~ x.covar + x.snp, id=famid, family="binomial", na.action=na.omit,corstr="independence"))
   	} else {
		gee.test.v <- try(geese(phen1 ~ x.snp, id=famid, family="binomial", na.action=na.omit,corstr="independence"))
   	}
        if(!"try-error" %in% class(gee.test.v)){ ######changed from if (class(gee.test) != "try-error")
              gee.test <- summary(gee.test.v)
        	if (gee.test$error!=0){
                	if(min(chisq.test(table(phen1,x.snp))$expected)<5)
                        warning="not converged and exp count<5" else warning="not converged"
        	}else{
                	if(min(chisq.test(table(phen1,x.snp))$expected)<5)
                        warning="exp count<5" else warning=NA

        	}
		chisq=try(t(matrix(gee.test$mean[c("x.snp1","x.snp2"),1],ncol=1))%*%
			solve(gee.test.v$vbeta[names(gee.test.v$beta)%in%c("x.snp1","x.snp2"),names(gee.test.v$beta)%in%c("x.snp1","x.snp2")])%*%
			matrix(gee.test$mean[c("x.snp1","x.snp2"),1],ncol=1))
   		
		
		if(class(chisq)!="try-error"){
			gee.out=matrix(c(count1,count.d,miss.0,miss.1,miss.diff.p,gee.test$mean["x.snp1",1],gee.test$mean["x.snp2",1],
			gee.test$mean["x.snp2",1]-gee.test$mean["x.snp1",1],
			gee.test$mean["x.snp1",2],
			gee.test$mean["x.snp2",2],
			sqrt(gee.test.v$vbeta[names(gee.test.v$beta)=="x.snp1",names(gee.test.v$beta)=="x.snp1"]
			+gee.test.v$vbeta[names(gee.test.v$beta)=="x.snp2",names(gee.test.v$beta)=="x.snp2"]
			-2*gee.test.v$vbeta[names(gee.test.v$beta)=="x.snp1",names(gee.test.v$beta)=="x.snp2"]),
			chisq,"2",model,warning,pchisq(chisq,2,lower.tail=F)), ncol=1)
		}else 
         		gee.out= matrix(c(count1,count.d,rep(NA,12),"try-error",NA),ncol=1)
	}else
		gee.out=matrix(c(count1,count.d,rep(NA,12),"try-error",NA),ncol=1)
  if (gee.test$error!=0) gee.out= matrix(c(count1,count.d,rep(NA,12),"not converged",NA),ncol=1)
  } 
  return(gee.out)
}



