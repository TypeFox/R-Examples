geepack.quant.batch <- function(phenfile,genfile,pedfile,phen,model="a",covars=NULL,outfile,col.names=T,sep.ped=",",sep.phe=",",sep.gen=","){
###########################################################
  if (!model %in% c("g","r","a","d")) 
	stop('please specify model as "a","g","r" or "d" only')
  
  cor.snp <- function(y,x){  
   if(!is.numeric(y))y<-as.numeric(as.factor(y)) 
   return(sd(y)==0 || abs(cor(y,x,use="complete"))>0.99999999 )} 

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
  assign("phen",phen,pos=-1,inherits=T)
  phensnp.dat <- read.in.data(phenfile,genfile,pedfile,sep.ped=sep.ped,sep.phe=sep.phe,sep.gen=sep.gen)
  snplist<-phensnp.dat$snps

  if (is.null(covars)) phenlist<-phensnp.dat$phen.name else 
     if (!is.null(covars) & sum(phensnp.dat$phen.name %in% covars)==length(covars)) phenlist<-phensnp.dat$phen.name[!phensnp.dat$phen.name %in% covars] else  
        stop('some covariates are not available')

  test.dat <- phensnp.dat$data
  test.dat <- test.dat[order(test.dat$famid),]
  assign("test.dat", test.dat, pos=-1,inherits=T)

  if (!is.null(covars) & sum(snplist %in% covars)>=1) {
     names(test.dat)[which(names(test.dat) %in% paste(snplist[snplist %in% covars],".x",sep=""))] <- snplist[snplist %in% covars]
     covars[covars %in% snplist] <- paste(covars[covars %in% snplist],".y",sep="")
  }

  if (!is.null(covars)) {
     covars.dat <- na.omit(test.dat[,covars])
     single.cov <- F
     if (length(covars)==1) single.cov <- var(covars.dat)==0 else {
        single.cov <- any(apply(covars.dat,2,var)==0)
        if (single.cov) stop(paste("Single category in covariates!"))
        for (i in covars){
            cov1 <- covars.dat[,i]
            if (!is.numeric(cov1)) cov1 <-as.numeric(as.factor(cov1))
            for (j in covars[covars!=i]){
                cov2 <- covars.dat[,j]
                if (!is.numeric(cov2)) cov2 <-as.numeric(as.factor(cov2))
                if (abs(cor(cov1,cov2))>0.99999999) stop(paste("Highly correlated covariates ",i," and ",j,"!!",sep=""))
            }
        }
     }
  }   
  
  idlab <- "famid"
  result <- NULL
  famid <- NULL; rm(famid)

      for (i in snplist) {
          assign("i",i,pos=-1,inherits=T)
          if (is.null(covars)) test2.dat <- na.omit(test.dat[,c(i,phen,idlab)]) else { 
             test2.dat <- na.omit(test.dat[,c(i,phen,idlab,covars)])
             x.covar<-as.matrix(test2.dat[,covars])
             assign("x.covar", x.covar, pos=-1,inherits=T)  
          } 
          assign("test2.dat", test2.dat, pos=-1,inherits=T)
                            
          count<-table(test2.dat[,i])
          gntps<-names(count)
          count1<-rep(0,3) # if count of certain gntp=0, match counts with gntps in ouput
          count1[as.numeric(gntps)+1]<-count

          if (!is.null(covars) & length(count)==1) colinear <- F else
             if (!is.null(covars) & length(covars)>1 & length(count)>1) colinear <- apply(x.covar,2,cor.snp,x=test2.dat[,i]) else 
                if (!is.null(covars) & length(covars)==1 & length(count)>1) colinear <- cor.snp(x.covar,test2.dat[,i]) else 
                   if (is.null(covars)) colinear <- F 

          if (sum(colinear)>0) { 
             if (model %in% c("a","r","d")) result<-rbind(result,c(phen,i,count1,rep(NA,6))) else  result<-rbind(result,c(phen,i,count1,rep(NA,10))) 
             } else {

             if (sort(count1)[1]+sort(count1)[2]<1 | (count1[3]==0 & model=="r")){ #if count of less common genotypes <10, not run
                if (model %in% c("a","r","d")) result<-rbind(result,c(phen,i,count1,rep(NA,6))) else  result<-rbind(result,c(phen,i,count1,rep(NA,10)))
                } else {  		
                  if (model=="a"){
                      if (is.null(covars)) gee.test.v<-try(geese(test2.dat[,phen]~test2.dat[,i],data=test2.dat,id=famid,corstr="independence")) else
                                gee.test.v<-try(geese(test2.dat[,phen]~test2.dat[,i]+x.covar,data=test2.dat,id=famid,corstr="independence"))
                      if (class(gee.test.v)!="try-error") {
                         gee.test <- summary(gee.test.v)
                         if (class(gee.test)!="try-error" & gee.test$error==0) tmp<-c(as.numeric(unlist(gee.test$mean[2,1:3])),1,"additive",pchisq(gee.test$mean[2,"wald"],1,lower.tail=F)) else tmp<-rep(NA,6)
                      } else tmp<-rep(NA,6)
                  } else if (model=="g"){
	                    if (min(count1)<1){ #run dominant model
                               snp.i<-test2.dat[,i]; snp.i[snp.i==2]<-1; assign("snp.i",snp.i,pos=-1,inherits=T);
		                 if (is.null(covars)) gee.test.v<-try(geese(test2.dat[,phen]~snp.i,data=test2.dat,id=famid,corstr="independence")) else
                                                    gee.test.v<-try(geese(test2.dat[,phen]~snp.i+x.covar,data=test2.dat,id=famid,corstr="independence"))
                               if (class(gee.test.v)!="try-error") {
                                  gee.test <- summary(gee.test.v)
                                  if (class(gee.test)!="try-error" & gee.test$error==0) tmp<-c(gee.test$mean[2,1],NA,NA,gee.test$mean[2,2], 
                                         NA,NA,gee.test$mean[2,3],1,"dominant",pchisq(gee.test$mean[2,3],1,lower.tail=F)) else tmp<-rep(NA,10)
                               } else tmp<-rep(NA,10)
  	   	             } else {
                                x.snp <- factor(test2.dat[,i])
                                if (is.null(covars)) gee.test.v<-try(geese(test2.dat[,phen]~x.snp,data=test2.dat,id=famid,corstr="independence")) else
                                   gee.test.v<-try(geese(test2.dat[,phen]~x.snp+x.covar,data=test2.dat,id=famid,corstr="independence"))
                                if (class(gee.test.v)!="try-error") {
                                   gee.test <- summary(gee.test.v)
                                   if (class(gee.test.v)!="try-error" & gee.test$error==0) {
                                      chisq=try(t(matrix(gee.test$mean[c("x.snp1","x.snp2"),1],ncol=1))%*%solve(gee.test.v$vbeta[names(gee.test.v$beta)%in%c("x.snp1","x.snp2"),names(gee.test.v$beta)%in%c("x.snp1","x.snp2")])%*%matrix(gee.test$mean[c("x.snp1","x.snp2"),1],ncol=1))                                  
                                      tmp<-c(as.numeric(unlist(gee.test$mean[2:3,1])),gee.test$mean[3,1]-gee.test$mean[2,1],gee.test$mean["x.snp1",2],gee.test$mean["x.snp2",2],
                                         sqrt(gee.test.v$vbeta[names(gee.test.v$beta)=="x.snp1",names(gee.test.v$beta)=="x.snp1"]+gee.test.v$vbeta[names(gee.test.v$beta)=="x.snp2",names(gee.test.v$beta)=="x.snp2"]-2*gee.test.v$vbeta[names(gee.test.v$beta)=="x.snp1",names(gee.test.v$beta)=="x.snp2"]),
                                         chisq,2,"general",pchisq(chisq,2,lower.tail=F))
                                   } else tmp<-rep(NA,10)
                                } else tmp<-rep(NA,10)
                              }
	          } else if (model=="d"){#11 12 22 code as 0 1 1 
		              snp.i<-test2.dat[,i]; snp.i[snp.i==2]<-1; assign("snp.i",snp.i,pos=-1,inherits=T);
       	              if (is.null(covars)) gee.test.v<-try(geese(test2.dat[,phen]~snp.i,data=test2.dat,id=famid,corstr="independence")) else
                                                   gee.test.v<-try(geese(test2.dat[,phen]~snp.i+x.covar,data=test2.dat,id=famid,corstr="independence"))
                            if (class(gee.test.v)!="try-error") {
                               gee.test <- summary(gee.test.v)
                               if (class(gee.test)!="try-error" & gee.test$error==0) tmp<-c(as.numeric(unlist(gee.test$mean[2,1:3])),1,"dominant",pchisq(gee.test$mean[2,"wald"],1,lower.tail=F)) else tmp<-rep(NA,6)
                            } else tmp<-rep(NA,6)
	          } else if (model=="r"){#11 12 22 code as 0 0 1
             		       snp.i<-test2.dat[,i];snp.i[snp.i==1]<-0;snp.i[snp.i==2]<-1; assign("snp.i",snp.i,pos=-1,inherits=T);
		              if (is.null(covars)) gee.test.v<-try(geese(test2.dat[,phen]~snp.i,data=test2.dat,id=famid,corstr="independence")) else
                                                   gee.test.v<-try(geese(test2.dat[,phen]~snp.i+x.covar,data=test2.dat,id=famid,corstr="independence"))
                            if (class(gee.test.v)!="try-error") {
                               gee.test <- summary(gee.test.v)
                               if (class(gee.test)!="try-error" & gee.test$error==0) tmp<-c(as.numeric(unlist(gee.test$mean[2,1:3])),1,"recessive",pchisq(gee.test$mean[2,"wald"],1,lower.tail=F)) else tmp<-rep(NA,6)
                            } else tmp<-rep(NA,6)
                 }
                 result <- rbind(result,c(phen,i,count1,tmp)) 
             }
          } 
          
      }    

  if (model %in% c("a","d","r")){
  	colnames(result)<-c("phen","snp","n0","n1","n2","beta","se","chisq","df","model","pval")  
  }else colnames(result)<-c("phen","snp","n0","n1","n2","beta10","beta20","beta21","se10",
			"se20","se21","chisq","df","model","pval")  

  write.table(result, outfile, quote=F,row.names=F, col.names=col.names,sep=",",na="",append=T)

}

