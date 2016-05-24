
### ********************************************** ###
coxph.ped <- function(phenfile,
                      phen,covars=NULL,
                      mafRange=c(0,0.05),
                      chr,
                      genfile,
                      pedfile,
                      snpinfoRdata,
                      sep.ped=",",
                      sep.phe=",",
                      sep.gen=" ",
                      time, 
                      aggregateBy="SKATgene",
                      maf.file,
                      snp.cor,
                      ssq.beta.wts=c(1,25),
                      singleSNP.outfile=F){  
  
  ### ********************************************** ###
  read.in.data <- function(phenfile,
                           genfile,
                           pedfile,
                           sep.ped = sep.ped,
                           sep.phe = sep.phe,
                           sep.gen = sep.gen,
                           snpinfo) {
    ### ********************************************** ###
    
    print("Reading in Data")
    
    ped.dat <- read.table(genfile,header=TRUE,na.strings="NA",sep=sep.gen, stringsAsFactors=FALSE, colClasses="integer")
    snp.names <- scan(genfile,nlines=1,sep=sep.gen,what="")[-1]
    snp.names <- snp.names[snp.names%in%snpinfo$Name] 
    colnames(ped.dat)[-1] <- snp.names
    
    pedigree <- read.table(pedfile,header=TRUE,sep=sep.ped,stringsAsFactors=FALSE, colClasses="integer" )
    colnames(pedigree)[colnames(pedigree)=="sex"] <- "sex.in.ped" 
    gntp.all <- merge(pedigree[,c("famid","id","sex.in.ped")],ped.dat,by="id") 
    
    phen.dat <- read.table(phenfile,header=TRUE,sep=sep.phe)  
    phen.name=colnames(phen.dat)[-1]
    
    phensnp.dat<-merge(gntp.all,phen.dat,by=c("id"))
    print("Done reading in data")
    return(list(data=phensnp.dat,snps=snp.names,phen.name=phen.name))
  }
  
  #library(survival); library(MASS); library(Matrix)
    
  print(paste("phenotype data = ", phenfile))
  print(paste("genotype data  = ", genfile ))
  print(paste("pedigree data  = ", pedfile ))
  
  if (is.null(covars)) print("Covariates = NONE") else print(paste("Covariate(s) =",covars,collapse=", "))
  
  print("Running COXPH")
  loadsnpinfo <- try(load(snpinfoRdata))  
  if (inherits(loadsnpinfo,"try-error")) stop(paste('SNP info Rdata does not exist at ',snpinfoRdata))
  snpinfo <- snpinfo[snpinfo$Chr==chr,]
  
  rsnpsingene.cor <- NULL
  loadsnpcor <- try(load(snp.cor))  
  if (inherits(loadsnpcor,"try-error")) stop(paste('SNP correlation matrix Rdata does not exist at ',snp.cor))
  
  maf <- read.csv(maf.file, stringsAsFactors=FALSE, colClasses=c("character","numeric") )
  
  snpinfo <- merge(snpinfo,maf,by="Name",all.x=T)  
  
  #### read in data
  phensnp.dat <- read.in.data(phenfile=phenfile,
                              genfile=genfile,
                              pedfile=pedfile,
                              sep.ped=sep.ped,
                              sep.phe=sep.phe,
                              sep.gen=sep.gen,
                              snpinfo=snpinfo)
  
  snplist<-phensnp.dat$snps
  if (!any(snplist%in%snpinfo$Name)) stop(paste('SNP name from genotype file do not match Name in SNP info Rdata!!'))
  test.dat<-phensnp.dat$data
    
  if (!is.null(covars) & sum(snplist %in% covars)>=1) {
    names(test.dat)[which(names(test.dat) %in% paste(snplist[snplist %in% covars],".x",sep=""))] <- snplist[snplist %in% covars]
    covars[covars %in% snplist] <- paste(covars[covars %in% snplist],".y",sep="")
  }
  test.dat<-test.dat[order(test.dat$famid),]
  
  #### collinearity checking
  if (!is.null(covars)) {
    
    covars.dat <- na.omit(test.dat[,covars])
    single.cov <- F
    
    if (length(covars)==1) {
      single.cov <- var(covars.dat)==0 
    }  else {
      single.cov <- any(apply(covars.dat,2,var)==0)
      if (single.cov) stop(paste("Single category in covariates!"))
      
      for (i in covars){
        cov1 <- covars.dat[,i]
               
        for (j in covars[covars!=i]){
          cov2 <- covars.dat[,j]
          if (abs(cor(cov1,cov2))>0.99999999) stop(paste("Highly correlated covariates ",i," and ",j,"!!",sep=""))
        } #end of j loop
      } # end of i loop
    } # end of if (length(covars)==1) 
  } # end of if (!is.null(covars))  
  
  
  snplist.len <- length(snplist)
  
  ### single variant analysis
  if (singleSNP.outfile==F) { 
    
    temp.out <- data.frame(ntotal = vector(length=snplist.len, mode="numeric"),
                           nmiss = vector(length=snplist.len, mode="numeric"),
                           maf_ntotal = vector(length=snplist.len, mode="numeric"),
                           beta = vector(length=snplist.len, mode="numeric"),
                           se = vector(length=snplist.len, mode="numeric"),
                           Z = vector(length=snplist.len, mode="numeric"),
                           remark = vector(length=snplist.len, mode="character"),
                           p = vector(length=snplist.len, mode="numeric"),
                           MAC = vector(length=snplist.len, mode="numeric"),
                           n0 = vector(length=snplist.len, mode="numeric"),
                           n1 = vector(length=snplist.len, mode="numeric"),
                           n2 = vector(length=snplist.len, mode="numeric"),
                           stringsAsFactors=FALSE) 
    
    if (is.null(covars)){
           
      for (i in 1:snplist.len) temp.out[i,] <- coxph.EC(snp=test.dat[,snplist[i]],
                                                        phen=phen,
                                                        test.dat=test.dat,
                                                        chr=chr,
                                                        time=time)  
      
    } else {

      for (i in 1:snplist.len) temp.out[i,] <- coxph.EC(snp=test.dat[,snplist[i]],
                                                        phen=phen,
                                                        test.dat=test.dat,
                                                        chr=chr,
                                                        covar=covars,
                                                        time=time)    
    } # end of  else for if (is.null(covars))
    
    temp.out$Name <- snplist           
    temp.out <- merge(temp.out,snpinfo[,c("Name","maf",aggregateBy)],by="Name",sort=F)
    temp.out <- temp.out[,c(aggregateBy,"Name","maf","ntotal","nmiss","maf_ntotal","beta","se","Z","remark","p","MAC","n0","n1","n2")]
    colnames(temp.out)[1] <- "gene"
    
    single.outfile <- paste(phen,"_singleSNP_chr_",chr,".txt",sep="")
    write.table(as.matrix(temp.out),single.outfile,quote=F,row.names=F,sep="\t",na="NA")
  } else {  
     temp.out <- read.table(paste(phen,"_singleSNP_chr_",chr,".txt",sep=""),header=T,sep="\t",as.is=T)  
     temp.out <- temp.out[temp.out$Name%in%snpinfo$Name,]
  }  
  
  ### Rdata & Burden test
  rdata.outfile <- paste(phen,"_chr_",chr,".RData",sep="")   
  gene.list <- na.omit(unique(snpinfo[snpinfo$Name%in%snplist,aggregateBy]))  
    
  gene.ssq.out  <- data.frame(gene=vector(mode="character",length=length(gene.list)),
                              SSQ=vector(mode="numeric", length=length(gene.list)),
                              cmafTotal=vector(mode="numeric", length=length(gene.list)),
                              cmafUsed=vector(mode="numeric", length=length(gene.list)),
                              nsnpsTotal=vector(mode="numeric", length=length(gene.list)),
                              nsnpsUsed=vector(mode="numeric", length=length(gene.list)),
                              nmiss=vector(mode="numeric", length=length(gene.list)),
                              df=vector(mode="numeric", length=length(gene.list)),
                              p=vector(mode="numeric", length=length(gene.list)),
                              stringsAsFactors=FALSE)
  
  gene.counts.out  <- data.frame(gene = vector(mode="character",length=length(gene.list)),
                                 beta = vector(mode="numeric", length=length(gene.list)),
                                 se = vector(mode="numeric", length=length(gene.list)),
                                 Z = vector(mode="numeric", length=length(gene.list)),
                                 cmafTotal = vector(mode="numeric", length=length(gene.list)),
                                 cmafUsed = vector(mode="numeric", length=length(gene.list)),
                                 nsnpsTotal = vector(mode="numeric", length=length(gene.list)),
                                 nsnpsUsed = vector(mode="numeric", length=length(gene.list)),
                                 nmiss = vector(mode="numeric", length=length(gene.list)),
                                 remark = vector(mode="numeric", length=length(gene.list)),
                                 p = vector(mode="numeric", length=length(gene.list)),
                                 stringsAsFactors=FALSE)
  
  gene.mbweights.out  <- data.frame(gene = vector(mode="character",length=length(gene.list)),
                                    beta = vector(mode="numeric", length=length(gene.list)),
                                    se = vector(mode="numeric", length=length(gene.list)),
                                    Z = vector(mode="numeric", length=length(gene.list)),
                                    cmafTotal = vector(mode="numeric", length=length(gene.list)),
                                    cmafUsed = vector(mode="numeric", length=length(gene.list)),
                                    nsnpsTotal = vector(mode="numeric", length=length(gene.list)),
                                    nsnpsUsed = vector(mode="numeric", length=length(gene.list)),
                                    nmiss = vector(mode="numeric", length=length(gene.list)),
                                    remark = vector(mode="numeric", length=length(gene.list)),
                                    p = vector(mode="numeric", length=length(gene.list)),
                                    stringsAsFactors=FALSE)
   
  Rdata <- rep(list(NA),length(gene.list))
  names(Rdata) <- gene.list
  
  is.na.p <- is.na(temp.out$p)
  
  for (i in 1:length(gene.list)){
    gene <- gene.list[i]; 
    
    if (singleSNP.outfile==F) { 
      
      snps.in.temp <- as.character(temp.out[temp.out$gene==gene,"Name"])
                                  
      na.snps.in.temp <- as.character(temp.out[temp.out$gene==gene & is.na.p,"Name"])   
      p.snps.in.temp <- as.character(temp.out[temp.out$gene==gene & !is.na.p,"Name"]) 
      
      scores <- sign(temp.out[temp.out$gene==gene,"beta"]) * sqrt(temp.out[temp.out$gene==gene,"Z"]^2)/(temp.out[temp.out$gene==gene,"se"])   
      scores[is.na(scores)] <- 0   
      icov <- matrix(0,length(scores),length(scores))   
      rownames(icov) <- colnames(icov) <- snps.in.temp   
      
      if (sum(temp.out$gene==gene & !is.na.p) >1 ) {
        
        diag.temp <- diag(1/temp.out[temp.out$gene==gene & !is.na.p,"se"])
        icov.base <- diag.temp %*% rsnpsingene.cor[[gene]][p.snps.in.temp,p.snps.in.temp] %*% diag.temp  
        rownames(icov.base) <- colnames(icov.base) <- p.snps.in.temp   
        icov[rownames(icov.base),colnames(icov.base)] <- icov.base[rownames(icov.base),colnames(icov.base)]
      } else {
        
        if (sum(temp.out$gene==gene & !is.na.p) == 1) {
          
          icov.base <- matrix(1/(temp.out[temp.out$gene==gene & !is.na.p,"se"]^2),1,1)
          rownames(icov.base) <- colnames(icov.base) <- p.snps.in.temp   
          icov[rownames(icov.base),colnames(icov.base)] <- icov.base[rownames(icov.base),colnames(icov.base)]
        }
      }
      
      n <- max(temp.out[temp.out$gene==gene,"ntotal"],na.rm=T) 
      af <- temp.out[temp.out$gene==gene,"maf_ntotal"]
      af[is.na(af)] <- 0
      Rdata[[i]] <- list(scores=scores,cov=Matrix(icov,sparse=T),n=n,maf=af,sey=1)
    }  
    
    snps.total <- snpinfo[snpinfo[,aggregateBy]==gene & snpinfo$Name%in%snplist,"Name"]
    cmafTotal <- sum(temp.out[temp.out$Name%in%snps.total,"maf_ntotal"],na.rm=T)
    nsnpsTotal <- sum(snpinfo[,aggregateBy]==gene & snpinfo$Name%in%snplist)   
    rsnps.in.gene <- as.character(temp.out[temp.out[,"Name"]%in%snps.total & ((!is.na(temp.out[,"maf_ntotal"])) & temp.out[,"maf_ntotal"]<max(mafRange) & temp.out[,"maf_ntotal"]>min(mafRange)),"Name"])   
    rsnps.maf <- temp.out[temp.out$Name%in%rsnps.in.gene,c("Name","maf_ntotal")]
    cmafUsed <- sum(rsnps.maf$maf_ntotal)   
    nsnpsUsed <- length(rsnps.in.gene)
    nmiss <- sum(temp.out[temp.out$Name %in% rsnps.in.gene,"nmiss"])
       
    #### ssq test
    if (length(rsnps.in.gene)==0) {
      
      gene.ssq.out[i,] <- list(gene,NA,NA,NA,NA,NA,NA,NA,NA)
    } else {       
      
      if (length(rsnps.in.gene)==1) {                     
        ssq.out <- temp.out[temp.out$Name==rsnps.in.gene,]
        gene.ssq.out[i,] <- list(gene, ssq.out[9]^2, cmafTotal, cmafUsed, nsnpsTotal, nsnpsUsed, nmiss, 1,ssq.out[11])
      } else {
        z.1 <- temp.out[temp.out$Name%in%rsnps.in.gene,c("Name","Z")]
        z.1 <- merge(z.1,rsnps.maf,by="Name",sort=F)
        wi <- dbeta(z.1$maf,ssq.beta.wts[1],ssq.beta.wts[2])^2
        ssq <- sum(z.1$Z ^2 * wi)
        
        rsnps.in.gene.cor <- rsnpsingene.cor[[gene]][rsnps.in.gene,rsnps.in.gene]
        
        diag.temp<-diag(sqrt(wi))  
        c <- eigen(diag.temp %*% rsnps.in.gene.cor %*% diag.temp)$values  
        a <- sum(c^3)/sum(c^2)
        b <- sum(c)-((sum(c^2))^2)/sum(c^3)
        d <- ((sum(c^2))^3)/((sum(c^3))^2)        
        p <- try(pchisq((ssq-b)/a,df=d,lower.tail=FALSE)) 

        if (!"try-error"%in%class(p)) gene.ssq.out[i,] <- list(gene,ssq,cmafTotal,cmafUsed,nsnpsTotal,nsnpsUsed,nmiss,d,p) else 
           gene.ssq.out[i,] <- list(gene,NA,NA,NA,NA,NA,NA,NA,NA)        
      } 
    }
    
    ### Burden tests      
    if (length(rsnps.in.gene)==0) {
      
      gene.counts.out[i,] <- list(gene,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
      gene.mbweights.out[i,] <- list(gene,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA) 
    } else {
      if (is.null(covars)) {
        rsnps.dat <- test.dat[,c("famid","id","sex.in.ped",rsnps.in.gene,phen,time)]

        ### T count test
        if (length(rsnps.in.gene)==1){
          counts.out <- temp.out[temp.out$Name==rsnps.in.gene,-(1:3)]
        } else {   
                  rsnps.dat$counts <- rowSums(rsnps.dat[,rsnps.in.gene],na.rm=TRUE)                        
                  counts.out <- coxph.EC(snp=rsnps.dat$counts,phen=phen,test.dat=rsnps.dat,chr=chr,time=time)   
        } 
        
        
        gene.counts.out[i,] <- list(gene,
                                    counts.out[4],counts.out[5],counts.out[6],
                                    cmafTotal,cmafUsed,nsnpsTotal,nsnpsUsed,nmiss,counts.out[7],counts.out[8])
        
        
        ### MB weight test
        for (rsnp in rsnps.in.gene) rsnps.dat[,rsnp] <- rsnps.dat[,rsnp]/(1-rsnps.maf[rsnps.maf$Name==rsnp,"maf_ntotal"])/rsnps.maf[rsnps.maf$Name==rsnp,"maf_ntotal"]
        
        if (length(rsnps.in.gene)==1) rsnps.dat$mbweights <- rsnps.dat[,rsnps.in.gene] else rsnps.dat$mbweights <- rowSums(rsnps.dat[,rsnps.in.gene],na.rm=TRUE)      
        mbweights.out <- coxph.EC(snp=rsnps.dat$mbweights,phen=phen,test.dat=rsnps.dat,chr=chr,time=time)        
        
        gene.mbweights.out[i,]  <- list(gene,
                                        mbweights.out[4],mbweights.out[5],mbweights.out[6],
                                        cmafTotal,cmafUsed,nsnpsTotal,nsnpsUsed,nmiss,
                                        mbweights.out[7] ,mbweights.out[8] )      
        
      } else {
        rsnps.dat <- test.dat[,c("famid","id","sex.in.ped",rsnps.in.gene,phen,covars,time)]  

        ### T count test       
        if (length(rsnps.in.gene)==1) counts.out <- temp.out[temp.out$Name==rsnps.in.gene,-(1:3)] else {   
          
          rsnps.dat$counts <- rowSums(rsnps.dat[,rsnps.in.gene], na.rm=TRUE)        
          counts.out <- coxph.EC(snp=rsnps.dat$counts,phen=phen,test.dat=rsnps.dat,covar=covars,chr=chr,time=time)           
        }    
        
        gene.counts.out[i,] <- list(gene,
                                    counts.out[4],counts.out[5],counts.out[6],
                                    cmafTotal,cmafUsed,nsnpsTotal,nsnpsUsed,nmiss,
                                    counts.out[7],counts.out[8])
        
        ### MB weight test
        for (rsnp in rsnps.in.gene) rsnps.dat[,rsnp] <- rsnps.dat[,rsnp]/(1-rsnps.maf[rsnps.maf$Name==rsnp,"maf_ntotal"])/rsnps.maf[rsnps.maf$Name==rsnp,"maf_ntotal"]
        
        if (length(rsnps.in.gene)==1) rsnps.dat$mbweights <- rsnps.dat[,rsnps.in.gene] else rsnps.dat$mbweights <- rowSums(rsnps.dat[,rsnps.in.gene],na.rm=TRUE)       
        mbweights.out <- coxph.EC(snp=rsnps.dat$mbweights,phen=phen,test.dat=rsnps.dat,covar=covars,chr=chr,time=time)
        
        gene.mbweights.out[i,] <- list(gene, 
                                       mbweights.out[4],mbweights.out[5],mbweights.out[6],
                                       cmafTotal,cmafUsed,nsnpsTotal,nsnpsUsed,nmiss,
                                       mbweights.out[7], mbweights.out[8])        
        
      }
    }
  }
  
  Tcounts.outfile <- paste(phen,"_T",max(mafRange),"_chr_",chr,".txt",sep="")  
  MBweights.outfile <- paste(phen,"_MB",max(mafRange),"_chr_",chr,".txt",sep="")  
  SSQ.outfile <- paste(phen,"_SSQ",max(mafRange),"_chr_",chr,".txt",sep="")  
  save(Rdata,file=rdata.outfile)  
  write.table(as.matrix(gene.counts.out),Tcounts.outfile,quote=F,row.names=F,sep="\t",na="NA")
  write.table(as.matrix(gene.mbweights.out),MBweights.outfile,quote=F,row.names=F,sep="\t",na="NA")
  write.table(as.matrix(gene.ssq.out),SSQ.outfile,quote=F,row.names=F,sep="\t",na="NA")
  if (singleSNP.outfile==F) save(Rdata,file=rdata.outfile)  
}