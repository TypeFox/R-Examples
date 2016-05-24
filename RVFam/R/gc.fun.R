
gc.fun <- function(path,phen,snpinfoRdata,snp.cor,mac,aggregateBy="SKATgene",maf.file,mafRange,ssq.beta.wts=c(1,25)){
  #library(MASS); library(Matrix)

  if (missing(mac)) stop(paste('Please input MAC threshold for GC correction!'))   
  loadsnpinfo <- try(load(snpinfoRdata))  
  if (inherits(loadsnpinfo,"try-error")) stop(paste('SNP info Rdata does not exist at ',snpinfoRdata))  

  rsnpsingene.cor <- NULL
  loadsnpcor <- try(load(snp.cor))  
  if (inherits(loadsnpcor,"try-error")) stop(paste('SNP correlation matrix Rdata does not exist at ',snp.cor))  

  maf <- read.csv(maf.file, stringsAsFactors=FALSE, colClasses=c("character","numeric") ) 
  snpinfo <- merge(snpinfo,maf,by.x="Name",all.x=T)

  setwd(path)
  temp.out <- NULL
  if (!all(paste(phen,"_singleSNP_chr_",c(1:22,"X"),".txt",sep="")%in%list.files(path=path,full.names=F))) stop(paste("Not all 23 chromosome single SNP result files available in ",path,"!!",sep="")) else
     for (chr in c(1:22,"X")) temp.out <- rbind(temp.out,read.table(paste(phen,"_singleSNP_chr_",chr,".txt",sep=""),header=T,sep="\t",as.is=T))
 
  lambda <- qchisq(median(na.omit(temp.out$p[!is.na(temp.out$MAC) & temp.out$MAC<mac])),df=1,lower.tail=F)/qchisq(0.5,1)
  temp.out$Z.gc <- temp.out$Z
  temp.out$p.gc <- temp.out$p
  if (lambda>1) {
     temp.out$Z.gc[!is.na(temp.out$MAC) & temp.out$MAC<mac] <- temp.out$Z[!is.na(temp.out$MAC) & temp.out$MAC<mac]/sqrt(lambda) 
     temp.out$p.gc[!is.na(temp.out$MAC) & temp.out$MAC<mac] <- pchisq((temp.out$Z.gc[!is.na(temp.out$MAC) & temp.out$MAC<mac])^2,1, lower.tail=F)
  } else stop("Lambda <= 1, no need for GC correction")

  write.table(temp.out,paste(phen,"_singleSNP_GC.txt",sep=""),quote=F,row.names=F,sep="\t",na="NA")

  gene.list <- na.omit(unique(snpinfo[snpinfo$Name%in%temp.out$Name,aggregateBy]))  
  Rdata <- rep(list(NA),length(gene.list))
  names(Rdata) <- gene.list

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

  is.na.p <- is.na(temp.out$p)

  for (i in 1:length(gene.list)){
      gene <- gene.list[i]; 
      snps.in.temp <- as.character(temp.out[temp.out$gene==gene,"Name"])   
      na.snps.in.temp <- as.character(temp.out[temp.out$gene==gene & is.na.p,"Name"])   
      p.snps.in.temp <- as.character(temp.out[temp.out$gene==gene & !is.na.p,"Name"])  
      scores <- sign(temp.out[temp.out$gene==gene,"beta"])*sqrt(temp.out[temp.out$gene==gene,"Z"]^2)/(temp.out[temp.out$gene==gene,"se"])/lambda  
      scores[is.na(scores)] <- 0   
      icov <- matrix(0,length(scores),length(scores))   
      rownames(icov) <- colnames(icov) <- snps.in.temp   

      if (sum(temp.out$gene==gene & !is.na.p)>1) {
         diag.temp <- diag(1/temp.out[temp.out$gene==gene & !is.na.p,"se"])
         icov.base <- diag.temp %*% rsnpsingene.cor[[gene]][p.snps.in.temp,p.snps.in.temp] %*% diag.temp/lambda
         rownames(icov.base) <- colnames(icov.base) <- p.snps.in.temp   
         icov[rownames(icov.base),colnames(icov.base)] <- icov.base[rownames(icov.base),colnames(icov.base)]
      } else 
      if (sum(temp.out$gene==gene & !is.na.p)==1) {
         icov.base <- matrix(1/(temp.out[temp.out$gene==gene & !is.na.p,"se"]^2)/lambda,1,1)
         rownames(icov.base) <- colnames(icov.base) <- p.snps.in.temp   
         icov[rownames(icov.base),colnames(icov.base)] <- icov.base[rownames(icov.base),colnames(icov.base)]
      }
      n <- max(temp.out[temp.out$gene==gene,"ntotal"],na.rm=T) 
      af <- temp.out[temp.out$gene==gene,"maf_ntotal"]
      af[is.na(af)] <- 0
      Rdata[[i]] <- list(scores=scores,cov=Matrix(icov,sparse=T),n=n,maf=af,sey=1)

      snps.total <- snpinfo[snpinfo[,aggregateBy]==gene & snpinfo$Name%in%temp.out$Name,"Name"]
      cmafTotal <- sum(temp.out[temp.out$Name%in%snps.total,"maf_ntotal"],na.rm=T)
      nsnpsTotal <- sum(snpinfo[,aggregateBy]==gene & snpinfo$Name%in%temp.out$Name)
                
      rsnps.in.gene <- as.character(temp.out[temp.out[,"Name"]%in%snps.total & ((!is.na(temp.out[,"maf_ntotal"])) & as.numeric(as.character(temp.out[,"maf_ntotal"]))<max(mafRange) & as.numeric(as.character(temp.out[,"maf_ntotal"]))>min(mafRange)),"Name"])
      rsnps.maf <- temp.out[temp.out$Name%in%rsnps.in.gene,c("Name","maf_ntotal")]
      cmafUsed <- sum(rsnps.maf$maf_ntotal)
      nsnpsUsed <- length(rsnps.in.gene)
      nmiss <- sum(temp.out[temp.out$Name %in% rsnps.in.gene,"nmiss"])

      if (length(rsnps.in.gene)==0) gene.ssq.out[i,] <- list(gene,NA,NA,NA,NA,NA,NA,NA,NA) else {   
         
         if (length(rsnps.in.gene)==1) {                     
            ssq.out <- temp.out[temp.out$Name==rsnps.in.gene,]
            gene.ssq.out[i,] <- list(gene, ssq.out[9]^2, cmafTotal, cmafUsed, nsnpsTotal, nsnpsUsed, nmiss, 1,ssq.out[11])
         } else {
            z.1 <- temp.out[temp.out$Name%in%rsnps.in.gene,c("Name","Z.gc")]
            colnames(z.1)[2] <- "Z"
            z.1 <- merge(z.1,rsnps.maf,by="Name",sort=F)
            wi <- dbeta(z.1$maf,ssq.beta.wts[1],ssq.beta.wts[2])^2
            ssq <- sum(z.1$Z^2*wi)
            rsnps.in.gene.cor <- rsnpsingene.cor[[gene]][rsnps.in.gene,rsnps.in.gene]

            diag.temp<-diag(sqrt(wi))
            c <- eigen(diag.temp%*%rsnps.in.gene.cor%*%diag.temp)$values
            a <- sum(c^3)/sum(c^2)
            b <- sum(c)-((sum(c^2))^2)/sum(c^3)
            d <- ((sum(c^2))^3)/((sum(c^3))^2)
            p <- try(pchisq((ssq-b)/a,df=d,lower.tail=FALSE)) 

            if (!"try-error"%in%class(p)) gene.ssq.out[i,] <- list(gene,ssq,cmafTotal,cmafUsed,nsnpsTotal,nsnpsUsed,nmiss,d,p) else 
               gene.ssq.out[i,] <- list(gene,NA,NA,NA,NA,NA,NA,NA,NA)        
         } 
      }
   }

  SSQ.outfile <- paste(phen,"_SSQ_GC.txt",sep="")
 
  write.table(as.matrix(gene.ssq.out),SSQ.outfile,quote=F,row.names=F,sep="\t",na="NA")
  rdata.outfile <- paste(phen,"_GC.RData",sep="")   
  save(Rdata,file=rdata.outfile)
}      
