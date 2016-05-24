
mixmapPI <-
  function(data.set,pval="pval",snp="snp",gene="gene",coord="coord",chr="chr",alpha=0.05){
    ############################
    #defining errors
    ############################
    #names must be specified
    #p-values git greg
    if (!pval%in%names(data.set)) stop(gettextf(paste0('Variable "',pval,'" not found in input data.frame.  Please specify variable name for p-values')))
    #snp
    if (!snp%in%names(data.set)) stop(gettextf(paste0('Variable "',snp,'" not found in input data.frame.  Please specify variable name for SNPs.')))
    #basepair
    if (!coord%in%names(data.set)) stop(gettextf(paste0('Variable "',coord,'" not found in input data.frame.  Please specify variable name for coordinate location.')))
    #chromosome
    if (!chr%in%names(data.set)) stop(gettextf(paste0('Variable "',chr,'" not found in input data.frame.  Please specify variable name for chromosome')))
    #gene
    if (!gene%in%names(data.set)) stop(gettextf(paste0('Variable "',gene,'" not found in input data.frame.  Please specify variable name for genes.')))
    
    #Lengths of input must be the same
    if (pval%in%names(data.set) & length(data.set[[pval]])!=length(data.set[[gene]])) stop(gettextf(paste('Lengths differ: Length of pval is ',length(data.set[[pval]]),'; Length of gene is ',length(data.set[[gene]]),sep="")))
    
    ############################
    #Warnings
    ############################
    #are pvalues numeric?
    if (!is.numeric(data.set[[pval]])) stop(gettextf('p-values must be numeric'))
    if (sum(is.na(data.set[[pval]]))>0) stop(gettextf('Some p-values are missing'))
    if (sum(is.na(data.set[[gene]]))>0) stop(gettextf('Some gene names are missing'))
    
    ############################
    #Pull out the subset of data that will be used
    ############################
    #Pull out the data that we need from the bigger data file
    datTemp<-data.frame(pvalTemp=data.set[[pval]],geneTemp=as.character(data.set[[gene]]),snpTemp=as.character(data.set[[snp]]))
    
    fret<-datTemp[datTemp$geneTemp!="",]
    
    #How many SNPs per gene?
    tab<-data.frame(table(datTemp$geneTemp))
    names(tab)<-c("gene","snpCount")
    
    #Inverse normal transformation of the p-values after ranking
    datTemp$probit.rank.transform<-qnorm((rank(datTemp$pvalTemp)-0.5)/length(datTemp$pvalTemp))
    
    #Run lmer function
    fm.rawg=lmer(probit.rank.transform ~ 1+(1|geneTemp),data=datTemp)
    aa=ranef(fm.rawg,condVar=TRUE)
    beta<-fixef(fm.rawg)
    post.est=aa$geneTemp[,1,]
    post.var=attr(aa[[1]],"postVar")[1,,]
    n.i<-as.vector(table(datTemp$geneTemp))
    sigma.sq.b<-VarCorr(fm.rawg)$geneTemp[1,1]
    sigma.sq<-(attr(VarCorr(fm.rawg),"sc"))^2
    #lambda<-sigma.sq.b/(sigma.sq.b+sigma.sq/n.i)
    
    var.out<-post.var
    
    
    ############################################################
    ##Calculate the prediction interval
    pred.upper<-post.est+sqrt(post.var)*qnorm(1-alpha/(length(post.est)))##Bonferroni correction
    ############################################################
    
    ############################
    #Defining Output
    ############################
    out<-data.frame(gene=as.character(rownames(aa$gene)),postEst=post.est,var=var.out,predUpper=as.numeric(as.character(pred.upper)))
    names(out)[1]<-"gene"
    out<-merge(out,tab,by.x="gene",by.y="gene",all.x=TRUE)
    
    data.set.g<-data.set[!duplicated(data.set[[gene]]),c(gene,chr,coord)]
    out<-merge(out,data.set.g,by.x="gene",by.y=gene,all.x=TRUE)
    out[[coord]]<-as.numeric(as.character(out[[coord]]))
    out[[chr]]<-as.numeric(as.character(out[[chr]]))
    
    #Calculate the gene level p-value
    out$MixMAP_pvalue<-pnorm(out$postEst/sqrt(out$var))
    out$MixMAP_pvalue_BonferroniAdjusted<-dim(out)[1]*pnorm(out$postEst/sqrt(out$var))
    out$MixMAP_pvalue_adj<-p.adjust(out$MixMAP_pvalue,method="BH")
    
    
    cutoff<-0
    num<-c("number detected"=sum(out$predUpper<cutoff),"total number of genes"=dim(out)[1])
    detected<-out[out$predUpper<cutoff,]
    
    ############################
    #If any genes detected
    ############################
    if (num[1]>0){
      genes.detect<-as.character(detected[["gene"]])
      snpTemp<-datTemp[datTemp$geneTemp%in%genes.detect,]
      snpTemp$geneTemp<-as.character(snpTemp$geneTemp)
      snpTemp$snpTemp<-as.character(snpTemp$snpTemp)
      
      #Pull out the min SNP 
      snp.tmp.list<-list()
      for (g in unique(snpTemp$geneTemp)){
        minSNP<-snpTemp[snpTemp$geneTemp==g,]
        snp.tmp.list[[g]]<-unlist(c(minSNP[min(minSNP$pvalTemp)==minSNP$pvalTemp,][1,-1],summary(minSNP$pvalTemp)))
      }
      
      snp.min<-data.frame(do.call(rbind,snp.tmp.list))
      names(snp.min)<-c("gene","minSNP","probit.rank.transform","pval.min","pval.Q1","pval.median","pval.mean","pval.Q3","pval.max")
      snp.min[["gene"]]<-as.character(snp.min[["gene"]])
      snp.min$minSNP<-as.character(snp.min$minSNP)
      snp.min$probit.rank.transform<-as.numeric(as.character(snp.min$probit.rank.transform))
      snp.min$pval.min<-as.numeric(as.character(snp.min$pval.min))
      snp.min$pval.Q1<-as.numeric(as.character(snp.min$pval.Q1))
      snp.min$pval.median<-as.numeric(as.character(snp.min$pval.median))
      snp.min$pval.mean<-as.numeric(as.character(snp.min$pval.mean))
      snp.min$pval.Q3<-as.numeric(as.character(snp.min$pval.Q3))
      snp.min$pval.max<-as.numeric(as.character(snp.min$pval.max))
      
      
      #merge on the min pvalue and SNP name
      detected<-merge(detected,snp.min,by.x="gene",by.y="gene",all.x=TRUE)
      
      #merge in the chromosome and Coordinate
      names(out)[c(1,6,7)]<-names(detected)[c(1,6,7)]<-c("gene","chr","coordinate")
      out[["gene"]]<-as.character(out[["gene"]])
      detected[["gene"]]<-as.character(detected[["gene"]])
      
      #detected.merg<-merge(detected,gene.location.file.num,by.x="gene",by.y="external_gene_id",all.x=TRUE)
      MixMAP.out<-new("MixMAP",output=out,num.genes.detected=num,detected.genes=detected,lmer.out=fm.rawg)
      
    }
    ############################
    #If no genes are detected
    ############################
    if (num[1]==0) {MixMAP.out<-new("MixMAP",output=out,num.genes.detected=num,detected.genes=detected,lmer.out=fm.rawg)}
    
    #return MixMAP object
    MixMAP.out
  }