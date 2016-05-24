PCgamma <-
function(formula,data,snpprefix="snp",gene,
                  PCpctVar=80,gammaShape=1,
                  STT=NULL,pheno.type=c("case.control","quantitative","survival"),
                  perm=T,n.perm=1000,seed=12212012) {  

    set.seed(seed)
    snplist <- names(data)[substring(names(data),1,nchar(snpprefix))==snpprefix]

    if(is.null(STT) & is.null(gammaShape)) stop("need either STT or gammaShape specified")
    if(!is.null(STT) & !is.null(gammaShape)) stop("specify only one of STT or gammaShape")
        
    if(!is.null(STT)) gammaShape<-STTtoShapeParameter(STT)
    phenotype <- as.character(formula)[2]
    if(grepl("Surv",substring(phenotype,1,4))==TRUE) {
        strings <- gsub(" ","",unlist(strsplit(phenotype,"\\(|\\)|\\,")))
        if(length(strings)==3) pheno <- Surv(data[,strings[2]],data[,strings[3]])
        if(length(strings)==4) pheno <- Surv(data[,strings[2]],data[,strings[3]],data[,strings[4]])
        if(class(pheno)!="Surv") stop("error: need survival object")
        } else pheno <- data[,phenotype]
    if(class(pheno) != "list") pheno<-list(pheno)
    if(!perm) n.perm<-0

    i<-2
    while(i<=n.perm+1) {
      if(pheno.type=="survival") {
        temp <- sample(c(1:nrow(pheno[[1]])),replace=F) 
        pheno[[i]] <- pheno[[1]][c(temp)]
      } else {
        pheno[[i]]<-sample(pheno[[1]],replace=F)
      }
        i<-i+1
    }
    rm(i)

    u.gene<-unique(gene)
    geno <- data[,snplist]
    pc.x<-lapply(u.gene,function(p1,p2,p3,p4) {
                        X<-p1[,p2==p3,drop=F]
                        X<-X[,!duplicated(t(X)),drop=F]
                        pc<-prcomp(X)
                        topXpct<-1
                        if(pc$sdev[1]!=0 && (pc$sdev[1]^2/sum(pc$sdev^2))<p4) {topXpct<-1:(max(which((cumsum(pc$sdev^2)/sum(pc$sdev^2))<p4))+1)}
                        pc<-list(pcs=pc$x[,topXpct,drop=F],npcs=max(topXpct),nsnps=dim(X)[2])
                                                return(pc)
                           },p1=geno,p2=gene,p4=PCpctVar/100)

        gene.info<-cbind(u.gene,do.call(rbind,lapply(pc.x,function(x) c(x$nsnps,x$npcs))))
        dimnames(gene.info)[[2]]<-c("geneName","nSNPs","nPCs")
        pc.x<-lapply(pc.x,function(x) x$pcs)
          
    #if(is.character(covariates)) {
    #    covariate.formula<-paste("+",covariates)
    #}else{
    #    covariate.formula<-ifelse(is.null(covariates),"",paste("+",deparse(substitute(covariates))))
    #}   
    covariate.formula <- gsub("\\s","",as.character(formula)[3])
    #attach(data)

    ### run analyses of phenotype vs. PCs based on genotypes within genes ###
    
    if(pheno.type=="case.control") {
        pvalues<-lapply(pheno,function(a1,a2,a3) {
                                unlist(lapply(a2,function(p1,p2,p3) {
                                        full<-glm(update(p1~p2,paste("~ .+",p3,sep="")),family="binomial",data=data)
                                        reduced<-glm(update(p1~1,paste("~ .+",p3,sep="")),family="binomial",data=data,subset=rowSums(is.na(p2))==0)
                                        lr.p<-pchisq(reduced$deviance - full$deviance,df=reduced$df.residual-full$df.residual,lower.tail=F)
                                        },p1=a1,p3=a3))
                             },a2=pc.x,a3=covariate.formula)

    }else{

        if(pheno.type=="quantitative") {
            pvalues<-lapply(pheno,function(a1,a2,a3) {
                                    unlist(lapply(a2,function(p1,p2,p3) {
                                                full<-lm(update(p1~p2,paste("~ .+",p3)),data=data)
                                                reduced<-lm(update(p1~1,paste("~ .+",p3)),data=data,subset=rowSums(is.na(p2))==0)
                                                out<-anova(full,reduced,test="F")[2,"Pr(>F)"]
                                                return(out)
                                            },p1=a1,p3=a3))
                                  },a2=pc.x,a3=covariate.formula)

        }else{
            if(class(pheno[[1]])!="Surv") stop("pheno indicated as survival but not a Surv object")
            pvalues<-lapply(pheno,function(a1,a2,a3) {
                                    unlist(lapply(a2,function(p1,p2,p3) {
                                                    full<-coxph(update(p1~p2,paste("~ .+",p3,sep="")),data=data)
                                                    reduced<-coxph(update(p1~1,paste("~ .+",p3,sep="")),data=data,subset=rowSums(is.na(p2))==0)
                                                    #if reduced model is null use regular likelihood ratio test??#
                                                    if(is.na(reduced$loglik[2])) {
                                                      blah <- summary(full)
                                                      lr.p <- blah$logtest[3]
                                                    }
                                                    else {
                                                    lr.p<-pchisq(-2*(reduced$loglik[2] - full$loglik[2]),df=length(full$coefficients)-length(reduced$coefficients),lower.tail=F)
                                                  }
                                                    },p1=a1,p3=a3))
                                  },a2=pc.x,a3=covariate.formula)


        }
    }

    pvalues<-do.call(rbind,pvalues)
    gene.info<-cbind(gene.info,pvalues[1,])
    dimnames(gene.info)[[2]][4]<-"geneP"
       
    gamma.stat<-rowSums(qgamma(1-signif(pvalues),shape=gammaShape),na.rm=T)
    if(perm) {
      gamma.perm.pvalue<-mean(gamma.stat[-1]>=gamma.stat[1])
      gene.perm.pvalue<-colMeans( signif(pvalues[-1,]) < matrix(signif(pvalues[1,]),n.perm,nrow(gene.info),byrow=TRUE))
      gene.info<-cbind(gene.info,gene.perm.pvalue)
      dimnames(gene.info)[[2]][5]<-"genePperm"      
    }else{
            gamma.perm.pvalue<-NULL
    }

    gamma.obs.pvalue<-pgamma(gamma.stat[1],shape=gammaShape*dim(pvalues)[2],lower.tail=F)
    return(list(gamma.pvalue=gamma.obs.pvalue,perm.pvalue=gamma.perm.pvalue,gene.info=gene.info))
}
