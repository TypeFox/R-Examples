####
##
## ANOVA function
##
####

makeAllContrasts <- function(X, annot){
    ## change annot into data.frame if vector or factor
    if(is.null(dim(annot))) annot <- as.data.frame(annot, row.names=rownames(X), stringsAsFactor=FALSE)
    ## cast columns to character
    for (i in 1:ncol(annot)) annot[,i] <- as.character(annot[,i])
    ## find unique combinations of factors
    uniqAnnot <- annot[!duplicated(annot), ,drop=FALSE]
    ## remove the na's 
    uniqAnnot <- uniqAnnot[!apply(uniqAnnot, 1, function(x) any(is.na(x))), ,drop=FALSE]

    nbVar <- ncol(annot)

    nbLevelsByVar <- apply(uniqAnnot, 2, function(x) length(unique(x)))

    ## compute the number of tests
    nbTests <- 0
    for (i in 1:nbVar) nbTests <- nbTests + choose(nbLevelsByVar[i], 2) * prod(nbLevelsByVar[-i])
    
    ## create the contrasts matrix
    C <- matrix(0, nrow=nbTests, ncol=ncol(X))
    colnames(C) <- colnames(X)
    rownames(C) <- 1:nbTests
    
    ## fill in the contrast matrix
    idTest <- 1
    for (i in 1:(nrow(uniqAnnot)-1)){
        id1 <- rownames(uniqAnnot)[i]
        id1Name <- paste(colnames(uniqAnnot), uniqAnnot[i,], collapse=".", sep="") 
        
        for (j in (i+1):nrow(uniqAnnot)){
            if(sum(uniqAnnot[i,] != uniqAnnot[j,]) == 1) {
                id2 <- rownames(uniqAnnot)[j]
                              
                id2Name <- paste(colnames(uniqAnnot), uniqAnnot[j,], collapse=".", sep="")
                ids <- c(id1,id2)[order(c(id1Name,id2Name))]

                C[idTest,] <- X[ids[2], ] - X[ids[1],]
                
                rownames(C)[idTest] <- paste(sort(c(id1Name, id2Name), decreasing=TRUE), collapse="-")

                idTest <- idTest + 1
            }
        }
    }
    C <- C[order(rownames(C)), ,drop=FALSE]
    return(C)
}

test.nested.model <- function(X,X0,Y) {
  r <- ncol(X)
  r0 <- ncol(X0)
  n <- nrow(X)
  theta <- solve(t(X)%*%X)%*%t(X)%*%Y
  Y.pred <- X%*%theta
  Y0.pred <- X0%*%solve(t(X0)%*%X0)%*%t(X0)%*%Y

  ## Numerator
  NUM <- Y.pred-Y0.pred
  NUM <- NUM^2
  NUM <- colSums(NUM)/(r-r0)

  ## Denominator
  DEN <- Y-Y.pred
  DEN <- DEN^2
  DEN <- colSums(DEN)/(n-r)

  ## F
  F <- NUM/DEN

  ## Pvalue
  pvalue <- 1-pf(F,r-r0,n-r)

  return(list(theta=theta, F=F, pvalue=pvalue, residu=Y-Y.pred, sigma2=DEN, X=X))
}

inverse <- function(X){
  res.svd <- svd(X)
  U <- res.svd$u
  V <- res.svd$v
  p <- length(res.svd$d)
  D <- matrix(0,p,p)
  diag(D) <- 1/res.svd$d
  return(V%*%D%*%t(U))
}


test.LC <- function(C, X, Y, global=FALSE, cor.multtest=TRUE, typeFDR="FDR-BH") {
    if(!is.logical(global))
      stop("Error : logical value expected for option 'global' ")
    if(!is.logical(cor.multtest))
      stop("Error : logical value expected for option 'cor.multtest' ")
    r <- ncol(X)
    q <- nrow(C)   
    n <- nrow(X)
    if (class(C)!="matrix") C=t(as.matrix(C))

    invX <- solve(t(X) %*% X)
    theta <- invX %*% t(X) %*% Y
    Y.pred <- X %*% theta
    DEN <- (Y - Y.pred)^2
    sigma2 <- colSums(DEN)/(n - r)
    phi <- C %*% theta

    if(global){
      pvalue=vector(mode="numeric",length=length(sigma2))
      if(is.matrix(Y) & !is.null(colnames(Y))) names(pvalue)=colnames(Y)
      for(i in 1:length(sigma2)){
            var.theta <- invX*sigma2[i]
            Fvalue <- t(phi[,i])%*%inverse(C%*%var.theta%*%t(C))%*%phi[,i]/q
            pvalue[i] <- 1-pf(Fvalue,q,n-r)
      }
    }

    if(!global){
      Fvalue <- phi^2/(as.matrix(rowSums((C %*% invX) * C)) %*% as.matrix(t(sigma2)))
      pvalue <- 1 - pf(Fvalue, 1, n - r)
      # Here there are 2 possible configurations : 
      # 1. Y is a vector, p-value is a matrix with 1 column and p lines (one per contrast), multiple correction applies to the p different contrasts. This is the situation where we test many (more than 10) linear combinations of parameters of a linear model.
      # 2. Y is a matrix (gene expression typical setting), p-value is a matrix with n column (one per gene) and p lines (one per contrast), multiple correction applies to the n genes. This is the situation where we test a few combinations of parameters of a linear model on many (more than 30) variables or genes.
      # We treat the 2 situations the same way, that's why if we are in the first one we first transform the one column matrix in a one row matrix
      if(ncol(pvalue)==1) pvalue=t(pvalue)
      if(cor.multtest==TRUE){
            for(i in 1:nrow(pvalue)){
                  pvalue[i,]=multiple.correction(pvalue[i,],typeFDR)
            }
      }
    }

    return(list(Estimate = phi, Fvalue = Fvalue, pvalue = pvalue, Y.pred = Y.pred, resid = Y-Y.pred, sigma2 = sigma2, theta = theta))
}


####
##
## Two samples NON PARAMETRIC comparison
##
####

runWilcox <-  function(data,labels,typeFDR="FDR-BH",q=0.05,plot=TRUE){

  list.exp <- colnames(data)
  list.probe <-rownames(data)
  if (is.null(list.probe)){
    list.probe<-1:nrow(data)
  }
  message("Launch Wilcoxon test")

  test<-mt.teststat.num.denum(data, labels, test="wilcoxon")
  n<-length(which(labels==0))
  m<-length(which(labels==1))
  
  message("Calculate pval")
  espU<-(n*m)/2
  pval.test <- 2*pwilcox(espU-abs(test$teststat.num), n, m, lower.tail=TRUE)
  
  message("Adjusted pval")

  W<-espU-test$teststat.num
  padj<-multiple.correction(pval.test,typeFDR)
  indexsignif<-which(padj<=q)

  ## Result
  out <- data.frame(list.probe, W, pval.test,padj)
  colnames(out)<-c("probeID", "Stat", "RawpValue", "AdjpValue")

  ## Graphical plot
  if (plot){
    col <- rep("black", length(W))
    col[indexsignif] <- "red"
    ## Empirical distribution
    par(font.lab=2, frame=FALSE)
    qqplot(W,rwilcox(length(W), m=m, n=n), col=col[order(W)], main=paste("QQplot. n=",n,"m=",m), xlab="X")
  }

  return (out)
}


runSAM <- function(data, labels, nbpermut=500, q=0.05, plot=TRUE, method="d.stat", var.equal=TRUE, include.zero=FALSE, paired=FALSE, seed=123){
  
    ## pour mettre labels en 0 et 1 si non paired
    if(paired == FALSE){
        labels <- as.integer(factor(labels))-1
    }
    
    message("SAM analysis is running...")
    if(method=="d.stat")
        output <- sam(data, cl=labels, B=nbpermut, method="d.stat", var.equal=var.equal, include.zero=include.zero, rand=seed, control=samControl(delta=seq(0.1,10,0.05), lambda=0.5), med=TRUE)
    if(method=="wilc.stat")
        output <- sam(data, cl=labels, method="wilc.stat")
    if(method=="chisq.stat")
        output <- sam(data, cl=labels, method="chisq.stat", B=nbpermut, rand=seed)
    if(method=="trend.stat")
        output <- sam(data, cl=labels, method="trend.stat", B=nbpermut, rand=seed)
    
        
    message("Create table for all genes...")
    if(length(unique(labels))==2 & (method=="d.stat" || method=="wilc.stat")){
        alllist <- data.frame(rownames(data), output@d, output@p.value, output@fold)## output@q.value, output@fold)
        colnames(alllist)<-c("probeID", "Stat", "RawpValue", "FoldChange")## "AdjpValue", "FoldChange")
    }else{
        alllist <- data.frame(rownames(data), output@d, output@p.value)##, output@q.value)
        colnames(alllist) <- c("probeID", "Stat", "RawpValue")##, "AdjpValue")
    }
    alllist$probeID=as.character(alllist$probeID)

    print(output@mat.fdr[which(output@mat.fdr[,'FDR']>0),c("Delta","p0","False","Called","FDR")])

    message("Find delta...")
    outdelta <- try(findDelta(output, fdr=q))
    if(class(outdelta) == "matrix"){
        delta <- outdelta[2,1] 
        print(paste("Delta : ", delta, sep=""))
        
        message("Create SAM plot...")
        if(plot){
            #plot(output)
            plot(output, delta)
        }
        
        alllist$Significant <- rep(FALSE, nrow(alllist))  
        if(outdelta[2,1]!=0){
            delta.sum <- summary(output, delta)
            ds.mat.sig <- delta.sum@mat.sig
            rows <- ds.mat.sig$Row
            print(paste("Find",length(rows),"significant genes ..."))
            alllist[rows, "Significant"] <- TRUE
        }
    }
    else
        message("Can't find delta with chosen FDR... No SAM plot... No Significant slots in out data.frame...")

    message("The adjusted pvalue is a qvalue, and is not related to the significant genes found with the SAM's FDR")
    return(alllist)
}

foldchange <- function (data, labels, unlog=TRUE){
    if(missing(data))
        stop("** FAILURE : 'data' is missing **")
    if(missing(labels))
        stop("** FAILURE : 'labels' is missing **")

    data1.ave = apply(data[,which(labels==0)], 1, mean, na.rm=TRUE)
    data2.ave = apply(data[,which(labels==1)], 1, mean, na.rm=TRUE)
    fold.change = data2.ave - data1.ave
    
    if(unlog){
        fold.change <- 2^fold.change
    }
    
    return(fold.change)
}


####
##
## Two samples parametric comparison
##
####


runTtest <-  function(data,labels,typeFDR="FDR-BH",algo="t", q=0.05, plot=TRUE){

  list.exp <- colnames(data)
  list.probe <-rownames(data)
  if (is.null(list.probe)){
    list.probe<-1:nrow(data)
  }
  
  ## Student / Welch
  message("Launch ",algo," test")
  test<-mt.teststat(data,labels, test=algo)

  ## ddl
  if (algo == "t.equalvar"){
    s1<-apply(data[,which(labels==0)],1,function(x){return (length(which(!is.na(x))))})	
    s2<-apply(data[,which(labels==1)],1,function(x){return (length(which(!is.na(x))))})
    ddl <- s1+s2-2
  }
  else if(algo == "pairt"){
      s1<-apply(data[,which(labels==0)],1,function(x){return (length(which(!is.na(x))))})
      ddl <- s1-1  
  }
  else if (algo == "t"){
    s1<-apply(data[,which(labels==0)],1,function(x){return (length(which(!is.na(x))))})	
    s2<-apply(data[,which(labels==1)],1,function(x){return (length(which(!is.na(x))))})
    w1<-apply(data[,which(labels==0)],1,var,na.rm=TRUE)/s1
    w2<-apply(data[,which(labels==1)],1,var,na.rm=TRUE)/s2
    
    ##Satterthwaite approximation
    ddl<-(w1+w2)^2/((w1^2/(s1-1))+(w2^2/(s2-1)))
  }

  ## pval
  message("Calculate pval")
  pval.test <- 2*(pt(-abs(test),ddl))

  ## Multiple test correction
  message("Adjusted pval")
  padj<-multiple.correction(pval.test,typeFDR)
  indexsignif<-which(padj<=q)
  message(length(indexsignif), "significant genes")

  out <- data.frame(list.probe, test, pval.test,padj)
  colnames(out)<-c("probeID",  "Stat", "RawpValue", "AdjpValue")

  if (plot){
    ## Graphical plot
    col <- rep("black", length(test))
    col[indexsignif] <- "red"
    qqnorm(test, col=col)
    qqline(test)
  }
  return (out)
}


runIndTest<-function (data, labels, gene.names = NULL, plot = TRUE, dirname= NULL, grp.name=c("Group1","Group2")) {
    if (is.null(gene.names)) {
        if (is.null(rownames(data))) {
            gene.names = 1:nrow(data)
        }
        else {
            gene.names = rownames(data)
        }
    }
    if (is.vector(data)) {
        grp1 <- matrix(data[which(labels == 0)], nrow = 1)
        grp2 <- matrix(data[which(labels == 1)], nrow = 1)
        cpt <- 1
    }else {
        grp1 <- data[, which(labels == 0)]
        grp2 <- data[, which(labels == 1)]
        cpt <- nrow(data)
    }
    out <- as.data.frame(matrix(NA, ncol = 3, nrow = cpt))
    colnames(out) <- c("Probeset", "StatisticalTest", "P-values")
    if (cpt>25){
        plot=FALSE
        warning("Cannot plot more than 25 genes")
    }
 
    for (i in 1:cpt) {
        if(sd(round(grp1[i, ]),3)!= 0 && sd(round(grp2[i, ]),3)!= 0){
            testnorm.grp1 <- shapiro.test(grp1[i, ])
            testnorm.grp2 <- shapiro.test(grp2[i, ])
            
            ##loi non normale
            if (testnorm.grp1$p.value < 0.05 && testnorm.grp2$p.value < 0.05) {
                test <- "Wilcoxon"
                testmoyenne <- wilcox.test(grp1[i, ], grp2[i, ])
            }
            else {
                testvariance <- var.test(grp1[i, ], grp2[i, ])
                ##var non egale
                if (testvariance$p.value < 0.05) {
                    test <- "Welch"
                    testmoyenne <- t.test(grp1[i, ], grp2[i, ], var.equal = FALSE)
                }
                ##var egale
                else {
                    test <- "Student"
                    testmoyenne <- t.test(grp1[i, ], grp2[i, ], var.equal = TRUE)
                }
            }
            pval <- round(testmoyenne$p.value, 3)
            if (plot){
                if (!is.null(dirname)){
                    dirnamepng=paste(dirname,"/",rownames(data)[i],".png",sep="")
                    bitmap(file=dirnamepng,type="png16m",taa=4, gaa=4, height = 6, width = 6, res=150)
                }else{dev.new()}
                
                hg1 <- hist(grp1[i, ], plot = FALSE)
                hg2 <- hist(grp2[i, ], plot = FALSE)
                par(mfrow = c(1, 3))
                ym <- max(c(hg1$counts, hg2$counts), na.rm=TRUE)+1
                plot(hg1, col="light blue", main=grp.name[1], ylim=c(0,ym), xlab="Expression Level")
                points(density(grp1[i, ], na.rm = TRUE), type = "l", lwd = 1, col = "red")
                boxplot(grp1[i, ], grp2[i, ], main = paste(gene.names[i],"\n",test,"-",pval),names = grp.name, col = c("light blue", "#2896C8"))
                plot(hg2, main=grp.name[2], col="#2896C8", xlab="Expression Level", ylim=c(0,ym))
                points(density(grp2[i, ], na.rm = TRUE), type = "l",lwd = 1, col = "red")
                if (!is.null(dirname)){
                    dev.off()
                }
            }
            out[i, ] <- c(gene.names[i], test, pval)
       }else {out[i,]<-c(gene.names[i], "NA", "NA")}
    }
    return(out)
}


ordinal.chisq <- function(x){
    stopifnot(is.matrix(x))
    
    cat("\tValue\tdf\tpvalue\n")
    ##Pearson Chi-square
    ct <- chisq.test(x)
    cat("Pearson Chisq\t",ct$statistic,"\t",ct$parameter,"\t",ct$p.value,"\n")
    
    ##Linera-by-Linear association = A chi-square for the linear effect 
    ordc <- 0:(ncol(x)-1)
    ordr <- 0:(nrow(x)-1)
    
    c1 <- rep(ordr,apply(x,1,sum))
    c2 <- as.vector(unlist(apply(x,1,function(x){rep(ordc,x)})))
    lc <- cor.test(c1,c2)
    
    lla <- (sum(x)-1)*(lc$estimate)^2
    cat("Linera-by-Linear association\t",lla,"\t1\t",lc$p.value,"\n")

    cat("Deviation from linearity\t",ct$statistic-lla,"\t",ct$parameter-1,1-pchisq(ct$statistic-lla, df=ct$parameter-1),"\n")   
}


####
##
## Multiple testing functions
##
####
multiple.correction <-   function (pval,typeFDR,q){
  if (missing(pval))
    stop("Error : raw pvalues")
  if (missing(typeFDR))
    stop("Error : no typeFDR")
  if(!typeFDR %in% c("FWER","FDR-BH","FDR-TST","qvalue") )  
    stop("Error : unknown typeFDR")

  print (paste("typeFDR=", typeFDR))
  if (typeFDR == "FWER")
    padj <- FWER.Bonf(pval)

  else if (typeFDR == "FDR-BH")
    padj <- FDR.BH(pval)
  
  else if (typeFDR == "FDR-TST")
    padj <- FDR.BH(pval, TST=TRUE,q) 

  else if (typeFDR == "qvalue"){
    ##siggenes qvalue (version=2 for robust=TRUE)
    pi0 <- pi0.est(pval)$p0
    padj <- qvalue.cal(pval, pi0, version=2) 
  }

  return (padj)
}


FWER.Bonf <-  function (pval){
  m <- length(pval)
  signif <- rep(FALSE,m)
  
  ## Bonferroni
  resFDR<-mt.rawp2adjp(pval, proc=c("Bonferroni"))
  adjp <- resFDR$adjp[order(resFDR$index),]

  return (adjp[,"Bonferroni"])
}

FDR.BH<-  function (pval, TST=FALSE,q){
  m <- length(pval)
  signif <- rep(FALSE,m)

  ## LSU de BH & BY (1995)
  ## m0/m under-estimation of the FDR
  resFDR<-mt.rawp2adjp(pval, proc=c("BH"))
  adjp <- resFDR$adjp[order(resFDR$index),]
  
  if (TST){
      ## Second step (TST)
      ## Estimation of m0 to control FDR at q level
      r<-mt.reject(adjp[,"BH"],(q/(1+q)))$r
      outFDR<-adjp[,"BH"]*((m-r)/m)
  }
  else{
      outFDR<-adjp[,"BH"]
  }
  
  return (outFDR)
}


####
##
## GSA functions
##
####
GSA.correlate.txt <- function (GSA.genesets.obj, genenames) {
  
  nsets <- length(GSA.genesets.obj$genesets)
  ngenes <- unlist(lapply(GSA.genesets.obj$genesets, length))
  allgenes <- unlist(GSA.genesets.obj$genesets)
  sets.in.exp <- match(unique(allgenes), genenames)
  exp.in.sets <- match(genenames, allgenes)
  
  message("Number of gene-sets :", nsets)
  message("Total number of unique genes in gene-set collection :", length(unique(allgenes)))
  message("Total number of unique genes in genenames list :", length(unique(genenames)))
  message("Number of unique genes in both collections :", sum(!is.na(sets.in.exp)))
  nn <- rep(NA, nsets)
  for (i in 1:nsets) {
    nn[i] <- sum(!is.na(match(GSA.genesets.obj$genesets[[i]],genenames)))
  }
  message("Quantiles of fraction coverage of gene-sets")
  print(quantile(nn/ngenes, seq(0, 1, by = 0.1)))
  
  return()
}

runGSA <- function(nData, labels, gmtfile, chip="hgu133plus2" ,np=1000, minsize=10, maxsize=800, resp.type="Two class unpaired", fdr=0.25){

  #pkg <- chip
  pkg <- paste(chip, "db", sep=".")
  require(pkg, character.only=TRUE)

  if(resp.type != "Two class paired" && is.element(c(0,1), unique(labels))[1] == TRUE && is.element(c(0,1), unique(labels))[2] == TRUE)
    labels <- labels+1
  
  geneset.obj <- GSA.read.gmt(gmtfile)
  genenames <- mget(rownames(nData), eval(as.name(paste(chip,"SYMBOL", sep=""))))
  genenames <- unlist(genenames)
  minsize<-as.numeric(minsize)
  maxsize<-as.numeric(maxsize)

  print("GSA Analysis is running...")
  GSA.obj <- GSA(nData,labels, genenames=genenames, genesets=geneset.obj$genesets, resp.type=resp.type, method="maxmean", nperms=np, minsize=minsize ,maxsize=maxsize)

  message("GSA parameters :")
  message("Number of permutations :", np)
  message("Type of analysis :", resp.type)
  message("Minimun number of genes in a gene set :", minsize)
  message("Maximun number of genes in a gene set :", maxsize)
  message("FDR chosen :", fdr)
  GSA.correlate.txt(geneset.obj, genenames)

  results <- GSA.listsets(GSA.obj,geneset.names=geneset.obj$geneset.names, FDRcut=fdr, maxchar=150)

  ## genes scores for negative set
  gscore.negative=list()
  if(dim(results$negative)[[1]]!=0){
      for (i in 1:dim(results$negative)[[1]]){
          gscore.negative[[results$negative[,"Gene_set_name"][i]]]=GSA.genescores(as.numeric(results$negative[,"Gene_set"][i]),geneset.obj$genesets,GSA.obj,as.vector(genenames),negfirst=TRUE)
      }
  }
  ## genes scores for positive sets
  gscore.positive=list()
  if(dim(results$positive)[[1]]!=0){
      for (i in 1:dim(results$positive)[[1]]){
          gscore.positive[[results$positive[,"Gene_set_name"][i]]]=GSA.genescores(as.numeric(results$positive[,"Gene_set"][i]),geneset.obj$genesets,GSA.obj,as.vector(genenames))
      }
  }

  results$genescore.positive=gscore.positive
  results$genescore.negative=gscore.negative

  par(mfrow=c(nr=2, nc=2))
  hist(GSA.obj$pvalues.lo, breaks=50)
  hist(GSA.obj$pvalues.hi, breaks=50)
  hist(GSA.obj$fdr.lo, breaks=50)
  hist(GSA.obj$fdr.hi, breaks=50)

  message("GSA Analysis OK !")
  return(results)
}

