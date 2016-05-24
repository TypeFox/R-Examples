"defectiveCVdisc" <-
function(x, cl, nfold=NULL, FUN=aovFbyrow,
           nfeatures=2, seed=31, funda=lda, foldids=NULL,
           subset=NULL, print.progress=TRUE){
    ## Option to omit one or more points
    if(!is.null(subset)) cl[!is.na(cl)][!subset] <- NA
    if(length(nfold)==1)nfold <- c(nfold,1)
    if(any(is.na(cl))){x <- x[,!is.na(cl)]
                       cl <- cl[!is.na(cl)]
                     }
    nobs <- dim(x)[2]
    ## Get fold information from foldids, if specified,
    ## else if nfold is not specified, use leave-one-out CV
    if(!is.null(foldids))
      nfold <- c(length(unique(foldids)), dim(foldids)[2])
    if(is.null(nfold)&is.null(foldids))nfold <- sum(!is.na(cl))
    else if(nfold[1]==nobs)foldids <- sample(1:nfold[1])
    else foldids <- sapply(1:nfold[2], function(x)
                     divideUp(cl, nset=nfold[1]))
    if(length(nfold)==1)nfold <- c(nfold,1)
    cl <- factor(cl)
    ngp <- length(levels(cl))
    genes <- rownames(x)
     if(is.null(genes)){
      genes <- paste(1:dim(x)[1])
      print("Input rows (features) are not named. Names")
      print(paste(1,":", dim(x)[1], " will be assigned.", sep=""))
      rownames(x) <- genes
    }
    if(!is.null(seed))set.seed(seed)
    Fcut <- NULL
    maxgenes <- max(nfeatures)

    stat <- FUN(x=x, cl)
    Fcut <- list(F=sort(stat, decreasing=TRUE)[nfeatures],
                 df=c(ngp-1, nobs-ngp))
    ord <- order(-abs(stat))[1:maxgenes]
    genes.ord <- genes[ord]
    selectonce.df <- data.frame(t(x[ord, , drop=FALSE]))
    acc.resub <- acc.sel1 <- numeric(maxgenes)
    if(nfold[1]==0)acc.sel1 <- NULL

    for(ng in nfeatures){
      resub.xda <- funda(cl~., data=selectonce.df[,1:ng,drop=FALSE])
      hat.rsb <- predict(resub.xda)$class
      tab.rsb <- table(hat.rsb, cl)
      acc.resub[ng] <- sum(tab.rsb[row(tab.rsb)==col(tab.rsb)])/sum(tab.rsb)
      if(nfold[1]==0)next
      if(nfold[1]==nobs){
        hat.sel1 <- funda(cl~., data=selectonce.df[,1:ng,drop=FALSE],
                          CV=TRUE)$class
        tab.one <- table(hat.sel1, cl)
        acc.sel1[ng] <- sum(tab.one[row(tab.one)==col(tab.one)])/sum(tab.one)
      } else
      {
      hat <- cl
      if(print.progress)cat(paste(ng,":",sep=""))
      for(k in 1:nfold[2])
      {
        foldk <- foldids[,k]
        ufold <- sort(unique(foldk))
        for(i in ufold){
          testset <- (1:nobs)[foldk==i]
          trainset <- (1:nobs)[foldk!=i]
          dfi <- selectonce.df[-testset, 1:ng, drop=FALSE]
          newdfi <- selectonce.df[testset, 1:ng, drop=FALSE]
          cli <- cl[-testset]
          xy.xda <- funda(cli~., data=dfi)
          subs <- match(colnames(dfi), rownames(df))
          newpred.xda <- predict(xy.xda, newdata=newdfi, method="debiased")
          hat[testset] <- newpred.xda$class
        }
        tabk <- table(hat,cl)
        if(k==1)tab <- tabk else tab <- tab+tabk
      }
      acc.sel1[ng] <- sum(tab[row(tab)==col(tab)])/sum(tab)
      }
    }
    if(print.progress)cat("\n")
    invisible(list(acc.resub=acc.resub, acc.sel1=acc.sel1, genes=genes.ord))
  }

