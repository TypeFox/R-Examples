`HierarchicalSparseCluster.wrapper` <-
  function(file,  method=c("average", "complete", "single", "centroid"), wbound=NULL, 
           silent=FALSE, cluster.features=FALSE, method.features=c("average", "complete", "single", "centroid"),output.cluster.files=TRUE,
           outputfile.prefix=NULL, maxnumgenes=5000, standardize.arrays=TRUE){
    method <- match.arg(method)
    method.features <- match.arg(method.features)
    dissimilarity <- "absolute.value"
    
    dat<-read.gct(file)
    x=t(dat$data)
    genedesc=dat$row.descriptions
    genenames=row.names(dat$data)
    dat <- NULL
    if(!is.null(maxnumgenes) && maxnumgenes<ncol(x)){
      vars <- apply(x,2,var, na.rm=TRUE)
      vars[apply(is.na(x),2,sum)>.5] <- 0
      whichers <- order(vars,decreasing=TRUE)[1:maxnumgenes]
      x <- x[,whichers]
      genenames <- genenames[whichers]
      genedesc <- genedesc[whichers]
    }
    perm.out <- NULL
    if(is.null(wbound)){
      perm.out<-HierarchicalSparseCluster.permute(x, nperms=10, wbounds=NULL, dissimilarity=dissimilarity, standardize.arrays=standardize.arrays)
      wbound=perm.out$bestw
    }

    if(is.null(outputfile.prefix)){outputfile.prefix=extract.prefix(file)}

    output <- HierarchicalSparseCluster(x=x,dists=perm.out$dists,method=method,wbound=wbound,dissimilarity=dissimilarity,
      silent=silent,cluster.features=cluster.features,method.features=method.features,output.cluster.files=output.cluster.files,
      outputfile.prefix=outputfile.prefix,genenames=genenames,genedesc=genedesc,standardize.arrays=standardize.arrays)
    if(!is.null(perm.out)) out <- list(perm.out=perm.out, out=output, x=x)
    if(is.null(perm.out)) out <- list(out=output,x=x)    
    if(!is.null(out$perm.out)){
      x <- out$perm.out
      cat("Results of running permutations to choose optimal wbound:", fill=TRUE, file=paste(outputfile.prefix,sep="",".out"))
      mat <- round(cbind(x$wbounds, x$nnonzerows, x$gaps, x$sdgaps),4)
      dimnames(mat) <- list(1:length(x$wbounds), c("Wbound", "# Non-Zero W's", "Gap Statistic", "Standard Deviation"))
      write.table(cbind(c("", rownames(mat)), rbind(colnames(mat),mat)),quote=FALSE,sep="\t",file=paste(outputfile.prefix,sep="",".out"), append=TRUE, col.names=FALSE, row.names=FALSE)
      cat("Tuning parameter that leads to largest Gap statistic: ", x$bestw, fill=TRUE,file=paste(outputfile.prefix,sep="",".out"), append=TRUE)
      cat(fill=TRUE, file=paste(outputfile.prefix,sep="",".out"), append=TRUE)
    }
    if(!is.null(out$perm.out)) cat("Results of running Sparse Clustering:", fill=TRUE, file=paste(outputfile.prefix,sep="",".out"), append=TRUE)
    if(is.null(out$perm.out)) cat("Results of running Sparse Clustering:", fill=TRUE, file=paste(outputfile.prefix,sep="",".out"), append=FALSE)
    x <- out$out
    cat("Wbound is ", x$wbound, ":", fill=TRUE, file=paste(outputfile.prefix,sep="",".out"), append=TRUE)
    cat("Number of non-zero weights: ", sum(x$ws!=0), fill=TRUE, file=paste(outputfile.prefix,sep="",".out"), append=TRUE)
    cat("Sum of weights: ", sum(x$ws), fill=TRUE, file=paste(outputfile.prefix,sep="",".out"), append=TRUE)
 
    # Print out gene names & weights
    ord <- order(out$out$ws,decreasing=TRUE)[1:sum(out$out$ws>0)]
    write.table(file=paste(outputfile.prefix,sep="",".txt"), rbind(c("Gene Name", "Gene Description", "Weight"), cbind(genenames[ord], genedesc[ord], out$out$ws[ord])), quote=FALSE, sep="\t\t", row.names=FALSE, col.names=FALSE)

    out$out$genenames <- genenames
    out$out$genedesc <- genedesc
    
    return(out)
  }

