pamr.listgenes <- function (fit,   data, threshold, fitcv=NULL, genenames = FALSE)  {
  x <- data$x[fit$gene.subset, fit$sample.subset]
if (genenames) {
    gnames <- data$genenames[fit$gene.subset]
  }
  if (!genenames) {
    gnames <- NULL
  }
  geneid <- data$geneid[fit$gene.subset]
  if(!is.null(fit$y)){
       nc <- length(fit$y)
      }
 if(is.null(fit$y)){
       nc <- ncol(fit$proby)
      }
 clabs <- colnames(fit$centroids)

  aa <- pamr.predict(fit, x, threshold = threshold, type = "nonzero")
  cen <- pamr.predict(fit, x, threshold = threshold, type = "centroid")
  d <- (cen - fit$centroid.overall)[aa,, drop=FALSE]/fit$sd[aa]
  
  gene.order <- order(-apply(abs(d), 1, max))
  d <- round(d, 4)
  g <- gnames[aa]
  g1 <- geneid[aa]
  if (is.null(gnames)) {
    gnhdr <- NULL
  }
  if (!is.null(gnames)) {
    gnhdr <- "name"
  }

if(!is.null(fitcv)){
nfold=length(fitcv$cv.objects)

ind=matrix(F,nrow=nrow(x),ncol=nfold)
ranks=NULL
for( ii in 1:nfold){
	cen=pamr.predict(fitcv$cv.objects[[ii]], x[,-fitcv$folds[[ii]]],threshold=0, type="centroid")
	 dtemp <- (cen - fitcv$cv.objects[[ii]]$centroid.overall)[,, drop=FALSE]/fitcv$cv.objects[[ii]]$sd
	  r <- apply(abs(dtemp), 1, max)
	ranks=cbind(ranks,rank(-abs(r)))

	junk=pamr.predict(fitcv$cv.objects[[ii]], x[,-fitcv$folds[[ii]]],threshold=threshold, type="nonzero")
	ind[junk,ii]=T
}

av.rank=apply(ranks,1,mean)
av.rank=round(av.rank[aa],2)
prop=apply(ind[aa,,drop=F],1,sum)/nfold
}

  options(width = 500)
  schdr <- paste(clabs, "score", sep = "-")

if(is.null(fitcv)){
res <- cbind(as.character(g1), g, d)[gene.order,,drop=F]
  dimnames(res) <- list(NULL, c("id", gnhdr, schdr))

}
if(!is.null(fitcv)){
  res <- cbind(as.character(g1), g, d, av.rank, prop)[gene.order,,drop=F]
  dimnames(res) <- list(NULL, c("id", gnhdr, schdr, "av-rank-in-CV", "prop-selected-in-CV"))
}
  print(res, quote = FALSE)
}

