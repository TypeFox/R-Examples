cnvBatches<-function(intensities, batches, threshold.0, threshold.k, common.pi = TRUE, ...){
  batches <- as.character(batches)
  bb <- sort(unique(batches))
  nb <- length(bb)
  pb <- prop.table(table(batches))
  means <- sds <- pi <- NULL
  out <- rep(NA, length(intensities))
  mixtures <- list()
  for (i in 1:length(bb)){
    yy.i <- intensities[batches==bb[i]]
    cnv.i <- cnvDefault(yy.i, , , , , threshold.0, threshold.k, ...)
    means <- rbind(means, attr(cnv.i,"means"))
    sds <- rbind(sds, attr(cnv.i,"sds"))
    pi <- rbind(pi, attr(cnv.i,"pi"))
    mixtures[[i]] <- attr(cnv.i, "mixture")
  }
  k <- attr(cnv.i, "k")
  if (!missing(threshold.0) & !missing(threshold.k) & k<4)
      return(cnv(intensities,threshold.0 = threshold.0, threshold.k = threshold.k, ...))
  if (common.pi){
    pi <- pb%*%pi
    pi <- matrix(rep(pi,nb),nrow=nb,byrow=TRUE)
  }
  pp <- matrix(0, nrow = length(intensities), ncol = k)
  for (i in 1:length(bb))
    pp[batches==bb[i]] <- sapply(1:k, function(j) dnorm(intensities[batches==bb[i]], means[i, j], sds[i, j]) * pi[i, j])
  if (!missing(threshold.0))
    pp[intensities<threshold.0,]<-matrix(rep(rep(1:0,c(1,k-1)),sum(intensities<threshold.0)),ncol=k,byrow=TRUE)
  if (!missing(threshold.k))
    pp[intensities>threshold.k,]<-matrix(rep(rep(0:1,c(k-1,1)),sum(intensities>threshold.k)),ncol=k,byrow=TRUE)
  pp <- pp/rowSums(pp)
  if (length(bb) == 1)
    mixtures <- mixtures[[1]]
  out <- apply(pp, 1, which.max)
  num.copies <- attr(cnv.i, "num.copies")
  out <- num.copies[out]
  attr(out, "probabilities") <- pp
  attr(out, "means") <- means
  attr(out, "sds") <- sds
  attr(out, "pi") <- pi
  rownames(attr(out, "pi")) <- rownames(attr(out, "means")) <- rownames(attr(out, "sds")) <- bb
  attr(out, "k") <- k
  attr(out, "meanRatio") <- intensities
  attr(out, "num.copies") <- num.copies
  attr(out, "batches") <- batches
  attr(out, "mixture") <- mixtures
  k <- attr(cnv.i, "k")
  class(out) <- "cnv"
  out
}
