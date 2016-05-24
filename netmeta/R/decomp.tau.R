decomp.tau <- function(x, tau.preset=0){
  
  
  nmak <- nma.krahn(x, tau.preset=tau.preset)
  ##
  design <- nmak$design
  studies <- nmak$studies
  V <- nmak$V
  V.studies <- nmak$V.studies
  X.obs <- nmak$X.obs
  multicomp <- nmak$multicomp
  
  
  multiarm <- any(studies$narms>2)
  
  
  designs <- unique(as.character(design$design))
  TE.net.minus <- matrix(NA, nrow=nrow(X.obs), ncol=length(designs))
  colnames(TE.net.minus) <- designs
  ##
  for (i in designs){
    Xb <- X.obs[!(rownames(X.obs)==i),]
    Vb <- V[!(rownames(V)==i), !(colnames(V)==i)]
    if (qr(Xb)$rank==ncol(X.obs)){
      TE.net.minus[,i] <- X.obs %*% solve(t(Xb) %*% solve(Vb) %*% Xb) %*%
        t(Xb) %*% solve(Vb) %*% design$TE.dir[!(design$design==i)]
    }
  }
  ##
  rownames(TE.net.minus) <- design$design
  
  
  freq.d  <- rep(NA, nrow(design))
  Q.het.design <- rep(NA, nrow(design))
  ##
  for (i in unique(sort(as.character(studies$design)))){
    Vs <- V.studies[rownames(V.studies)==i, colnames(V.studies)==i]
    ##
    TE.i <- studies$TE[as.character(studies$design)==i]
    TE.dir.i <- studies$TE.dir[as.character(studies$design)==i]
    ##
    Q.het.design[which(design$design==i)[1]] <- t(TE.i-TE.dir.i) %*% solve(Vs) %*% (TE.i-TE.dir.i)
    freq.d[which(design$design==i)[1]] <- length(TE.i)
  }
  ##
  NAfd <- is.na(freq.d)
  design$narms[NAfd] <- NA
  design$freq[NAfd] <- NA
  
  
  df.het.design <- freq.d - (design$narms-1)
  pval.het.design <- 1-pchisq(Q.het.design, df=df.het.design)
  pval.het.design[df.het.design==0] <- NA
  ##
  Q.het <- t(studies$TE-studies$TE.dir) %*% solve(V.studies) %*% (studies$TE-studies$TE.dir)
  ## Formula (7) in Krahn et al. (2013)
  Q.net <- t(studies$TE-studies$TE.net) %*% solve(V.studies) %*% (studies$TE-studies$TE.net)
  ## Formula (8) in Krahn et al. (2013)
  Q.inc <- Q.net - Q.het
  ##
  df.het <- sum(df.het.design, na.rm=TRUE)
  df.net <- sum(freq.d, na.rm=TRUE) - ncol(X.obs)
  df.inc <- df.net - df.het
  ##
  Qid <- t(design$TE.dir-design$TE.net) %*% solve(V) * (design$TE.dir-design$TE.net)
  nam <- colnames(Qid)
  Q.inc.design <- as.vector(Qid)
  names(Q.inc.design) <- nam
  
  
  residuals <- apply(TE.net.minus, 2, function(x) design$TE.dir-x)
  residuals <- residuals[, rownames(residuals)]
  diag(residuals) <- rep(0, nrow(residuals))
  ##
  if (multiarm)
    for (i in 1:length(multicomp))
      residuals[rownames(residuals)==multicomp[i], colnames(residuals)==multicomp[i]] <- 0
  
  
  Q.inc.detach <- apply(residuals, 2, function(x) t(x) %*% solve(V) %*% x)
  Q.inc.detach[NAfd] <- NA
  ##
  df.inc.detach <- nrow(X.obs) - (ncol(X.obs) + (design$narms-1))
  df.inc.detach[df.inc.detach<0] <- NA
  df.inc.detach[NAfd] <- NA
  ##
  pval.inc.detach <- 1-pchisq(Q.inc.detach, df=df.inc.detach)
  pval.inc.detach[df.inc.detach==0] <- NA
  
  
  Q.decomp <- data.frame(Q=c(Q.net, Q.het, Q.inc),
                         df=c(df.net, df.het, df.inc),
                         pval=1-c(pchisq(Q.net, df=df.net),
                           pchisq(Q.het, df=df.het),
                           pchisq(Q.inc, df=df.inc)))
  ##
  Q.decomp$pval[c(df.net, df.net, df.inc)==0] <- NA
  ##
  rownames(Q.decomp) <- c("Whole network",
                          "Within designs",
                          "Between designs")
  
  
  Q.het.design <- data.frame(design=design$design,
                             Q=round(Q.het.design, 15),
                             df=df.het.design,
                             pval=pval.het.design)
  ##
  Q.inc.detach <- data.frame(design=design$design,
                             Q=Q.inc.detach,
                             df=df.inc.detach,
                             pval=pval.inc.detach)
  ##
  if (multiarm){
    Q.het.design <- Q.het.design[!(duplicated(Q.het.design$design)),]
    Q.inc.detach <- Q.inc.detach[!(duplicated(Q.inc.detach$design)),]
  }
  ##
  if (any(is.na(Q.inc.detach$Q)))
    Q.inc.detach <- Q.inc.detach[!(is.na(Q.inc.detach$Q)),]
  
  res <- list(
    Q.decomp=Q.decomp,
    Q.het.design=Q.het.design,
    Q.inc.detach=Q.inc.detach,
    Q.inc.design=Q.inc.design,
    residuals.inc.detach=residuals
    )
  
  res
}
