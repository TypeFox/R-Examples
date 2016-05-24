nma.krahn <- function(x, tau.preset=0){


  if (!inherits(x, "netmeta"))
    stop("Argument 'x' must be an object of class \"netmeta\"")


  n <- x$n


  if (x$reference.group=="")
    trts <- colnames(x$A.matrix)
  else
    trts <- c(x$reference.group,
              colnames(x$A.matrix)[colnames(x$A.matrix)!=x$reference.group])


  studies.pre <- data.frame(studlab=x$studlab,
                            treat1=x$treat1, treat2=x$treat2,
                            TE=-x$TE, seTE=sqrt(x$seTE^2+tau.preset^2),
                            narms=x$narms[match(x$studlab, x$studies)],
                            stringsAsFactors=FALSE)
  ##
  studies <- studies.pre <- studies.pre[order(studies.pre$studlab),]


  twoarm   <- any(studies$narms==2)
  multiarm <- any(studies$narms>2)
  selmulti <- studies$narms>2


  sel <- studies.pre$treat2==x$reference.group
  ##
  studies$treat1[sel] <- studies.pre$treat2[sel]
  studies$treat2[sel] <- studies.pre$treat1[sel]
  studies$TE[sel] <- -studies.pre$TE[sel]
  studies <- data.frame(studies,
                        comparison=paste(studies$treat1, studies$treat2, sep=":"))


  comparison.num.poss <- n*(n-1)/2
  comparisons <- levels(factor(as.character(studies$comparison)))
  comparison.num <- length(comparisons)


  trts.poss <- rep(NA, comparison.num.poss)
  k <- 1
  for (i in 1:(n-1))
    for (j in (i+1):n){
      trts.poss[k] <- paste(trts[i], trts[j], sep=":")
      k <- k+1
    }


  direct <- matrix(NA, nrow=comparison.num, ncol=6)
  ##
  colnames(direct) <- c("comparison", "TE", "seTE",
                        "TE.2arm", "seTE.2arm", "n.2arm")
  ##
  direct <- data.frame(direct)
  j <- 0
  ##
  for (i in names(table(studies$comparison))){
    j <- j+1
    ##
    TE.i <- studies$TE[studies$comparison==i]
    seTE.i <- studies$seTE[studies$comparison==i]
    m1 <- metagen(TE.i, seTE.i, sm=x$sm)
    ##
    direct$comparison[j] <- i
    direct$TE[j] <- m1$TE.fixed
    direct$seTE[j] <- m1$seTE.fixed
    ##
    if (sum(studies$comparison==i & !selmulti) > 0) {
      TE.i <- studies$TE[studies$comparison==i&studies$narms==2]
      seTE.i <- studies$seTE[studies$comparison==i&studies$narms==2]
      m2 <- metagen(TE.i, seTE.i, sm=x$sm)
      ##
      direct$TE.2arm[j] <- m2$TE.fixed
      direct$seTE.2arm[j] <- m2$seTE.fixed
      direct$n.2arm[j] <- m2$k
    }
  }


  if (multiarm){
    multistudies <- split(studies[selmulti,], as.character(studies$studlab[selmulti]))
    multistudies <- lapply(multistudies,
                           function(x) x[which(x$treat1==names(which.max(table(x$treat1)))),])
    multistudies <- lapply(multistudies,
                           function(x) x[order(x$treat2),])
    ##
    des <- lapply(multistudies, function(x) paste(c(x$treat1[1], x$treat2), collapse=":"))
    multistudies <- data.frame(unsplit(multistudies,
                                       rep(names(multistudies),
                                           unlist(lapply(multistudies, function(x) nrow(x))))),
                               design=unsplit(des, rep(names(multistudies),
                                 unlist(lapply(multistudies, function(x) nrow(x))))))
    ##
    multistudies <- data.frame(multistudies,
                               des=paste(multistudies$comparison, multistudies$design, sep="_"))
    ##
    row.names(studies) <- NULL
    multistudies2 <- split(studies[selmulti,], as.character(studies$studlab[selmulti]))
    multistudies2 <- lapply(multistudies2, function(x) x[do.call(order,x[,c("treat1","treat2")]),])
    multistudies2 <- lapply(multistudies2,
                            function(x)
                            rbind(
                              x[x$treat1==names(which.max(table(x$treat1))),],
                              x[x$treat1 %in% names(table(x$treat1)[-which.max(table(x$treat1))]),]))
    multistudies2 <- unsplit(multistudies2,
                             rep(names(multistudies2),
                                 unlist(lapply(multistudies2, function(x) nrow(x)))))
  }

  studies <- data.frame(studies, design=studies$comparison)
  if (multiarm & sum(is.na(direct$seTE.2arm))>0)
    direct2 <- data.frame(direct[!is.na(direct$seTE.2arm),])
  else
    direct2 <- direct
  ##
  direct2 <- data.frame(direct2)
  
  
  V.design <- diag(direct2$seTE.2arm^2,
                   nrow=length(direct2$seTE.2arm),
                   ncol=length(direct2$seTE.2arm))
  
  
  if (multiarm){
    sp <- split(multistudies2, multistudies2$studlab)
    armM <- unlist(lapply(split(multistudies2$narms, multistudies2$studlab), function(x) x[1]))
    ##
    covs <- lapply(sp,
                   function(x){
                     n <- x$narms[1]
                     k <- 0
                     m <- matrix(NA, nrow=n-1, ncol=n-1)
                     for (i in 1:(n-2)){
                       for (j in (i+1):(n-1)){
                         m[i, j] <- (x$seTE[i]^2 + x$seTE[j]^2 - x$seTE[n+k]^2)/2
                         m[j, i] <- (x$seTE[i]^2 + x$seTE[j]^2 - x$seTE[n+k]^2)/2
                         k <- k+1
                       }
                     }
                     diag(m) <- x$seTE[1:(n-1)]^2
                     m
                   })
    ##
    V3 <- NA
    ##
    for (i in 1:length(covs))
      V3 <- magic::adiag(V3, covs[[i]])
    ##
    V3 <- V3[-1, -1]
    ##
    V.studies <- magic::adiag(diag(studies$seTE[!selmulti]^2), V3)
    colnames(V.studies) <- c(as.character(studies$design[!selmulti]),
                             as.character(multistudies$design))
    rownames(V.studies) <- c(as.character(studies$design[!selmulti]),
                             as.character(multistudies$design))
    ##
    multicomp <- names(which(table(multistudies$design)>0))
    V3.agg <- NA
    TE.agg <- NA
    ##
    for (i in 1:length(multicomp)){
      studlabM <- unique(multistudies$studlab[multistudies$design==multicomp[i]])
      ncovs <- covs[names(covs) %in% studlabM]
      l <- sapply(ncovs, solve)
      dim <- multistudies$narms[multistudies$studlab==studlabM[1]][1]-1
      covs3 <- solve(matrix(apply(l, 1, sum), nrow=dim))
      V3.agg <- magic::adiag(V3.agg, covs3)
      m <- matrix(NA, nrow=dim, ncol=length(studlabM))
      for (j in 1:length(studlabM))
        m[,j] <- matrix(l[, j], nrow=dim) %*%
          multistudies$TE[multistudies$studlab==studlabM[j]]
      ##
      TE.agg <- c(TE.agg,covs3 %*% apply(m, 1, sum))
    }

    V3.agg <- V3.agg[-1,-1]
    TE.agg <- TE.agg[-1]

    V <- magic::adiag(V.design, V3.agg)
    ##
    nam <- rep(multicomp, unlist(lapply(split(multistudies, multistudies$design),
                                        function(x) x$narms[1]))-1)
    ##
    if (any(twoarm))
      rownames(V) <- colnames(V) <- c(direct2$comparison, nam)
    else
      rownames(V) <- colnames(V) <- nam
    ##
    TE.dir <- c(direct2$TE.2arm, TE.agg)
  }
  else{
    V <- magic::adiag(V.design)
    rownames(V) <- direct2$comparison
    colnames(V) <- direct2$comparison
    TE.dir <- direct2$TE.2arm
    V.studies <- diag(studies$seTE[!selmulti]^2)
    colnames(V.studies) <- rownames(V.studies) <- as.character(studies$comparison[!selmulti])
  }
  ##
  if (min(eigen(V, only.values = TRUE)$values)<0)
    stop("Covariance matrix is not nnd")


  fX <- function(n){
    possK <- n*(n-1)/2
    X <- matrix(0, nrow=possK, ncol=n-1)
    X[1:(n-1), 1:(n-1)] <- diag(rep(-1, n-1))
    X[n*(n-1)/2, (n-2):(n-1)] <- cbind(1, -1)
    if (n*(n-1)/2-(n-1)>1){
      l <- n
      j <- n-2
      u <- n+j-1
      for (k in 1:(n-3)){
        X[l:u, k:(n-1)] <- cbind(1, diag(rep(-1, n-k-1)))
        j <- j-1
        l <- u+1
        u <- l+j-1
      }
    }
    X
  }
  ##
  X.full <- fX(n)
  rownames(X.full) <- trts.poss
  colnames(X.full) <- trts.poss[1:n-1]
  ##
  X.obs2.design <- X.full[direct2$comparison,, drop=FALSE]


  if (multiarm){
    num.basics.design <- unlist(lapply(split(multistudies,multistudies$design),
                                       function(x) x$narms[1]))-1
    ##
    basics <- lapply(split(multistudies, multistudies$design),
                     function(x) split(x, x$studlab)[[1]]$comparison)
    basics <- unsplit(basics, rep(1:length(multicomp), num.basics.design))
    ##
    X.obs3.design <- X.full[as.character(basics),]
    rownames(X.obs3.design) <- rep(multicomp, num.basics.design)
    X.obs <- rbind(X.obs2.design, X.obs3.design)
  }
  else
    X.obs <- X.obs2.design
  ##
  H <- X.full %*% solve(t(X.obs) %*% solve(V) %*% X.obs) %*% t(X.obs) %*% solve(V)
  TE.net <- H %*% TE.dir


  covTE.net.base <- solve(t(X.obs) %*% solve(V) %*% X.obs)
  co <- NA
  for (i in 1:(n-2)){
    for (j in 2:(n-1)){
      if (i != j && i<j){
        co <- c(co,
                diag(covTE.net.base)[i] +
                diag(covTE.net.base)[j] -
                2*covTE.net.base[i,j])
      }
    }
  }
  ##
  covTE.net <- c(diag(covTE.net.base), co[-1])


  comps <- as.character(studies$comparison[!selmulti])
  studlabs <- as.character(studies$studlab[!selmulti])
  ##
  if (multiarm){
    comps <- c(comps, as.character(multistudies$comparison))
    studlabs <- c(studlabs, as.character(multistudies$studlab))
  }


  X.obs.studies <- X.full[comps,]


  H.studies <- X.full %*%
    solve(t(X.obs.studies) %*% solve(V.studies) %*% X.obs.studies) %*%
      t(X.obs.studies) %*% solve(V.studies)
  ##
  colnames(H.studies) <- studlabs


  network <- data.frame(TE=TE.net, seTE=sqrt(covTE.net))


  if (multiarm){
    len.designs <- c(rep(1, length(direct2$comparison)),
                     unlist(lapply(strsplit(multicomp, ":"),
                                   function(x){length(x)-1})))
    freq <- rep(c(direct2$n.2arm,
                  unlist(lapply(split(multistudies, multistudies$design),
                                function(x) length(names(table(x$studlab)))))),
                len.designs)
    narms <- rep(c(rep(2, nrow(direct2)),
                   unlist(lapply(strsplit(multicomp, ":"), function(x) length(x)))),
                 len.designs)
    ##
    design <- data.frame(design=c(direct2$comparison,
                           rep(multicomp, unlist(lapply(strsplit(multicomp,":"),
                                                        function(x) length(x)))-1)),
                         comparison=c(direct2$comparison,
                           as.character(unlist(lapply(split(multistudies,
                                                            multistudies$design),
                                                      function(x) as.character(unlist(
                                                        split(x,x$studlab)[[1]]["comparison"])))))),
                         narms=narms,
                         freq=freq,
                         TE.dir=TE.dir,
                         seTE.dir=sqrt(diag(V)))
  }
  else{
    len.designs <- c(rep(1, length(direct2$comparison)))
    freq <- rep(c(direct2$n.2arm), len.designs)
    narms <- rep(c(rep(2, nrow(direct2))), len.designs)
    ##
    design <- data.frame(design=colnames(V),
                         comparison=colnames(V),
                         narms=rep(2, length(direct2$comparison)),
                         freq=direct2$n.2arm,
                         TE.dir=TE.dir,
                         seTE.dir=sqrt(diag(V)))
  }
  ##
  rownames(design) <- NULL
  ##
  design <- data.frame(design,
                       TE.net=network[as.character(design$comparison), "TE"],
                       seTE.net=network[as.character(design$comparison), "seTE"])


  if (multiarm)
    studies <- rbind(studies[!selmulti,], multistudies[, 1:8])
  ##
  studies <- studies[, c("studlab", "treat1", "treat2",
                         "TE", "seTE", "narms",
                         "design", "comparison")]
  ##
  studies <- merge(studies, design[, names(design) != "narms"],
                   by=c("design", "comparison"))
  ##
  studies <- studies[, c("studlab", "design", "comparison", "treat1", "treat2",
                         "narms", "freq", "TE", "seTE",
                         "TE.dir", "seTE.dir", "TE.net", "seTE.net")]
  ##
  studies_lim <- studies[which(studies$narms==2),]
  studies_mult <- studies[which(studies$narms>2),]
  studies <- rbind(studies_lim[order(studies_lim$studlab),],
                   studies_mult[do.call(order, studies_mult[,c("studlab","treat1","treat2")]),])

  res <- list(n=n,
              k=x$k,
              d=length(unique(design$design)),
              trts=trts,
              comparisons=comparisons,
              studies=studies,
              direct=direct,
              network=network,
              design=design,
              multicomp=if (multiarm) multicomp else NULL,
              X.obs=X.obs,
              X.full=X.full,
              V=V,
              V.studies=V.studies,
              H=H,
              H.studies=H.studies)

  class(res) <- "nma.krahn"

  res
}
