calinhara <- function(x,clustering,cn=max(clustering)){
  x <- as.matrix(x)
  p <- ncol(x)
  n <- nrow(x)
  cln <- rep(0,cn)
  W <- matrix(0,p,p)
  for (i in 1:cn)
    cln[i] <- sum(clustering==i)
#  print(cln)
  for (i in 1:cn) {
    clx <- x[clustering==i,]
    cclx <- cov(as.matrix(clx))
#    print(cclx)
    if (cln[i] < 2) 
            cclx <- 0 
    W <- W + ((cln[i] - 1) * cclx)
  }
  S <- (n - 1) * cov(x)
  B <- S - W
  out <- (n-cn)*sum(diag(B))/((cn-1)*sum(diag(W)))
  out
}

dudahart2 <- function(x,clustering,alpha=0.001){
  x <- as.matrix(x)
  p <- ncol(x)
  n <- nrow(x)
  cln <- rep(0,2)
  W <- matrix(0,p,p)
  for (i in 1:2)
    cln[i] <- sum(clustering==i)
#  print(cln)
  for (i in 1:2) {
    clx <- x[clustering==i,]
    cclx <- cov(as.matrix(clx))
#    print(cclx)
    if (cln[i] < 2) 
            cclx <- 0 
    W <- W + ((cln[i] - 1) * cclx)
  }
  W1 <- (n-1)*cov(as.matrix(x))
  dh <- sum(diag(W))/sum(diag(W1))
  z <- qnorm(1-alpha)
  compare <- 1-2/(pi*p)-z*sqrt(2*(1-8/(pi^2*p))/(n*p))
  qz <- (-dh+1-2/(pi*p))/sqrt(2*(1-8/(pi^2*p))/(n*p))
  p.value <- 1-pnorm(qz)
  cluster1 <- dh>=compare
  out <- list(p.value=p.value,
              dh=dh,compare=compare,cluster1=cluster1,alpha=alpha,z=z)
  out
}
 

# # L1 equivalent to calinhara, large is good
# clarach <- function(x,claraout=NULL,partition=NULL,clara1=NULL){
#   require(cluster)
#   n <- nrow(x)
#   if (is.null(clara1)) clara1 <- clara(x,1)
#   if (is.null(claraout)){
#     cn <- max(partition)
#     co <- 0
#     for (i in 1:cn){
#       clara11 <- clara(x[partition==i,],1)
#       co <- co+sum(partition==i)*clara11$objective
# #      print("part")
# #      print(clara11$objective)
# #      print(mean(dist(x[partition==i,])))
# #      print(co)
#     }
#   }
#   else{  
#     co <- n*claraout$objective
#     cn <- max(claraout$clustering)
# #    print("co")
# #    print(co)
# #    print(mean(dist(x[claraout$clustering==1])))
#   }
#   c1 <- n*clara1$objective
#   out <- (n-cn)*(c1-co)^2/((cn-1)*co^2)
# #  print(c1)
# #  print(out)
#   out
# }


kmeansruns <- function(data,krange=2:10,criterion="ch",
                       iter.max=100,runs=100,
                       scaledata=FALSE,alpha=0.001,
                       critout=FALSE,plot=FALSE,...){
  data <- as.matrix(data)
  if (criterion=="asw") sdata <- dist(data)
  if (scaledata) data <- scale(data)
  cluster1 <- 1 %in% krange
  crit <- numeric(max(krange))
  km <- list()
  for (k in krange){
    if (k>1){
      minSS <- Inf
      kmopt <- NULL
      for (i in 1:runs){
        options(show.error.messages = FALSE)
        repeat{
#          cat(k," ",i,"before \n")
          kmm <- try(kmeans(data,k,iter.max=iter.max,...))
#          str(kmm)
          if (class(kmm) != "try-error") break
#         cat(k," ",i,"\n")
        }
        options(show.error.messages = TRUE)
        swss <- sum(kmm$withinss)
#        print(calinhara(data,kmm$cluster))
        if (swss<minSS){
          kmopt <- kmm
          minSS <- swss
        }
        if (plot){
          par(ask=TRUE)
          pairs(data,col=kmm$cluster,main=swss)
        }
      } # for i
      km[[k]] <- kmopt
#      print(km[[k]])
#      print(calinhara(data,km[[k]]$cluster))
      crit[k] <- switch(criterion,
             asw=cluster.stats(sdata,km[[k]]$cluster)$avg.silwidth,
             ch=calinhara(data,km[[k]]$cluster))
      if (critout)
          cat(k," clusters ",crit[k],"\n")
    } # if k>1
  } # for k
  if (cluster1)
    cluster1 <- dudahart2(data,km[[2]]$cluster,alpha=alpha)$cluster1
  k.best <- which.max(crit)
  if (cluster1)
    k.best <- 1
#  print(crit)
#  print(k.best)
#  print(km[[k.best]])
  km[[k.best]]$crit <- crit
  km[[k.best]]$bestk <- k.best
  out <- km[[k.best]]
  out 
}


pamk <- function (data, krange = 2:10, criterion = "asw", usepam = TRUE, 
    scaling = FALSE, alpha = 0.001, diss = inherits(data, "dist"), 
    critout = FALSE, ns=10, seed=NULL, ...) 
{
    ddata <- as.matrix(data)
#    require(cluster)
    if (!identical(scaling,FALSE))
        sdata <- scale(ddata, scale = scaling)
    else sdata <- ddata
    cluster1 <- 1 %in% krange
    critval <- numeric(max(krange))
    pams <- list()
    for (k in krange) {
        if (usepam) 
            pams[[k]] <- pam(sdata, k, diss=diss, ...)
        else pams[[k]] <- clara(sdata, k, ...)
        if (k != 1) 
            critval[k] <- switch(criterion, asw = pams[[k]]$silinfo$avg.width,
                             multiasw=distcritmulti(sdata,pams[[k]]$clustering,seed=seed,ns=ns)$crit.overall,
                ch = ifelse(diss, cluster.stats(sdata, pams[[k]]$clustering)$ch, 
                  calinhara(sdata, pams[[k]]$clustering)))
        if (critout) 
            cat(k, " clusters ", critval[k], "\n")
    }
    k.best <- (1:max(krange))[which.max(critval)]
    if (cluster1){
      if(diss)
        cluster1 <- FALSE
      else{
        cxx <- dudahart2(sdata, pams[[2]]$clustering, alpha = alpha)
        critval[1] <- cxx$p.value
        cluster1 <- cxx$cluster1
      }
    }
    if (cluster1) 
        k.best <- 1
    out <- list(pamobject = pams[[k.best]], nc = k.best, crit=critval)
    out
}



cluster.stats <- function (d = NULL, clustering, alt.clustering = NULL,
                           noisecluster=FALSE,
                              silhouette = TRUE, G2 = FALSE, G3 = FALSE,
                              wgap=TRUE, sepindex=TRUE, sepprob=0.1,
                              sepwithnoise=TRUE,
                              compareonly = FALSE,
                              aggregateonly = FALSE)
{
    if (!is.null(d)) 
        d <- as.dist(d)
    cn <- max(clustering)
    clusteringf <- as.factor(clustering)
    clusteringl <- levels(clusteringf)
    cnn <- length(clusteringl)
    if (cn != cnn) {
        warning("clustering renumbered because maximum != number of clusters")
        for (i in 1:cnn) clustering[clusteringf == clusteringl[i]] <- i
        cn <- cnn
    }
    n <- length(clustering)
    noisen <- 0
    cwn <- cn
    if (noisecluster){
      noisen <- sum(clustering==cn)
      cwn <- cn-1
    }
    diameter <- average.distance <- median.distance <- separation <- average.toother <- cluster.size <- within.dist <- between.dist <- numeric(0)
    for (i in 1:cn) cluster.size[i] <- sum(clustering == i)
    pk1 <- cluster.size/n
    pk10 <- pk1[pk1 > 0]
    h1 <- -sum(pk10 * log(pk10))
    corrected.rand <- vi <- NULL
    if (!is.null(alt.clustering)) {
        choose2 <- function(v) {
            out <- numeric(0)
            for (i in 1:length(v)) out[i] <- ifelse(v[i] >= 2, 
                choose(v[i], 2), 0)
            out
        }
        cn2 <- max(alt.clustering)
        clusteringf <- as.factor(alt.clustering)
        clusteringl <- levels(clusteringf)
        cnn2 <- length(clusteringl)
        if (cn2 != cnn2) {
            warning("alt.clustering renumbered because maximum != number of clusters")
            for (i in 1:cnn2) alt.clustering[clusteringf == clusteringl[i]] <- i
            cn2 <- cnn2
        }
        nij <- table(clustering, alt.clustering)
        dsum <- sum(choose2(nij))
        cs2 <- numeric(0)
        for (i in 1:cn2) cs2[i] <- sum(alt.clustering == i)
        sum1 <- sum(choose2(cluster.size))
        sum2 <- sum(choose2(cs2))
        pk2 <- cs2/n
        pk12 <- nij/n
        corrected.rand <- (dsum - sum1 * sum2/choose2(n))/((sum1 + 
            sum2)/2 - sum1 * sum2/choose2(n))
        pk20 <- pk2[pk2 > 0]
        h2 <- -sum(pk20 * log(pk20))
        icc <- 0
        for (i in 1:cn) for (j in 1:cn2) if (pk12[i, j] > 0) 
            icc <- icc + pk12[i, j] * log(pk12[i, j]/(pk1[i] * 
                pk2[j]))
        vi <- h1 + h2 - 2 * icc
    }
    if (compareonly) {
        out <- list(corrected.rand = corrected.rand, vi = vi)
    }
    else {
#        if (silhouette) 
#            require(cluster)
        dmat <- as.matrix(d)
        within.cluster.ss <- 0
        overall.ss <- nonnoise.ss <- sum(d^2)/n
        if (noisecluster)
          nonnoise.ss <- sum(as.dist(dmat[clustering<=cwn,
                                          clustering<=cwn])^2)/
                                            sum(clustering<=cwn)
        ave.between.matrix <-
          separation.matrix <- matrix(0, ncol = cn, nrow = cn)
        di <- list()
        for (i in 1:cn) {
            cluster.size[i] <- sum(clustering == i)
            di <- as.dist(dmat[clustering == i, clustering == 
                i])
            if (i<=cwn){
              within.cluster.ss <- within.cluster.ss + sum(di^2)/cluster.size[i]
              within.dist <- c(within.dist, di)
            }
            if (length(di) > 0) 
                diameter[i] <- max(di)
            else diameter[i] <- NA
            average.distance[i] <- mean(di)
            median.distance[i] <- median(di)
            bv <- numeric(0)
            for (j in 1:cn) {
                if (j != i) {
                  sij <- dmat[clustering == i, clustering == 
                    j]
                  bv <- c(bv, sij)
                  if (i < j) {
                    separation.matrix[i, j] <- separation.matrix[j, 
                      i] <- min(sij)
                    ave.between.matrix[i, j] <- ave.between.matrix[j, i] <-
                        mean(sij)
                    if (i<=cwn & j<=cwn)
                      between.dist <- c(between.dist, sij)
                  }
                }
            }
            separation[i] <- min(bv)
            average.toother[i] <- mean(bv)
        }
        average.between <- mean(between.dist)
        average.within <- mean(within.dist)
        nwithin <- length(within.dist)
        nbetween <- length(between.dist)
        between.cluster.ss <- nonnoise.ss - within.cluster.ss
        ch <- between.cluster.ss * (n - noisen - cwn)/(within.cluster.ss * 
            (cwn - 1))
        clus.avg.widths <- avg.width <- NULL
        if (silhouette) {
            sii <- silhouette(clustering, dmatrix = dmat)
            sc <- summary(sii)
            clus.avg.widths <- sc$clus.avg.widths
            if (noisecluster)
              avg.width <- mean(sii[clustering<=cwn,3])
            else
              avg.width <- sc$avg.width
        }
        g2 <- g3 <- cn2 <- cwidegap <- widestgap <- sindex <- NULL
        if (G2) {
            splus <- sminus <- 0
            for (i in 1:nwithin) {
                splus <- splus + sum(within.dist[i] < between.dist)
                sminus <- sminus + sum(within.dist[i] > between.dist)
            }
            g2 <- (splus - sminus)/(splus + sminus)
        }
        if (G3) {
            sdist <- sort(c(within.dist, between.dist))
            sr <- nwithin + nbetween
            dmin <- sum(sdist[1:nwithin])
            dmax <- sum(sdist[(sr - nwithin + 1):sr])
            g3 <- (sum(within.dist) - dmin)/(dmax - dmin)
        }
        pearsongamma <- cor(c(within.dist, between.dist), c(rep(0, 
            nwithin), rep(1, nbetween)))
        dunn <- min(separation[1:cwn])/max(diameter[1:cwn],na.rm=TRUE)
        acwn <- ave.between.matrix[1:cwn,1:cwn]
        dunn2 <- min(acwn[upper.tri(acwn)])/
          max(average.distance[1:cwn],na.rm=TRUE)
        if (wgap){
          cwidegap <- rep(0,cwn)
          for (i in 1:cwn)
            if (sum(clustering==i)>1)
              cwidegap[i] <- max(hclust(as.dist(dmat[clustering==i,
                                                   clustering==i]),
                                      method="single")$height)
          widestgap <- max(cwidegap)
        }
        if (sepindex){
          psep <- rep(NA,n)
          if (sepwithnoise | !noisecluster){
            for (i in 1:n)
              psep[i] <- min(dmat[i,clustering!=clustering[i]])
            minsep <- floor(n*sepprob)
          }
          else{
            dmatnn <- dmat[clustering<=cwn,clustering<=cwn]
            clusteringnn <- clustering[clustering<=cwn]
            for (i in 1:(n-noisen))
              psep[i] <- min(dmatnn[i,clusteringnn!=clusteringnn[i]])
            minsep <- floor((n-noisen)*sepprob)
          }
          sindex <- mean(sort(psep)[1:minsep])
        }
        if (!aggregateonly)
          out <- list(n = n, cluster.number = cn, cluster.size = cluster.size,
                      min.cluster.size = min(cluster.size[1:cwn]),
                      noisen=noisen,
            diameter = diameter, average.distance = average.distance, 
            median.distance = median.distance, separation = separation, 
            average.toother = average.toother, separation.matrix = separation.matrix,
                      ave.between.matrix=ave.between.matrix,
            average.between = average.between, average.within = average.within, 
            n.between = nbetween, n.within = nwithin,
                      max.diameter = max(diameter[1:cwn],na.rm=TRUE),
                      min.separation = sepwithnoise*min(separation)+
                      (!sepwithnoise)*min(separation[1:cwn]),
                      within.cluster.ss = within.cluster.ss, 
            clus.avg.silwidths = clus.avg.widths, avg.silwidth = avg.width, 
            g2 = g2, g3 = g3, pearsongamma = pearsongamma, dunn = dunn,
                      dunn2=dunn2,
            entropy = h1, wb.ratio = average.within/average.between, 
            ch = ch, cwidegap=cwidegap, widestgap=widestgap,
                    sindex=sindex,
                    corrected.rand = corrected.rand, vi = vi)
        else
          out <- list(n = n, cluster.number = cn,
                      min.cluster.size = min(cluster.size[1:cwn]),
                      noisen=noisen,
            average.between = average.between, average.within = average.within, 
                      max.diameter = max(diameter[1:cwn],na.rm=TRUE),
                      min.separation = sepwithnoise*min(separation)+
                      (!sepwithnoise)*min(separation[1:cwn]), 
            ave.within.cluster.ss = within.cluster.ss/(n-noisen), 
            avg.silwidth = avg.width, 
            g2 = g2, g3 = g3, pearsongamma = pearsongamma, dunn = dunn,
                      dunn2=dunn2,
            entropy = h1, wb.ratio = average.within/average.between, 
            ch = ch, widestgap=widestgap,
                    sindex=sindex,
                    corrected.rand = corrected.rand, vi = vi)
    }
    out
}

   
distcritmulti <- function(x,clustering,part=NULL,ns=10,criterion="asw",
                    fun="dist",metric="euclidean",
                     count=FALSE,seed=NULL,...){
  if (!is.null(seed)) set.seed(seed)
  n <- length(clustering)
  if (is.null(part)){
    pn1 <- n %/% ns
    pn2 <- n %% ns
    part <- rep(pn1,ns)
    part[ns] <- part[ns]+pn2
  }
  np <- sum(part)
  ns <- length(part)
  n <- length(clustering)
  npsam <- sample(n,np)
  cp <- cumsum(part)
  ss <- list()
  ss[[1]] <- npsam[1:cp[1]]
  asw <- numeric(0)
  for (i in 2:ns)
    ss[[i]] <- npsam[(cp[i-1]+1):cp[i]]
  for (i in 1:ns){
    if (count) cat("Subset ",i,"\n")
    if (fun=="dist")
      dx <- dist(x[ss[[i]],],method=metric,...)
    else
      dx <- daisy(x[ss[[i]],],metric=metric,...)
    if (criterion=="asw")
      asw[i] <- summary(silhouette(clustering[ss[[i]]],dx))$avg.width
    if (criterion=="pearsongamma")
      asw[i] <- cluster.stats(dx,clustering[ss[[i]]],silhouette=FALSE)$pearsongamma
  }
  aswav <- sum(part*asw)/np
  crit.sd <- sd(asw)
  out <- list(crit.overall=aswav,crit.sub=asw,crit.sd=crit.sd,subsets=ss)
  out
}
