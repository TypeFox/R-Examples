# # d is a distance object or matrix, clustering is assumed to be numbered
# # 1 to clusternumber
# # If alt.clustering is another clustering, the corrected rand index will
# # be computed.
# # silhouette, G2 and G3 indicate if the corresponding statistics are
# # computed. Silhouette requires library cluster, G2 and G3 may take very long
# # for large n.
# cluster.stats <- function(d,clustering,alt.clustering=NULL,
#                           silhouette=TRUE,G2=FALSE,G3=FALSE,
#                           compareonly=FALSE){
#   cn <- max(clustering)
#   n <- length(clustering)
#   diameter <- average.distance <- median.distance <- separation <-
#       average.toother <- 
#       cluster.size <- within.dist <- between.dist <- numeric(0)
#   for (i in 1:cn)
#     cluster.size[i] <- sum(clustering==i)
#   pk1 <- cluster.size/n
#   pk10 <- pk1[pk1>0]
#   h1 <- -sum(pk10*log(pk10))
#   corrected.rand <- vi <- NULL
#   if (!is.null(alt.clustering)){
#     choose2 <- function(v){
#       out <- numeric(0)
#       for (i in 1:length(v))
#         out[i] <- ifelse(v[i]>=2,choose(v[i],2),0)
#       out
#     }
#     cn2 <- max(alt.clustering)
#     nij <- table(clustering,alt.clustering)
#     dsum <- sum(choose2(nij))
#     cs2 <- numeric(0)
#     for (i in 1:cn2)
#       cs2[i] <- sum(alt.clustering==i)
#     sum1 <- sum(choose2(cluster.size))
#     sum2 <- sum(choose2(cs2))
#     pk2 <- cs2/n
#     pk12 <- nij/n
#     corrected.rand <- (dsum-sum1*sum2/choose2(n))/
#       ((sum1+sum2)/2-sum1*sum2/choose2(n))
#     pk20 <- pk2[pk2>0]
#     h2 <- -sum(pk20*log(pk20))
#     icc <- 0
#     for (i in 1:cn)
#       for (j in 1:cn2)
#         if (pk12[i,j]>0)
#           icc <- icc+pk12[i,j]*log(pk12[i,j]/(pk1[i]*pk2[j]))
# #    print(icc)
#     vi <- h1+h2-2*icc 
#   }
#   if (compareonly){
#     out <- list(corrected.rand=corrected.rand,vi=vi)
#   }
#   else{
#     if (silhouette) require(cluster)
#     dmat <- as.matrix(d)
#     within.cluster.ss <- 0
#     separation.matrix <- matrix(0,ncol=cn,nrow=cn)
#     di <- list()
#     for (i in 1:cn){
#       cluster.size[i] <- sum(clustering==i)
#       di <- as.dist(dmat[clustering==i,clustering==i])
#       within.cluster.ss <- within.cluster.ss+sum(di^2)/cluster.size[i]
#       within.dist <- c(within.dist,di)
#       if (sum(clustering==i)>1)
#         diameter[i] <- max(di)
#       else
#         diameter[i] <- 0        
#       average.distance[i] <- mean(di)
#       median.distance[i] <- median(di)
#       bv <- numeric(0)
#       for (j in 1:cn){
#         if (j!=i){
#           sij <- dmat[clustering==i,clustering==j]
#           bv <- c(bv,sij)
#           if (i<j){
#             separation.matrix[i,j] <- separation.matrix[j,i] <- min(sij)
#             between.dist <- c(between.dist,sij)
#           }
#         }
#       }
#       separation[i] <- min(bv)
#       average.toother[i] <- mean(bv)
#     }
#     average.between <- mean(between.dist)
#     average.within <- mean(within.dist)
#     nwithin <- length(within.dist)
#     nbetween <- length(between.dist)
#     clus.avg.widths <- avg.width <- NULL
#     if (silhouette){
#       sc <- summary(silhouette(clustering,dmatrix=dmat))
#       clus.avg.widths <- sc$clus.avg.widths
#       avg.width <- sc$avg.width
#     }
#     g2 <- g3 <- cn2 <- NULL
#     if (G2){
#           splus <- sminus <- 0
#           for (i in 1:nwithin) {
#              splus  <- splus  + sum(within.dist[i]<between.dist)
#              sminus <- sminus + sum(within.dist[i]>between.dist) 
#           }
#           g2 <- (splus - sminus)/(splus + sminus)
#     }
#     if (G3){
#       sdist <- sort(c(within.dist,between.dist))
#       sr <- nwithin+nbetween
#       dmin <- sum(sdist[1:nwithin])
#       dmax <- sum(sdist[(sr-nwithin+1):sr])
#       g3 <- (sum(within.dist)-dmin)/(dmax-dmin)
#     }
#     hubertgamma <- cor(c(within.dist,between.dist),c(rep(0,nwithin),
#                                                      rep(1,nbetween)))
#     dunn <- min(separation)/max(diameter)
#     out <- list(n=n,
#                 cluster.number=cn,
#                 cluster.size=cluster.size, # vector of cluster sizes
#                 diameter=diameter, # vector of cluster diameters
#                 average.distance=average.distance,
#                                           # vector of within cl. av. dist.
#                 median.distance=median.distance,
#                                           # vector of within cl. median dist.
#                 separation=separation, # vector of min. clusterwise between dist.
#                 average.toother=average.toother, 
#                                        # vector of mean clusterwise between dist.
#                 separation.matrix=separation.matrix,
#                                        # clusterwise matrix of min. between dist.
#                 average.between=average.between, # mean between cl. distance
#                 average.within=average.within, # mean within cl. distance
#                 n.between=nbetween, # number of between cl. distances
#                 n.within=nwithin, # number of within cl. distances
#                 within.cluster.ss=within.cluster.ss,
#                 clus.avg.silwidths=clus.avg.widths,
#                                   # vector of cluster avg. silhouette widths
#                 avg.silwidth=avg.width, # average silhouette width
#                 g2=g2, # Goodman and Kruskal coefficient, see Gordon p. 62
#                 g3=g3, # G3 index, see Gordon p. 62
#                 hubertgamma=hubertgamma, # Correlation between distances and
#                                           # 0-1-vector same/different cluster
#                 dunn=dunn, # Dunn index, see Halkidi et al. (2002)
#                            # Min. sepatation / max. diameter
#                 entropy=h1,
#                 wb.ratio=average.within/average.between,
#                 corrected.rand=corrected.rand, vi=vi) # Corrected rand index between
#                                           # clustering and alt.clustering
#   #  class(out) <- "cluster.stats"
#   }
#   out
# }
# 
# 
#
# face benchmark dataset by M. Maechler and C. Hennig
#
## MM: -  function(n, p)  {where `n' was a bit tricky}
##     -  return grouping as well

rFace <- function(n, p = 6, nrep.top = 2, smile.coef = 0.6,
                  dMoNo = 1.2, dNoEy = 1)
{
    ## Purpose: Generate random "Face" data set -- to be "Hard for Clustering"
    ## -------------------------------------------------------------------------
    ## Arguments: (n,p)    : dimension of result
    ##            nrep.top : #{repetitions} of the top point / "hair tip"
    ##            dMoNo    : distance{Mouth, Nose}
    ##            dNoEy    : distance{Nose,  Eyes} {vertically only}
    ##
    ## MM- TODO's:   o  provide more arguments (face shape)
    ##			    --> separation of clusters
    ## -------------------------------------------------------------------------
    ## Author: Christian Hennig & Martin Maechler, 26 Jun 2002

    if((p <- as.integer(p)) < 2) stop("number of variables p must be at least 2")
    if((n <- as.integer(n)) < 10) stop("number of points  n  must be at least 10")
    if((nrep.top <- as.integer(nrep.top)) < 1) stop("`nrep.top' must be positive")

    ntips <- nrep.top + 2
    n0 <- n - ntips
    m <- n0 %/% 5
    n.5 <- n0 %/% 2
    n.7 <- n.5+m
    n.9 <- n.7+m
    ## shouldn't happen:
    if(m < 1 || n.9 >= n0) stop("number of points n is too small")

    ## Indices of the different groups :
    Gr <- list(chin = 1:m,
               mouth= (m+1):n.5,
               nose = (n.5+1):n.7,
               rEye = (n.7+1):n.9,
               lEye = (n.9+1):n0,
               tips = (n0+1):n)

    face <- matrix(nrow = n, ncol = p)
    ## chin :
    face[Gr$chin, 1] <- U <- runif(m, -3, 3)
    face[Gr$chin, 2] <- rnorm(m, mean = U^2, sd=0.1)
    ## mouth:
    m0m <- 3 # lower mouth mean
    face[Gr$mouth, 1] <- Z <- rnorm(n.5-m, sd= 0.5)
    face[Gr$mouth, 2] <- rnorm(n.5-m, mean = m0m + smile.coef * Z^2, sd= 0.2)
    ## nose :
    face[Gr$nose, 1] <-  rnorm(m, mean=0, sd=0.2)
    yEye <- 17
    ##face[Gr$nose, 2] <- rnorm(m, mean=9, sd= 2.5)
    face[Gr$nose, 2] <- pmin(yEye - dNoEy,
                             (m0m + dMoNo) + rgamma(m, shape= 0.8, scale= 2.5))
    ## right eye  U[circle] :
    rangle <- runif(m, 0, 2*pi)
    rpos   <- runif(m)
    face[Gr$rEye, 1] <- rpos*cos(rangle) +  2
    face[Gr$rEye, 2] <- rpos*sin(rangle) + yEye
    ## left eye:
    face[Gr$lEye, 1] <- rnorm(n0-n.9, mean= -2,   sd= 0.5)
    face[Gr$lEye, 2] <- rnorm(n0-n.9, mean= yEye, sd= 0.5)
    ## `ntips'  ``hair tips'', the last one `nrep.top' times:
    face[n0+1,1:2] <- c(-4.5, 25)
    face[n0+2,1:2] <- c( 4.5, 25)
    for(k in 1:nrep.top)
        face[n-k+1, 1:2] <- c(0,32)

    ##-- Extra coordinates with noise ---

    if(p >= 3) {
        face[, p] <- rexp(n)
        if(p >= 4) {
            face[, p-1] <- rt(n,df=1)
            if(p >= 5)
                for(k in 3:(p-2))
                    face[,k] <- rnorm(n)
        }
    }
    gr <- character(n)
    ng <- names(Gr)
    for(i in seq(ng)) gr[Gr[[i]]] <- ng[i]
    structure(face, grouping = as.factor(gr), indexlist = Gr)
}

