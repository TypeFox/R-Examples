B.EmpDistrib <-function(net,n.seeds,n.neigh,sam.size=1,n.boot,otherNetParameters=FALSE){
  #sam.size (==1 for LSMI1) is the number of different samples taken from the network for each i and j
  #otherNetParameters is true if intervals and fallins for the rest of the parmeters
  #  (other than mean) are required.
  Obs.distrib.out <- w.p0s <- nw.p0sEkb <- nw.p0sEks <- as.list(rep(NA, length(n.seeds)*length(n.neigh)))
  counter <- 1
  for(i in n.seeds){
    for(j in n.neigh){

      if(j==0){
        n.dist <- 1
        # ^ n.dist is the number of different emp distr.
      }else{n.dist<-3}

      if(j==0){
        Obs.distrib<-Oempdegreedistrib(net, n.seeds=i, n.neigh=j, num.sam=sam.size)
        TMP <- Obs.distrib$seeds1
      }else{
        Obs.distrib<-Oempdegreedistrib(net, n.seeds=i,n.neigh=j, num.sam=sam.size, seeds=NULL)
      }
      Oparam<-OparametersEst(Obs.distrib)
      #B.distrib<-bootdeg(Obs.distrib, num.sam=sam.size,n.boot=n.boot)
      #return(B.distrib)
      #browser()
      tmp <- bootdeg(Obs.distrib, num.sam=sam.size,n.boot=n.boot)$empd[[1]]
      Obs.distrib.out[[counter]] <-Obs.distrib
      w.p0s[[counter]] <- tmp$empd.nw.p0sEkb
      nw.p0sEkb[[counter]] <- tmp$empd.nw.p0sEks
      nw.p0sEks[[counter]] <- tmp$empd.w.p0s
      counter <- counter+1
    }
    #return(B.distrib)
  }
  return(list(Obs.distrib.out=Obs.distrib.out, w.p0s=w.p0s, nw.p0sEkb=nw.p0sEkb, nw.p0sEks=nw.p0sEks))
}


estimable_k <- function(bootEmpD){
  res <- intersect(dimnames(bootEmpD$w.p0s[[1]])[[2]],dimnames(bootEmpD$w.p0s[[2]])[[2]])
  for(i in 3:length(bootEmpD$w.p0s)){
    res <- intersect(res,dimnames(bootEmpD$w.p0s[[i]])[[2]])
  }
  res
}

combineLSMINodes <- function(bootEmpD){
  #function combines the unodes(or seed1) from the elements of
  #Obs.distrib.out object, which is inside the bootEmpD list
  if("unodes"%in%names(bootEmpD$Obs.distrib.out[[1]])){
    nodes=bootEmpD$Obs.distrib.out[[1]]$unodes
  } else nodes=bootEmpD$Obs.distrib.out[[1]]$seeds1


  for(i in 2:length(bootEmpD$Obs.distrib.out)){
    if("unodes"%in%names(bootEmpD$Obs.distrib.out[[i]])){
      tmp=bootEmpD$Obs.distrib.out[[i]]$unodes
    } else tmp=bootEmpD$Obs.distrib.out[[i]]$seeds1
    nodes=c(nodes,tmp)
  }
  unlist(nodes)
}
closestCoverNDX <- function(x,coverage=.95){
  which(abs(x-coverage)==min(abs(x-coverage)),arr.ind=T)
}

# A function that sorts a matrix of tied optimal seed-wave combinations so that
# the first is the largest number of seeds among the smallest waves.
#
# The function takes a matrix containing tied optimal seed-wave combinations
# from the training proxy and sorts the rows of this matrix so that the
# largest number of seeds for the lowest number of waves is on top
# @param inMat is a matrix of tied seed-wave indices.
# @return A matrix with rows such that the indices for the largest number of
# seeds (first column) for the lowest number of waves (second column) is on top.
# @rdname sort_tied_opti
sort_tied_opti <- function(inMat){
  if(dim(inMat)[1] == 1)
    return(inMat)
  largestSeedFirst <- inMat[order(inMat[,1],decreasing = T),]
  outMat <- largestSeedFirst[order(largestSeedFirst[,2]),]
  outMat
}
# Alternative is to get the smallest sample size which is estimated with mean(degree)
# ##1.2. Find approximate sample size
# SS <- matrix(NA, length(n.neigh), length(n.seeds))
# SS[1,] <- n.seeds
# for(i in 1:5){
#       SS[(i+1),] <- SS[i,] + SS[1,]*MU*(MU-1)^(i-1)
# }
# #SS

#' A function that uses cross-validation to select seed-wave combination for
#' estimation of a degree's frequency.
#'
#' The function's inputs are a network, a vector of possible seed sample-sizes,
#' a vector of possible waves, and a few tuning parameters. The output will
#' contain the best seed-wave combination for each degree and the width of the
#' 95 percent bootstrap confidence intervals at each degree for
#' the best seed-wave combination.
#' @note Only one LSMI per seed-wave combination is currently supported.
#' @references Efron, B. (1979). Bootstrap methods: another look at the
#'  jackknife. The annals of Statistics, 1-26.
#' @references Thompson, M. E., Ramirez Ramirez, L. L., Lyubchich, V. and
#'  Gel, Y. R. (2015), Using the bootstrap for statistical inference
#'  on random graphs. Can J Statistics. doi: 10.1002/cjs.11271
#' @param network A network object that is list containing:
#'  \describe{
#'    \item{edges}{The edgelist of the network. A two column
#'      \code{matrix} where each row is an edge.}
#'    \item{degree}{The degree sequence of the network, which is
#'      an \code{integer} vector of length n.}
#'    \item{n}{The network order.}
#'  }
#'    The object can be created by \code{\link{local.network.MR.new5}} or
#'    it can be imported.
#' @param n.seeds A numeric vector for the different sample sizes of seed to use
#'  in cross-validation.
#' @param n.neigh A numeric vector for the different waves to use
#'  in cross-validation.
#' @param n.boot The number of bootstrap sample.
#' @param kmax The largest degree to preform cross-validation on.
#' @param proxyRep The number of time to sample a proxy. Default is 19.
#' @param proxyOrder The size of the proxy sample. Default is 30.
#' @return A list consisting of
#'  \item{selected_seed_wave}{A list of 3 matrices (one per estimation method.
#'    See supporting documentation \code{\link{bootdeg}}). Each matrix provides
#'    the best seed-wave combinations (obtained via cross-validation) for
#'    the respective estimation method.}
#'  \item{selected_seed_wave}{A list of 3 matrices (one per estimation method.
#'    See supporting documentation \code{\link{bootdeg}}). Each matrix provides
#'    the 95 percent bootstrap confidence intervals for the estimated degree frequency
#'    using the best seed-wave combinations (see above).}
#' @export
#' @examples
#' net <- artificial_networks[[1]]
#' a <- cross_validation(network = net, n.seeds = c(10, 20, 30), n.neigh = c(1, 2),
#'  n.boot = 200, kmax = 30)

cross_validation <- function(network, n.seeds, n.neigh, n.boot,
                             kmax, proxyRep = 19, proxyOrder = 30){
  sam.size = 1
  n.seeds <- sort(n.seeds)
  n.neigh <- sort(n.neigh)
  net_order <- network$n

    #make bootEmpD list for seed-wave combos
    bootEmpD=B.EmpDistrib(network,n.seeds,n.neigh,sam.size,n.boot)
    estimable_k_from_boot <- estimable_k(bootEmpD)
    used <- unique(combineLSMINodes(bootEmpD))
    count <- 1
    fallin.proxy.w.p0s <- fallin.proxy.nw.p0sEkb <- fallin.proxy.nw.p0sEks <-
      array(0, c(length(n.seeds), length(n.neigh), proxyRep, kmax))
    for(i in 1:length(n.seeds)){
      # i=1
      for(j in 1:length(n.neigh)){
        # j=1
        # build proxy from bootEmpD$Obs.empd.out
        tmp.w.p0s <- apply(bootEmpD$w.p0s[[count]], 2, stats::quantile, probs=c(0.025, 0.975))
        tmp.w.p0s <- tmp.w.p0s[,match(1:kmax,dimnames(tmp.w.p0s)[[2]], nomatch = NA)]
        dimnames(tmp.w.p0s)[[2]] <- 1:kmax

        tmp.nw.p0sEkb <- apply(bootEmpD$nw.p0sEkb[[count]], 2, stats::quantile, probs=c(0.025, 0.975))
        tmp.nw.p0sEkb <- tmp.nw.p0sEkb[,match(1:kmax,dimnames(tmp.nw.p0sEkb)[[2]], nomatch = NA)]
        dimnames(tmp.nw.p0sEkb)[[2]] <- 1:kmax

        tmp.nw.p0sEks <- apply(bootEmpD$nw.p0sEks[[count]], 2, stats::quantile, probs=c(0.025, 0.975))
        tmp.nw.p0sEks <- tmp.nw.p0sEks[,match(1:kmax,dimnames(tmp.nw.p0sEks)[[2]], nomatch = NA)]
        dimnames(tmp.nw.p0sEks)[[2]] <- 1:kmax

        tmp.w.p0s[is.na(tmp.w.p0s)] <- 0
        tmp.nw.p0sEkb[is.na(tmp.nw.p0sEkb)] <- 0
        tmp.nw.p0sEks[is.na(tmp.nw.p0sEks)] <- 0

        for(k in 1:proxyRep){
          # k=1
          proxyNodes <- sample(used, proxyOrder, replace = F)
          proxyDegrees <- network$degree[proxyNodes]
          proxyPMF <- table(proxyDegrees)/proxyOrder
          #trasform into vector with "0" where no obs occur
          PMFvector <- as.vector(proxyPMF[match(1:kmax, dimnames(proxyPMF)[[1]], nomatch=NA)])
          PMFvector[is.na(PMFvector)] <- 0
          #take first 5
          firstK <- PMFvector[1:kmax]

          fallin.proxy.w.p0s[i, j, k, ] <- ((tmp.w.p0s[1, 1:kmax]<firstK) & (firstK<tmp.w.p0s[2, 1:kmax]))
          fallin.proxy.nw.p0sEkb[i, j, k, ] <- ((tmp.nw.p0sEkb[1, 1:kmax]<firstK) & (firstK<tmp.nw.p0sEkb[2, 1:kmax]))
          fallin.proxy.nw.p0sEks[i, j, k, ] <- ((tmp.nw.p0sEks[1, 1:kmax]<firstK) & (firstK<tmp.nw.p0sEks[2, 1:kmax]))
        }
        count <- count+1
      }
    }
    coverage.proxy.w.p0s <- apply(fallin.proxy.w.p0s, c(1, 2, 4), mean, na.rm=T)
    coverage.proxy.nw.p0sEkb <- apply(fallin.proxy.nw.p0sEkb, c(1, 2, 4), mean, na.rm=T)
    coverage.proxy.nw.p0sEks <- apply(fallin.proxy.nw.p0sEks, c(1, 2, 4), mean, na.rm=T)

    opti.cover.w.p0s <- apply(coverage.proxy.w.p0s, 3, closestCoverNDX)
    opti.cover.nw.p0sEkb <- apply(coverage.proxy.nw.p0sEkb, 3, closestCoverNDX)
    opti.cover.nw.p0sEks <- apply(coverage.proxy.nw.p0sEks, 3, closestCoverNDX)
    #output Matrices

    opti.CI.w.p0s <- opti.CI.nw.p0sEkb <- opti.CI.nw.p0sEks <-
      matrix(nrow = 2, ncol = length(estimable_k_from_boot),
             dimnames = list(c("LB", "UB"), estimable_k_from_boot))

    opti.seed_wave.w.p0s <- opti.seed_wave.nw.p0sEkb <- opti.seed_wave.nw.p0sEks <-
      matrix(nrow = 2, ncol = length(estimable_k_from_boot),
             dimnames = list(c("Seeds", "Waves"), estimable_k_from_boot))
    #test if matrix or list
    if (is.list(opti.cover.w.p0s)){

      for(kDegree in estimable_k_from_boot){
        kDegree_num <- as.numeric(kDegree)
        sortedOpti.cover.w.p0s <- sort_tied_opti(opti.cover.w.p0s[[kDegree_num]])
        #browser()
        optimalSeedNDX <-  sortedOpti.cover.w.p0s[1, ][1]
        optimalNeighNDX <- sortedOpti.cover.w.p0s[1, ][2]
        opti.seed_wave.w.p0s[1, kDegree] <- n.seeds[optimalSeedNDX]
        opti.seed_wave.w.p0s[2, kDegree] <- n.neigh[optimalNeighNDX]
        listLocation <- (optimalSeedNDX-1)*length(n.neigh)+optimalNeighNDX
        opti.CI.w.p0s[,kDegree] <- stats::quantile(bootEmpD$w.p0s[[listLocation]][, kDegree], probs=c(0.025, 0.975))

      }
    }else if(is.matrix(opti.cover.w.p0s)){
      for(kDegree in estimable_k_from_boot){
        kDegree_num <- as.numeric(kDegree)
        optimalSeedNDX = opti.cover.w.p0s[1, kDegree_num]
        optimalNeighNDX <- opti.cover.w.p0s[2, kDegree_num]
        opti.seed_wave.w.p0s[1, kDegree] <- n.seeds[optimalSeedNDX]
        opti.seed_wave.w.p0s[2, kDegree] <- n.neigh[optimalNeighNDX]
        listLocation <- (optimalSeedNDX-1)*length(n.neigh)+optimalNeighNDX
        opti.CI.w.p0s[,kDegree] <- stats::quantile(bootEmpD$w.p0s[[listLocation]][, kDegree], probs=c(0.025, 0.975))
        }
    }else print("unknown data type output from optimal function")


    if (is.list(opti.cover.nw.p0sEkb)){

      for(kDegree in estimable_k_from_boot){
        kDegree_num <- as.numeric(kDegree)
        sortedOpti.cover.nw.p0sEkb <- sort_tied_opti(opti.cover.nw.p0sEkb[[kDegree_num]])
        #browser()
        optimalSeedNDX <-  sortedOpti.cover.nw.p0sEkb[1, ][1]
        optimalNeighNDX <- sortedOpti.cover.nw.p0sEkb[1, ][2]
        opti.seed_wave.nw.p0sEkb[1, kDegree] <- n.seeds[optimalSeedNDX]
        opti.seed_wave.nw.p0sEkb[2, kDegree] <- n.neigh[optimalNeighNDX]
        listLocation <- (optimalSeedNDX-1)*length(n.neigh)+optimalNeighNDX
        opti.CI.nw.p0sEkb[,kDegree] <- stats::quantile(bootEmpD$nw.p0sEkb[[listLocation]][, kDegree], probs=c(0.025, 0.975))

      }
    }else if(is.matrix(opti.cover.nw.p0sEkb)){
      for(kDegree in estimable_k_from_boot){
        kDegree_num <- as.numeric(kDegree)
        optimalSeedNDX = opti.cover.nw.p0sEkb[1, kDegree_num]
        optimalNeighNDX <- opti.cover.nw.p0sEkb[2, kDegree_num]
        opti.seed_wave.nw.p0sEkb[1, kDegree] <- n.seeds[optimalSeedNDX]
        opti.seed_wave.nw.p0sEkb[2, kDegree] <- n.neigh[optimalNeighNDX]
        listLocation <- (optimalSeedNDX-1)*length(n.neigh)+optimalNeighNDX
        opti.CI.nw.p0sEkb[,kDegree] <- stats::quantile(bootEmpD$nw.p0sEkb[[listLocation]][, kDegree], probs=c(0.025, 0.975))
      }
    }else print("unknown data type output from optimal function")


    if (is.list(opti.cover.nw.p0sEks)){

      for(kDegree in estimable_k_from_boot){
        kDegree_num <- as.numeric(kDegree)
        sortedOpti.cover.nw.p0sEks <- sort_tied_opti(opti.cover.nw.p0sEks[[kDegree_num]])
        #browser()
        optimalSeedNDX <-  sortedOpti.cover.nw.p0sEks[1, ][1]
        optimalNeighNDX <- sortedOpti.cover.nw.p0sEks[1, ][2]
        opti.seed_wave.nw.p0sEks[1, kDegree] <- n.seeds[optimalSeedNDX]
        opti.seed_wave.nw.p0sEks[2, kDegree] <- n.neigh[optimalNeighNDX]
        listLocation <- (optimalSeedNDX-1)*length(n.neigh)+optimalNeighNDX
        opti.CI.nw.p0sEks[,kDegree] <- stats::quantile(bootEmpD$nw.p0sEks[[listLocation]][, kDegree], probs=c(0.025, 0.975))

      }
    }else if(is.matrix(opti.cover.nw.p0sEks)){
      for(kDegree in estimable_k_from_boot){
        kDegree_num <- as.numeric(kDegree)
        optimalSeedNDX = opti.cover.nw.p0sEks[1, kDegree_num]
        optimalNeighNDX <- opti.cover.nw.p0sEks[2, kDegree_num]
        opti.seed_wave.nw.p0sEks[1, kDegree] <- n.seeds[optimalSeedNDX]
        opti.seed_wave.nw.p0sEks[2, kDegree] <- n.neigh[optimalNeighNDX]
        listLocation <- (optimalSeedNDX-1)*length(n.neigh)+optimalNeighNDX
        opti.CI.nw.p0sEks[,kDegree] <- stats::quantile(bootEmpD$w.p0s[[listLocation]][, kDegree], probs=c(0.025, 0.975))
      }
    }else print("unknown data type output from optimal function")
    CI_selected_seed_wave <- list(opti.CI.w.p0s = opti.CI.w.p0s,
                                  opti.CI.nw.p0sEkb = opti.CI.nw.p0sEkb,
                                  opti.CI.nw.p0sEks = opti.CI.nw.p0sEks)
    selected_seed_wave <-
      list(opti.seed_wave.w.p0s = opti.seed_wave.w.p0s,
           opti.seed_wave.nw.p0sEkb = opti.seed_wave.nw.p0sEkb,
             opti.seed_wave.nw.p0sEks = opti.seed_wave.nw.p0sEks)
    res <- list(selected_seed_wave = selected_seed_wave,
                CI_selected_seed_wave = CI_selected_seed_wave)
    res
}
