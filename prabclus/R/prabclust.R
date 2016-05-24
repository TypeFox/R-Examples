
# Reorders components of clustsummary
prabclust <- function (prabobj, mdsmethod = "classical", mdsdim = 4,
                       nnk = ceiling(prabobj$n.species/40), 
    nclus = 0:9, modelid = "all", permutations = 0) 
{
#    require(MASS)
#    require(mclust)
    oregions <- order(prabobj$specperreg)
    prabo1 <- prabobj$prab[oregions, ]
    ospecies <- do.call("order", as.data.frame(t(prabo1)))
    dm <- prabobj$distmat[ospecies, ospecies]
    if (mdsmethod != "classical") {
        mindm <- min(dm[dm > 0])/10
        for (i in 1:(prabobj$n.species - 1)) for (j in (i + 1):prabobj$n.species) if (dm[i, 
            j] < mindm) 
            dm[i, j] <- dm[j, i] <- mindm
    }
#    print(dm)
#    print(mdsdim)
    mdsbest <- mdsout <- switch(mdsmethod, classical = cmdscale(dm, 
        k = mdsdim), kruskal = isoMDS(dm, k = mdsdim), sammon = sammon(dm, 
        k = mdsdim))
    if (mdsmethod == "classical") 
        mds <- mdsout
    else mds <- mdsout$points
    permchange = FALSE
    operm <- NULL
    if (permutations > 0) {
        if (mdsmethod != "classical") {
            for (i in 1:permutations) {
                inumbers <- sample(1:prabobj$n.species, prabobj$n.species)
                dmperm <- dm[inumbers, inumbers]
                mdsout <- switch(mdsmethod, kruskal = isoMDS(dmperm, 
                  k = mdsdim), sammon = sammon(dm, k = mdsdim))
                if (mdsout$stress < mdsbest$stress) {
                  mdsbest <- mdsout
                  operm <- inumbers
                  permchange <- TRUE
                }
            }
            mdsout <- mdsbest
            if (!is.null(operm)) {
                mdsout$points[operm, ] <- mdsbest$points
                mds <- mdsout$points
            }
            operm <- NULL
        }
    }
    if (nnk == 0) {
        if (0 == nclus[1]) 
            nclus <- nclus[-1]
        if (identical(modelid,"all")) 
            kem <- mclustBIC(mds, G = nclus, modelNames = c("EII", 
                "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "EEV", 
                "VEV", "VVV"))
        else {
            if (identical(modelid,"noVVV")) 
                kem <- mclustBIC(mds, G = nclus, modelNames = c("EII", 
                  "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "EEV", 
                  "VEV"))
            else kem <- mclustBIC(mds, G = nclus, modelNames = modelid)
        }
    }
    else {
        kn <- NNclean(mds, k = nnk)
        if (identical(modelid,"all")) 
            kem <- mclustBIC(mds, G = nclus, modelNames = c("EII", 
                "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "EEV", 
                "VEV", "VVV"), initialization = list(noise = as.logical(1 - 
                kn$z)))
        else {
            if (identical(modelid,"noVVV")) 
                kem <- mclustBIC(mds, G = nclus, modelNames = c("EII", 
                  "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "EEV", 
                  "VEV"), initialization = list(noise = as.logical(1 - 
                  kn$z)))
            else kem <- mclustBIC(mds, G = nclus, modelNames = modelid, 
                initialization = list(noise = as.logical(1 - 
                  kn$z)))
        }
    }
    skembest <- skem <- summary(kem, mds)
    if (permutations > 0) {
        for (i in 1:permutations) {
            inumbers <- sample(1:prabobj$n.species, prabobj$n.species)
            mdsperm <- mds[inumbers, ]
            if (nnk == 0) {
                if (0 == nclus[1]) 
                  nclus <- nclus[-1]
                if (identical(modelid,"all")) 
                  kem <- mclustBIC(mdsperm, G = nclus, modelNames = c("EII", 
                    "VII", "EEI", "VEI", "EVI", "VVI", "EEE", 
                    "EEV", "VEV", "VVV"))
                else {
                  if (identical(modelid,"noVVV")) 
                    kem <- mclustBIC(mdsperm, G = nclus, modelNames = c("EII", 
                      "VII", "EEI", "VEI", "EVI", "VVI", "EEE", 
                      "EEV", "VEV"))
                  else kem <- mclustBIC(mdsperm, G = nclus, modelNames = modelid)
                }
            }
            else {
                kn <- NNclean(mdsperm, k = nnk)
                if (identical(modelid,"all")) 
                  kem <- mclustBIC(mdsperm, G = nclus, modelNames = c("EII", 
                    "VII", "EEI", "VEI", "EVI", "VVI", "EEE", 
                    "EEV", "VEV", "VVV"), initialization = list(noise = as.logical(1 - 
                    kn$z)))
                else {
                  if (identical(modelid,"noVVV")) 
                    kem <- mclustBIC(mdsperm, G = nclus, modelNames = c("EII", 
                      "VII", "EEI", "VEI", "EVI", "VVI", "EEE", 
                      "EEV", "VEV"), initialization = list(noise = as.logical(1 - 
                      kn$z)))
                  else kem <- mclustBIC(mdsperm, G = nclus, modelNames = modelid, 
                    initialization = list(noise = as.logical(1 - 
                      kn$z)))
                }
            }
            skem <- summary(kem, mdsperm)
            if (skem$bic > skembest$bic) {
                skembest <- skem
                operm <- inumbers
                permchange <- TRUE
            }
        }
    }
    skem <- skembest
    if (!is.null(operm)) {
        skem$classification[operm] <- skembest$classification
        skem$z[operm, ] <- skembest$z
        if (!is.null(attr(skembest, "initialization")$noise)) 
            attr(skem, "initialization")$noise[operm] <- attr(skembest, 
                "initialization")$noise
        skembest <- skem
    }
    mdsr <- mds
    mds[ospecies, ] <- mdsr
    skem$classification[ospecies] <- skembest$classification
    skem$z[ospecies, ] <- skembest$z
    if (!is.null(attr(skembest, "initialization")$noise)) 
        attr(skem, "initialization")$noise[ospecies] <- attr(skembest, 
            "initialization")$noise
    uclustering <- skem$classification
    ncl <- max(uclustering)
    nc <- ncl + 1
    noisec <- min(uclustering) == 0
    clustering <- uclustering
    csum <- function(n, cv) {
        out <- c()
        for (i in 1:length(n)) out[i] <- sum(cv == n[i])
        out
    }
    cs <- csum(1:ncl, clustering)
    ocs <- order(-cs)
    prob <- skem$parameters$pro
    pmean <-  skem$parameters$mean
    psigma <- skem$parameters$variance$sigma
    pcsigma <- skem$parameters$variance$cholsigma
    pscale <- skem$parameters$variance$scale
    pshape <- skem$parameters$variance$shape
    pz <- skem$z
#    print(prob)
#    print(pmean)
#    print(psigma)
#    print(pcsigma)
#    print(pscale)
#    print(pshape)
#    print(skem)
    for (i in 1:ncl){
      clustering[uclustering == ocs[i]] <- i
      skem$classification <-  clustering
      skem$parameters$pro[i] <- prob[ocs[i]]
      skem$parameters$mean[,i] <- pmean[,ocs[i]]
      skem$parameters$variance$sigma[,,i] <- psigma[,,ocs[i]]
      skem$parameters$variance$cholsigma[,,i] <- pcsigma[,,ocs[i]]
      if (length(skem$parameters$variance$scale)>1)
        skem$parameters$variance$scale[i] <- pscale[ocs[i]]
      if (is.matrix(skem$parameters$variance$shape))
        skem$parameters$variance$shape[,i] <- pshape[,ocs[i]]
      skem$z[,i] <- pz[,ocs[i]]
    }
    if (noisec & nc == 1) 
        clsym <- rep("N", prabobj$n.species)
    else {
        if (noisec) {
            symbols <- c("N", sapply(1:ncl, toString))
            clsym <- symbols[clustering + 1]
        }
        else {
            symbols <- sapply(1:ncl, toString)
            clsym <- symbols[clustering]
        }
    }
    for (i in 1:ncl) if (sum(clustering == i) < 2) 
        clsym[clustering == i] <- "N"
    plot(mds, pch = clsym)
    out <- list(clustering = clustering, clustsummary = skem, 
        bicsummary = kem, points = mds, nnk = nnk, mdsdim = mdsdim, 
        mdsmethod = mdsmethod, symbols = clsym, permutations = permutations, 
        permchange = permchange)
    class(out) <- "prabclust"
    out
}

# "prabclust" <- function (prabobj, mdsmethod = "classical", mdsdim = 4,
#                        nnk = ceiling(prabobj$n.species/40), 
#     nclus = 0:9, modelid = "all", permutations=0) 
# {
#     require(MASS)
#     require(mclust)
#     # "Data-alphabetical" ordering 
#     oregions <- order(prabobj$specperreg)
#     prabo1 <- prabobj$prab[oregions,]
#     ospecies <- do.call("order",as.data.frame(t(prabo1)))    
#     dm <- prabobj$distmat[ospecies,ospecies]
#     if (mdsmethod != "classical") {
#         mindm <- min(dm[dm > 0])/10
#         for (i in 1:(prabobj$n.species - 1))
#           for (j in (i + 1):prabobj$n.species) if (dm[i, j] < mindm) 
#             dm[i, j] <- dm[j, i] <- mindm
#     }
#     mdsbest <- mdsout <-
#       switch(mdsmethod, classical = cmdscale(dm, k = mdsdim), 
#         kruskal = isoMDS(dm, k = mdsdim), sammon = sammon(dm, 
#             k = mdsdim))
#     if (mdsmethod=="classical") mds <- mdsout
#     else mds <- mdsout$points
#     permchange=FALSE
#     operm <- NULL
#     if (permutations>0){
#       if (mdsmethod!="classical"){
#         for (i in 1:permutations){
#           inumbers <- sample(1:prabobj$n.species,prabobj$n.species)
#           dmperm <- dm[inumbers,inumbers]
#           mdsout <- switch(mdsmethod,  
#             kruskal = isoMDS(dmperm, k = mdsdim),
#                         sammon = sammon(dm, k = mdsdim))
#           if (mdsout$stress<mdsbest$stress){
#             mdsbest <- mdsout
#             operm <- inumbers
#             permchange <- TRUE
#           } # if mds$stress<beststress
#         } # for i
#         mdsout <- mdsbest
#         if (!is.null(operm)){
#           mdsout$points[operm,] <- mdsbest$points
#           mds <- mdsout$points
#         }
#         operm <- NULL
#       } # if mdsmethod!="classical"
#     } # if permutations>0  
#     if (nnk==0){
#       if (0==nclus[1])
#         nclus <- nclus[-1]
#       if (modelid == "all") 
#         kem <- mclustBIC(mds, G = nclus, modelNames=c("EII", 
#                 "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "EEV", 
#                 "VEV","VVV"))
#       else {
#         if (modelid == "noVVV") 
#             kem <- mclustBIC(mds, G = nclus, modelNames = c("EII", 
#                 "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "EEV", 
#                 "VEV"))
#         else kem <- mclustBIC(mds, G = nclus, modelNames = modelid)
#       }
#     }
#     else{  
#       kn <- NNclean(mds, k = nnk)
#       if (modelid == "all") 
#         kem <- mclustBIC(mds, G = nclus, modelNames=c("EII", 
#                 "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "EEV", 
#                 "VEV","VVV"),
#                          initialization=list(noise = as.logical(1 - kn$z)))
#       else {
#         if (modelid == "noVVV") 
#             kem <- mclustBIC(mds, G = nclus, modelNames = c("EII", 
#                 "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "EEV", 
#                 "VEV"), initialization=list(noise = as.logical(1 - kn$z)))
#         else kem <- mclustBIC(mds, G = nclus, modelNames = modelid, 
#             initialization=list(noise =as.logical( 1 - kn$z)))
#       }
#     }
#     skembest <- skem <- summary(kem, mds)
#     if (permutations>0){
#       for (i in 1:permutations){
#         inumbers <- sample(1:prabobj$n.species,prabobj$n.species)
#         mdsperm <- mds[inumbers,]
#         if (nnk==0){
#           if (0==nclus[1])
#             nclus <- nclus[-1]
#           if (modelid == "all") 
#             kem <- mclustBIC(mdsperm, G = nclus, modelNames=c("EII", 
#                 "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "EEV", 
#                 "VEV","VVV"))
#           else {
#             if (modelid == "noVVV") 
#                 kem <- mclustBIC(mdsperm, G = nclus, modelNames = c("EII", 
#                     "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "EEV", 
#                     "VEV"))
#             else kem <- mclustBIC(mdsperm, G = nclus, modelNames = modelid)
#           } # else (modelid!="all")
#         } # if nnk==0
#         else{  
#           kn <- NNclean(mdsperm, k = nnk)
#           if (modelid == "all") 
#             kem <- mclustBIC(mdsperm, G = nclus, modelNames=c("EII", 
#                     "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "EEV", 
#                     "VEV","VVV"),initialization=
#                              list(noise = as.logical(1 - kn$z)))
#           else {
#             if (modelid == "noVVV") 
#               kem <- mclustBIC(mdsperm, G = nclus, modelNames = c("EII", 
#                     "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "EEV", 
#                     "VEV"), initialization=list(noise = as.logical(1 - kn$z)))
#             else kem <- mclustBIC(mdsperm, G = nclus, modelNames = modelid, 
#                 initialization=list(noise =as.logical(1 - kn$z)))
#           } # else - modelid!="all"
#         } # else - nnk>0
#         skem <- summary(kem, mdsperm)
#         if (skem$bic>skembest$bic){
#           skembest <- skem
#           operm <- inumbers
#           permchange <- TRUE
#         }
#       } # for i
#     } # if permutations>0    
#     skem <- skembest
#     if (!is.null(operm)){
#       skem$classification[operm] <- skembest$classification
#       skem$z[operm,] <- skembest$z
#       if (!is.null(attr(skembest,"initialization")$noise))
#         attr(skem,"initialization")$noise[operm] <-
#           attr(skembest,"initialization")$noise
#       skembest <- skem
#     }
#     mdsr <- mds
#     mds[ospecies,] <- mdsr
#     skem$classification[ospecies] <- skembest$classification
#     skem$z[ospecies,] <- skembest$z
#     if (!is.null(attr(skembest,"initialization")$noise))
#       attr(skem,"initialization")$noise[ospecies] <-
#         attr(skembest,"initialization")$noise
#     uclustering <- skem$classification
#     ncl <- max(uclustering)
#     nc <- ncl+1
#     noisec <- min(uclustering)==0
#     clustering <- uclustering
#     csum <- function(n, cv) {
#         out <- c()
#         for (i in 1:length(n)) out[i] <- sum(cv == n[i])
#         out
#     }
# #    print(skem$z)
# #    str(skem)
# #    print(ncl)
# #    print(nc)
#     cs <- csum(1:ncl, clustering)
#     ocs <- order(-cs)
#     for (i in 1:ncl) clustering[uclustering == ocs[i]] <- i
#     clusyms <- sapply(1:9,toString)
#     if (ncl>9)
#       clusyms <- c(clusyms, intToUtf8(97:(97+ncl-10),multiple=TRUE)) 
#     if (noisec & nc == 1)
#       clsym <- rep("N",prabobj$n.species)
#     else {
#       if (noisec){ 
#          symbols <- c("N",clusyms)
#          clsym <- symbols[clustering+1]
#        }
#        else{
#          symbols <- clusyms
#          clsym <- symbols[clustering]
#        }
#     }
#     for (i in 1:ncl) if (sum(clustering == i) < 2) 
#         clsym[clustering == i] <- "N"
#     plot(mds, pch = clsym)
#     out <- list(clustering = clustering, clustsummary = skem, 
#         bicsummary = kem, points = mds, nnk = nnk, mdsdim = mdsdim, 
#         mdsmethod = mdsmethod, symbols = clsym, permutations=permutations,
#                 permchange=permchange, csreorder=ocs)
#     class(out) <- "prabclust"
#     out
# }

#   \item{csreorder}{integer vector. This gives the numbering of the
#     components in \code{clustsummary} relative to
#     \code{clustering}. Usually, \code{clustering} and \code{symbols}
#     will be used, but in order to use the information in
#     \code{clustsummary} (parameter values, posterior assignment
#     probabilities etc.), it has to be taken into account that cluster
#     no. 1 in \code{clustering} corresponds to cluster
#     no. \code{csreorder[1]} in \code{clustsummary} and so on. Noise, if
#     present, is numbered 0 in \code{clustering} as well as
#     \code{clustsummary}.}
