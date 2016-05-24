mmtfa.parallel <-function (x, Gs, Qs, clas, init, scale, models, 
                           dfstart, dfupdate, gauss, eps, known, numcores)
{
  requireNamespace("parallel", quietly=TRUE)
  origmod <- models
  modgen <- modelgen()
  oldmods <- modgen$modold
  allmods <- modgen$allmodels
  models <- oldmods[match(models, allmods)]
  
  modrep <- rep(models,times=length(Gs)*length(Qs))
  grep <- rep(Gs, each=length(Qs)*length(models))
  qrep<- rep(Qs,times=length(Gs), each=length(models),)
  if(length(grep[grep == 1]) > 0) {
      mod <- NA
      k <- 1
      UUUUgroup <- c("UUUU", "UUUC", "UCUU","UCUC", "CUUU", "CUUC", 
                     "CCUU", "CCUC", "Mt1U","Mt1C", "Mt2U", "Mt2C", 
                     "Mt3U", "Mt3C", "Mt4U", "Mt4C")
      uuudum <- models[models %in% UUUUgroup]
      if(length(uuudum) > 0){
        mod[k] <- uuudum[1]
        k <- k + 1
      }
      UUCUgroup <- c("UUCU", "UUCC", "UCCU","UCCC","CUCU", "CUCC", "CCCU", "CCCC")
      uucdum <- models[models %in% UUCUgroup]
      if (length(uucdum) >0){
        mod[k] <-uucdum[1]
        k <- k + 1
      }
      modrep <- c(rep(mod,times=length(Qs)), modrep[!grep == 1])
      qrep <- c(rep(Qs,each=length(mod)), qrep[!grep == 1])
      grep <- c(rep(1,times=length(mod)*length(Qs)), grep[!grep == 1])  
  }
  modrep <- allmods[match(modrep,oldmods)]
  runvec <- 1:length(modrep)
  #instal
  if(is.null(numcores)){ numcores <- detectCores() }
  clus <- makeCluster(numcores)
#   clusterEvalQ(clus, library(mclust))
#   clusterEvalQ(clus, library(e1071))
  clusterEvalQ(clus, library(parallel))
  clusterEvalQ(clus, library(mvnfast))
  clusterEvalQ(clus, library(matrixStats))
  clusterEvalQ(clus, library(mmtfa))
#  clusterExport(clus, )
  clusterExport(clus, ls(environment()), envir = environment())
  runlist <- clusterApplyLB(clus,runvec, function(g) mmtfaEM(x, grep[g], qrep[g], clas = clas, init = init, 
              scale = scale, models=modrep[g], dfstart = dfstart, gauss = gauss, dfupdate=dfupdate, eps=eps, known=known))
  stopCluster(clus)
  #
  bic.vec <- unlist(lapply(runlist, function(g) g$bic))
  icl.vec <- unlist(lapply(runlist, function(g) g$iclresults$icl))
 #modvec <- unlist(lapply(runlist, function(g) g$bestmod))
  #gvec <- unlist(lapply(runlist, function(g) max(g$bestg))
  #qvec <- unlist(lapply(runlist, function(g) g$bestq))
  modvec <- unlist(lapply(runlist, function(g) g$modelname[1]))
  gvec <- unlist(lapply(runlist, function(g) max(g$G)))
  qvec <- unlist(lapply(runlist, function(g) max(g$Q)))
  ### when maxG is larger than length of G error
  biccube <- iclcube <- array(-Inf,dim=c(length(models),length(Qs),length(Gs)),dimnames=list(origmod,paste("Q=",Qs,sep=""),paste("G=",Gs,sep="")))
  for(i in 1:length(bic.vec)){
    biccube[modvec[i], qvec[i], gvec[i]] <- bic.vec[i]
    iclcube[modvec[i], qvec[i], gvec[i]] <- icl.vec[i]
  }
  final <- runlist[[which.max(bic.vec)]]
  final$iclresults$allicl <- iclcube
 parallel.store <- final
 parallel.store[["allbic"]] <- biccube
 parallel.store
}