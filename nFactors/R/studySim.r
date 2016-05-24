studySim <- function(var, nFactors, pmjc, loadings, unique, N, repsim, reppar,
                     stats=1, quantile=0.5, model="components", r2limen=0.75,
                     all=FALSE, dir=NA, trace=TRUE) {
 nsubjects <- N
 result    <- NULL
 id        <- 0
 nid       <- length(nFactors) * length(loadings) * length(pmjc) * length(var) * length(unique) * length(nsubjects)
 for (i in 1:length(nFactors))  {
  for (j in 1:length(loadings))  {
   for (l in 1:length(pmjc))      {
    for (n in 1:length(var))       {
     for (k in 1:length(unique))    {
      for (m in 1:length(nsubjects)) {
       id    <- id + 1
       kid  <- paste(id,"/",nid,sep="")
       ident <- c(nFactors=nFactors[i], loadings=loadings[j], unique=unique[k], quantile=quantile,
                  pmjc=pmjc[l], nsubjects=nsubjects[m], var=var[n], reppar=reppar,
                  repsim=repsim, id=kid, model=model)
       if (trace == TRUE) print(ident)
       fStruct <- generateStructure(var=var[n], mjc=nFactors[i], pmjc=pmjc[l], loadings=loadings[j], unique=unique[k])
       fSim    <- structureSim(fload=as.matrix(fStruct), reppar=reppar, repsim=repsim, details=FALSE, all=all,
                               N=nsubjects[m], quantile=quantile, model=model, r2limen=r2limen)[[2]][stats,]
       if (length(stats) == 1) {
        fSim  <- data.frame(var=var[n], nsubjects=nsubjects[m], nfactors=nFactors[i], pmjc=pmjc[l],
                            loadings=loadings[j], unique=unique[k], t(fSim), repsim=repsim, reppar=reppar)
        }
       if (length(stats) > 1) {
        ls    <- length(stats)
        info  <-  data.frame(stats   =rownames(fSim),       id       =rep(id, ls),
                             var     =rep(var[n], ls),      nsubjects=rep(nsubjects[m], ls),
                             nfactors=rep(nFactors[i], ls), pmjc     =rep(pmjc[l], ls),
                             loadings=rep(loadings[j], ls), unique   =rep(unique[k], ls),
                             repsim  =rep(repsim, ls),      reppar   =rep(reppar, ls))
         fSim <- data.frame(info, fSim)
         }
        result           <- rbind(result, fSim)
        rownames(result) <- 1:dim(result)[1]
        fString          <- paste("RES_", paste(ident,"_", sep="", collapse=""), sep="")
        # if (!is.na(dir)) save("fSim", file=paste(dirPack, fString,".Rdata", sep=""))  # Old erroneous code
        if (!is.na(dir)) save("fSim", file=paste(dir, fString,".Rdata", sep=""))
  }}}}}}
 return(result)
 }
 
 