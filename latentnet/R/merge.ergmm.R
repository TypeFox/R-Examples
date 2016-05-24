## A routine to merge replications of ergmm calls with same inputs
## but different random seeds and/or starting values.

merge.ergmm<-function(x,y,...,verbose=FALSE){
  # Pairwise is may not be the most efficient way to do it, but it is the simplest.
  if(verbose) cat(".")
  object<-combine.2ergmm(x,y)
  for(mergewith in list(...)){
    if(verbose) cat(".")
    object<-combine.2ergmm(object,mergewith)
  }
  if(verbose) cat("\n")
  object<-statsreeval.ergmm(object,rerun=TRUE)
  object
}

combine.2ergmm<-function(fit1,fit2){
  if(!isTRUE(all.equal(fit1[["prior"]],fit2[["prior"]]))) stop("The prior distributions for the runs were not the same.")
  # Networks are never identical, so...
  net1<-fit1[["model"]][["Yg"]]
  fit1[["model"]][["Yg"]]<-NULL
  net2<-fit2[["model"]][["Yg"]]
  fit2[["model"]][["Yg"]]<-NULL
  if(!isTRUE(all.equal(fit1[["model"]],fit2[["model"]]))) stop("The models for the runs were not the same.")
  fit1[["model"]][["Yg"]]<-net1
  fit2[["model"]][["Yg"]]<-net2
  
  ## OK, the two fits are with the same inputs.
  ## Now, begin merging.
  ## [["model"]] and [["prior"]] are identical.

  # MCMC MLE
  if(fit2[["mcmc.mle"]][["lpY"]]>fit1[["mcmc.mle"]][["lpY"]])
    fit1[["mcmc.mle"]]<-fit2[["mcmc.mle"]]

  # MCMC posterior mode
  if(lpsum(fit2[["mcmc.pmode"]])>
     lpsum(fit1[["mcmc.pmode"]]))
    fit1[["mcmc.pmode"]]<-fit2[["mcmc.pmode"]]

  # Burnin start
  if(!("burnin.starts" %in% names(fit1)))
    fit1[["burnin.starts"]]<-list(fit1[["burnin.start"]])
  if(!("burnin.starts" %in% names(fit2)))
    fit2[["burnin.starts"]]<-list(fit2[["burnin.start"]])
  fit1[["burnin.starts"]]<-c(fit1[["burnin.starts"]],fit2[["burnin.starts"]])

  # Start
  if(!("starts" %in% names(fit1)))
    fit1[["starts"]]<-list(fit1[["start"]])
  if(!("starts" %in% names(fit2)))
    fit2[["starts"]]<-list(fit2[["start"]])

  # Sampling start (already a list)
  fit1[["sampling.start"]]<-c(fit1[["sampling.start"]],fit2[["sampling.start"]])

  # MCMC sample itself
  fit1[["sample"]]<-.stack.ergmm.par.list.list(c(unstack.ergmm.par.list(fit1[["sample"]]),unstack.ergmm.par.list(fit2[["sample"]])))

  # Control
  if(!("controls" %in% names(fit1)))
    fit1[["controls"]]<-list(fit1[["control"]])
  if(!("controls" %in% names(fit2)))
    fit2[["controls"]]<-list(fit2[["control"]])
  fit1[["controls"]]<-c(fit1[["controls"]],fit2[["controls"]])

  fit1[["control"]][["sample.size"]]<-sum(sapply(seq(along=fit1[["controls"]]),function(i)fit1[["controls"]][[i]][["sample.size"]]))
  fit1[["control"]][["interval"]]<-mean(sapply(seq(along=fit1[["controls"]]),function(i)fit1[["controls"]][[i]][["interval"]]))
  fit1[["control"]][["burnin"]]<-mean(sapply(seq(along=fit1[["controls"]]),function(i)fit1[["controls"]][[i]][["burnin"]]))
  fit1[["control"]][["threads"]]<-sum(sapply(seq(along=fit1[["controls"]]),function(i)fit1[["controls"]][[i]][["threads"]]))

  fit1
}
