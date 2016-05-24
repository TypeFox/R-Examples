pamr.train <-
function(data, gene.subset=NULL, sample.subset=NULL,
         threshold = NULL, n.threshold = 30,
        scale.sd = TRUE, threshold.scale = NULL, se.scale = NULL, offset.percent = 50, hetero=NULL,
         prior = NULL,  remove.zeros = TRUE, sign.contrast="both", ngroup.survival=2)

{
  
# modified aug 2003 to add survival analysis facilities
#
# possibilities for different outcomes in data object "data":
#         y= class variable => classification problem
#         proby= matrix of class probabilities => "soft classification"
#              (not currently used by Excel interface)
#        survival time, censoring status present => survival analysis problem
#
# here is how the two  problem types  are passed to nsc:
#
#       class: y is class label, proby, prob.ytest not used
#       surv.km: proby are soft class labels, computed from kaplan-meier
  
        this.call <- match.call()

 if(!is.null(data$y)){problem.type <- "class"}
  if(!is.null(data$survival.time)) {problem.type <- "surv.km"}

        
          if(!is.null(data$proby) & !is.null(data$y)) {
           stop("Can't have both proby and y present in data object")
         }
        
          if(!is.null(data$y) & !is.null(data$survival.time)) {
    stop("Can't have both class label y and survival.time present in data object")
  }
        
        if(!is.null(data$y) & !is.null(data$censoring.status)) {
           stop("Can't have both class label y and censoring status present in data object")
         }
         if(!is.null(data$survival.time) & is.null(data$censoring.status)) {
           stop("Survival time specified but censoring status missing")
         }
          if(is.null(data$survival.time) & !is.null(data$censoring.status)) {
           stop("Censoring status specified but survival time missing")
         }
  
        
        y <- data$y
        proby <- data$proby
        ytest <- NULL
        xtest <- NULL
        prob.ytest <- NULL
       
        
        cutoffs.survival <- NULL
        
#       estimate class probabilities via Kaplan-Meier
#        use Cox score cutoff of 2.0 or 20th largest score, whichever is smaller
#       this ensures at least 20 genes are used for clustering
        
        if(!is.null(data$survival.time)){
          junk <- pamr.surv.to.class2(data$survival.time, data$censoring.status,
             n.class=ngroup.survival)
                 
            proby <- junk$prob
            cutoffs.survival <- junk$cutoffs
        }

        
# ytemp is just used for computation of the prior

 if(!is.null(y)){ ytemp <- y}
 if(is.null(y) & !is.null(proby)){ ytemp <- apply(proby,1, which.is.max)}
 if(is.null(sample.subset)){sample.subset <-1:ncol(data$x)}
 if(is.null(gene.subset)){gene.subset <-1:nrow(data$x)}

# for survival analysis, make default prior the equal prior

  if(is.null(prior) & !is.null(data$survival.time) ){
         prior <- rep(1/ngroup.survival, ngroup.survival)
}

        if(is.null(prior) & is.null(data$survival.time) )
          {prior <- table(ytemp[sample.subset])/length(ytemp[sample.subset])
           prior <- prior[prior!=0]
        }
      
    if(!is.null(y)){y <-  factor(y[sample.subset])}
        
    if(!is.null(proby)){
        proby <- proby[sample.subset,]}
        junk <- nsc(data$x[gene.subset, sample.subset], y=y, proby=proby,
 xtest=xtest, ytest=ytest, prob.ytest=prob.ytest,
          offset.percent=offset.percent,  threshold = threshold, hetero=hetero,
        n.threshold = n.threshold,  scale.sd= scale.sd,
                    threshold.scale=threshold.scale,
           se.scale= se.scale, prior=prior, remove.zeros=remove.zeros,
            sign.contrast=sign.contrast, problem.type=problem.type)

        junk$call <- this.call
        junk$gene.subset <- gene.subset
        junk$sample.subset <- sample.subset
        junk$survival.time <- data$survival.time
        junk$censoring.status <- data$censoring.status
       junk$cutoffs.survival <- cutoffs.survival
        junk$ngroup.survival <- ngroup.survival
        junk$problem.type <- problem.type
        class(junk)="pamrtrained"
        junk
}
