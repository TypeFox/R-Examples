pamr.cv <-
function(fit, data, nfold = NULL, folds = NULL ,...)
{
        x <- data$x[fit$gene.subset, fit$sample.subset]

        if( !is.null(data$y) & !is.null(data$proby)){
           stop("Must have exactly one of y and  proby  present in the data object")
         }
        
        y <- NULL
        proby <- NULL
        
        if(!is.null(fit$y)){
           y<-  factor(fit$y[fit$sample.subset])
         }
        
        if(!is.null(fit$proby)){
           proby<-  fit$proby[fit$sample.subset,]
         }
        
        this.call <- match.call()
        
# three possibilities, 
# problem.type= class: y are class labels, proby=NULL
#               surv.km: y=NULL, proby are soft class probabilities from KM
#               surv.latent: y are latent class labels,
#                   proby are soft class probabilities from KM
# note; problem type is in fit$problem.type
        
        junk <- nsccv(x, y=y, proby=proby, object = fit, nfold=nfold, folds=folds, 
survival.time=data$survival.time, censoring.status = data$censoring.status, 
ngroup.survival=fit$ngroup.survival, problem.type=fit$problem.type, ...)

        junk$call <- this.call
        
        junk$sample.subset <- fit$sample.subset
        class(junk)="pamrcved"
        junk
}

