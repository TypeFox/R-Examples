estimateEffects <-
function(frontier.object, formula, prop.estimated = 1, mod.dependence.formula, continuous.vars = NA, seed = 1,
         means.as.cutpoints = FALSE){
    
    set.seed(seed)
    
    # These are the points that we'll estimate
    point.inds <- sort(sample(1:length(frontier.object$frontier$Xs),
                              round(length(frontier.object$frontier$Xs) * prop.estimated)))

    coefs <- vector(mode="list", length = length(point.inds))
    CIs <- vector(mode="list", length= length(point.inds))
    mod.dependence <- vector(mode="list", length= length(point.inds))

    treatment <- frontier.object$treatment

    if(!is.na(continuous.vars[1])){
        if(means.as.cutpoints){
            cutpoints <- lapply(continuous.vars, function(x) mean(frontier.object$dataset[[x]]))
            names(cutpoints) <- continuous.vars
        }
        cutpoints <- getCutpointList(frontier.object$dataset, mod.dependence.formula, continuous.vars)
    }

    covs <- strsplit(as.character(mod.dependence.formula[3]), '\\+')
    covs <- unlist(lapply(covs, trim))
    covs <- covs[!(covs %in% treatment)]
    
    print(cutpoints)
    pb <- txtProgressBar(min = 1, max = length(point.inds), style = 3)

    for(i in 1:length(point.inds)){
        this.dat.inds <- unlist(frontier.object$frontier$drop.order[point.inds[i]:length(frontier.object$frontier$drop.order)])
        dataset <- frontier.object$dataset[this.dat.inds,]

        if(frontier.object$ratio == 'variable'){
            w <- makeWeights(dataset, treatment)
            dataset$w <- w            
            results <- lm(formula, dataset, weights = w)
        } else {
            results <- lm(formula, dataset)
        }

        tryCatch( 
            this.mod.dependence <- modelDependence(dataset,
                                                   treatment,
                                                   mod.dependence.formula,
                                                   verbose = FALSE, cutpoints = cutpoints),
            error = function(e) this.mod.dependence <- NA
        )

        if(!is.na(this.mod.dependence[1])){           
            this.sig.hat <- this.mod.dependence

        } else{
            this.sig.hat <- NA
        }
        
        
        coefs[i] <- coef(results)[frontier.object$treatment]
        CIs[[i]] <- confint(results)[frontier.object$treatment,]       
        mod.dependence[i] <- this.sig.hat
        
        setTxtProgressBar(pb, i)
    }
    close(pb)
    
    return(list(Xs = frontier.object$frontier$Xs[point.inds], coefs = unlist(coefs), CIs = CIs, mod.dependence = unlist(mod.dependence)))
}

