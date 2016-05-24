opt.nested.crossval <-
    function(outerfold=10,nprocessors=1,cl=NULL,...){
        extra.vars <- list(...)
        getfolds <-
            function(N,nfolds){
                evenly <- nfolds * N%/%nfolds
                extras <- N - evenly
                groupsize <- evenly/nfolds
                evenly.folds <- rep(1:nfolds,groupsize)
                evenly.folds <- sample(evenly.folds,size=length(evenly.folds))
                extra.folds <- sample(1:nfolds,size=extras)
                folds <- c(evenly.folds,extra.folds)
                return(folds)
            }
        folds <- getfolds(nrow(extra.vars$penalized),nfolds=outerfold)
        clusterIsSet <- "cluster" %in% class(cl)
        if(nprocessors>1 | clusterIsSet){
            if(!clusterIsSet){
                nprocessors <- as.integer(round(nprocessors))
                cl <- makeCluster(nprocessors, type="PSOCK")
            }
            clusterSetRNGStream(cl, iseed=NULL)
            ##do the nested cross-validation
            output.all <- try(parLapply(cl,unique(folds),function(iFold,extra.vars){
                extra.vars$testset <- which(folds==iFold)
                output <- do.call(opt.splitval,args=extra.vars)
                return(output)
            },extra.vars=extra.vars),silent=TRUE)
            if(!clusterIsSet){
                stopCluster(cl)
            }
        }else{
            ##one processor
            ##do the nested cross-validation
            output.all <- try(lapply(unique(folds),function(iFold,extra.vars){
                extra.vars$testset <- which(folds==iFold)
                output <- do.call(opt.splitval,args=extra.vars)
                return(output)
            },extra.vars=extra.vars),silent=FALSE)
        }
        if(class(output.all)!="try-error"){
            output.all <- unlist(output.all)
            output.all <- output.all[match(rownames(extra.vars$penalized),names(output.all))]
        }else{
            output.all <- rep(NA,nrow(extra.vars$penalized))
            names(output.all) <- rownames(extra.vars$penalized)
        }
        return(output.all)
    }



