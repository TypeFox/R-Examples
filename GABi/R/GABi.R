if(getRversion() >= "2.15.1")  utils::globalVariables(c("x","fitnessArgs","diffRows"))

GABi <- function(x,nSols=0,convergenceGens=40,popsize=256,mfreq=1,xfreq=0.5,maxNgens=200,keepBest=FALSE,identityThreshold=0.75,nsubpops=4,experiod=10,diffThreshold=0.9,verbose=FALSE,maxLoop=1,fitnessArgs=list(consistency=0.8,featureWeights=rowMeans(x,na.rm=TRUE)),fitnessFun=getFitnesses.entropy,featureSelFun=featureSelection.basic){

        fitnessArgs$tabu <- list(features=0,samples=0)

        row.mins <- apply(x,MARGIN=1,min,na.rm=TRUE)
        row.maxes <- apply(x,MARGIN=1,max,na.rm=TRUE)

        invariant.rows <- which(mapply(function(x,y)x==y,x=row.mins,y=row.maxes))

        diffRows <- setdiff(1:nrow(x),invariant.rows)

        low.rows <- diffRows[which(((rowMeans(x[diffRows,],na.rm=TRUE)-row.mins[diffRows])/(row.maxes[diffRows]-row.mins[diffRows]))<(1- diffThreshold))]
        hi.rows <- diffRows[which(((rowMeans(x[diffRows,],na.rm=TRUE)-row.mins[diffRows])/(row.maxes[diffRows]-row.mins[diffRows]))> diffThreshold)]

        cat(paste(length(low.rows),"constantly low features and",length(hi.rows),"constantly high\n"))
        diffRows <- setdiff(1:nrow(x),c(low.rows,hi.rows,invariant.rows))
        cat(paste(length(diffRows),"sufficiently variable features of",nrow(x),"\n"))
        if(length(diffRows)<2){
                warning("no sufficiently variable features in dataset")
                return(list())
        }

        subpopsize <- ceiling( popsize/ nsubpops)
        demes <- list(rep(NA, nsubpops))

        # define feature selection and fitness function in current environment (so internal variables passed appropriately)
        featureSelFun <- featureSelFun
        fitnessFun <- fitnessFun
        environment(featureSelFun) <- environment()
        environment(fitnessFun) <- environment()

        allSols <- list()
        loopCount <- 0

        while(!length(allSols)>nSols & loopCount<maxLoop){

        cat(paste("initialising population for GA loop",loopCount+1,"\n"))
        for(j in 1:nsubpops){
                demes[[j]] <- array(0,dim=c(subpopsize,ncol(x)))
                for(i in 1:subpopsize){
                        demes[[j]][i,sample(ncol(x),2)]<-1
                }
        }

        bestSols <- list()

        genNo <- 0
        keptSols <- list()
        for(j in 1:nsubpops){
                bestSols[[j]] <- rep(0,ncol(demes[[j]]))
        }
        convergenceCounter <- 0
        notConverged <- TRUE
        transfergen <- experiod

        fitnesshash <- hash()

        cat("evolving GA solutions\n")

        while(notConverged & (genNo<maxNgens)){
                if(verbose) cat(paste("trying generation no.",genNo,"\n"))
                genNo <- genNo+1
                all_fitnesses <- list()
                stagnancy <- rep(FALSE,nsubpops)
                for(j in 1:nsubpops){
                        if(verbose) cat("looking up previously-calculated fitnesses\n")
                        fitnesses <- rep(0,nrow(demes[[j]]))
                        chrkeys <- apply(demes[[j]],MARGIN=1,function(x)paste(as.character(x),collapse=""))
                        hfitnesses <- lapply(chrkeys,function(x)fitnesshash[[x]])
                        prevSeen <- which(!unlist(lapply(hfitnesses,is.null)))
                        if(verbose) cat(paste(length(prevSeen),"solutions previously evaluated\n"))
                        fitnesses[prevSeen] <- unlist(hfitnesses)
                        newSols <- setdiff(1:nrow(demes[[j]]),prevSeen)
                        if(length(newSols)>1){
                                if(verbose) cat(paste("evaluating new fitnesses for",length(newSols),"sols\n"))
                                fitnesses[newSols] <- apply(demes[[j]][newSols,],MARGIN=1,fitnessFun)
                                fitnesses[is.infinite(fitnesses)] <- 0
                                fitnesses[is.na(fitnesses)] <- 0
                                for(k in 1:length(newSols)){
                                        fitnesshash[[chrkeys[newSols[k]]]] <- fitnesses[[newSols[k]]]
                                }
                        }
                        if(length(newSols)==1){
                                if(verbose) cat(paste("evaluating new fitnesses for 1 new sol\n"))
                                fitnesses[newSols] <- fitnessFun(demes[[j]][newSols,])
                                fitnesses[is.infinite(fitnesses)] <- 0
                                fitnesses[is.na(fitnesses)] <- 0
                                fitnesshash[[chrkeys[newSols]]] <- fitnesses[[newSols]]
                        }
                        all_fitnesses[[j]] <- fitnesses
                        if(verbose) cat("fitnesses evaluated\n")
                        new_bestSol <- demes[[j]][which(fitnesses==max(fitnesses))[1],]
                        identical <- as.numeric(unlist(lapply(bestSols,function(x){identical(x,new_bestSol)})))
                        if(verbose) {
                                cat("similarity of bestsol to previous bestsols:")
                                cat(identical)
                                cat("\n")
                        }
                        if(sum(identical)>0){
                                stagnancy[j] <- TRUE
                        }
                        else{
                                convergenceCounter <- 0
                        }
                        bestSols[[j]] <- new_bestSol
                        if(convergenceCounter < convergenceGens){
                                retrycount <- 0
                                while(sum(fitnesses)==0){
                                        retrycount <- retrycount+1
                                        if(verbose) cat("re-initialising population\n")
                                        demes[[j]] <- array(0,dim=c(subpopsize,ncol(x)))
                                        for(i in 1:subpopsize){
                                                demes[[j]][i,sample(ncol(x),2)]==1
                                        }
                                        if(verbose) cat("looking up previously-calculated fitnesses\n")
                                        fitnesses <- rep(0,nrow(demes[[j]]))
                                        chrkeys <- apply(demes[[j]],MARGIN=1,function(x)paste(as.character(x),collapse=""))
                                        hfitnesses <- lapply(chrkeys,function(x)fitnesshash[[x]])
                                        prevSeen <- which(!unlist(lapply(hfitnesses,is.null)))
                                        if(verbose) cat(paste(length(prevSeen),"solutions previously evaluated\n"))
                                        fitnesses[prevSeen] <- unlist(hfitnesses)
                                        newSols <- setdiff(1:nrow(demes[[j]]),prevSeen)
                                        if(length(newSols)>1){
                                                if(verbose) cat(paste("evaluating new fitnesses for",length(newSols),"sols\n"))
                                                fitnesses[newSols] <- apply(demes[[j]][newSols,],MARGIN=1,fitnessFun)
                                                fitnesses[is.infinite(fitnesses)] <- 0
                                                fitnesses[is.na(fitnesses)] <- 0
                                                for(k in 1:length(newSols)){
                                                        fitnesshash[[chrkeys[newSols[k]]]] <- fitnesses[[newSols[k]]]
                                                }
                                        }
                                        if(length(newSols)==1){
                                                if(verbose) cat(paste("evaluating new fitnesses for 1 new sol\n"))
                                                fitnesses[newSols] <- fitnessFun(demes[[j]][newSols,])
                                                fitnesses[is.infinite(fitnesses)] <- 0
                                                fitnesses[is.na(fitnesses)] <- 0
                                                fitnesshash[[chrkeys[newSols]]] <- fitnesses[[newSols]]
                                        }
                                        all_fitnesses[[j]] <- fitnesses
                                        if(retrycount > 100){
                                                diffThreshold <- diffThreshold+0.1
                                                if(diffThreshold>1){
                                                        warning("cannot find biclusters")
                                                        return(list())
                                                }
                                                if(verbose) cat(paste("no biclusters, lowering diffThreshold to",diffThreshold,"\n"))
                                                diffRows <- setdiff(1:nrow(x),invariant.rows)

                                                low.rows <- diffRows[which(((rowMeans(x[diffRows,],na.rm=TRUE)-row.mins[diffRows])/(row.maxes[diffRows]-row.mins[diffRows]))<(1-diffThreshold))]
                                                hi.rows <- diffRows[which(((rowMeans(x[diffRows,],na.rm=TRUE)-row.mins[diffRows])/(row.maxes[diffRows]-row.mins[diffRows]))>diffThreshold)]
                                                if(verbose) cat(paste(length(low.rows),"constantly low features and",length(hi.rows),"constantly high\n"))
                                                diffRows <- setdiff(1:nrow(x),c(low.rows,hi.rows,invariant.rows))
                                                if(verbose) cat(paste(length(diffRows),"differentially-expressed features of",nrow(x),"\n"))
                                                retrycount <- 0
                                        }
                                }
                                if(verbose){
                                        cat("fitnesses:\n")
                                        cat(fitnesses)
                                        cat("\n")
                                }
                                f <- fitnesses/mean(fitnesses[fitnesses>0])
                                if(keepBest & sum(fitnesses)>0){
                                        if(verbose) cat("saving best solutions\n")
                                        goodsols <- which(fitnesses>0)
                                        for(i in 1:length(goodsols)){
                                                identity <- lapply(keptSols,function(x,comp){sum(colMeans(rbind(x$chr,comp))>0.5)/min(c(sum(x$chr),sum(comp)))},comp=demes[[j]][goodsols[i],])
                                                toRemove <- which(identity>identityThreshold)
                                                newSols <- list()
                                                if(length(toRemove)>0){
                                                        newSols <- keptSols[toRemove]
                                                        keptSols <- keptSols[-toRemove]
                                                }
                                                newSols[[length(newSols)+1]] <- list(fitness=fitnesses[goodsols[i]],chr=demes[[j]][goodsols[i],])
                                                removeScores <- unlist(lapply(newSols,function(x){x$fitness}))
                                                ReplaceSol <- which(removeScores==max(removeScores))[1]
                                                keptSols[[length(keptSols)+1]] <- newSols[[ReplaceSol]]
                                        }
                                }
                                if(verbose) cat("performing selection\n")
                                intpop <- fps(population=demes[[j]],fitnesses=f,elitism=TRUE)
                                if(verbose) cat("performing reproduction\n")
                                demes[[j]] <- reproduction(intpop,xfreq,mfreq/ncol(intpop),xoverpoints=1,pinvert=0,elitism=TRUE)
                        }
                        else{
                                cat(paste("population convergence after",genNo,"generations\n"))
                                notConverged <- FALSE
                        }
                }
                if(sum(as.numeric(stagnancy))==nsubpops){
                        convergenceCounter <- convergenceCounter + 1
                }
                if(genNo==transfergen) {
                        demes <- exchangeSols(demes,all_fitnesses,TRUE,FALSE)
                        transfergen <- transfergen + experiod
                }
        }

        population <- demes[[1]]
        if(nsubpops>1){
                for(j in 2:nsubpops){
                        population <- rbind(population,demes[[j]])
                }
        }
        cat("returning final output\n")
        fitnesses <- apply(population,MARGIN=1,fitnessFun)
        fitnesses[is.infinite(fitnesses)] <- 0
        fitnesses[is.na(fitnesses)] <- 0
        goodsols <- which(fitnesses>0)
        if(length(goodsols)>0){
                if(!keepBest){
                        keptSols <- list()
                        keptSols[[1]] <- list(chr=population[goodsols[1],],fitness=fitnesses[1])
                }
                for(i in 1:length(goodsols)){
                        identity <- lapply(keptSols,function(x,comp){sum(colMeans(rbind(x$chr,comp))>0.5)/min(c(sum(x$chr),sum(comp)))},comp=population[goodsols[i],])
                        toRemove <- which(identity>identityThreshold)
                        newSols <- list()
                        if(length(toRemove)>0){
                                newSols <- keptSols[toRemove]
                                keptSols <- keptSols[-toRemove]
                        }
                        newSols[[length(newSols)+1]] <- list(fitness=fitnesses[goodsols[i]],chr=population[goodsols[i],])
                        removeScores <- unlist(lapply(newSols,function(x){x$fitness}))
                        ReplaceSol <- which(removeScores==max(removeScores))[1]
                        keptSols[[length(keptSols)+1]] <- newSols[[ReplaceSol]]
                }
        }
        if(length(keptSols)>0){
                keptSols <- lapply(keptSols,function(x){list(samples=which(x$chr==1),features=featureSelFun(which(x$chr==1)),score=fitnessFun(x$chr))})
        }
        allSols <- c(allSols,keptSols)
        cat(paste(length(allSols),"biclusters discovered\n"))

        newTabuFeatures <- unlist(lapply(keptSols,function(x)x$features))
        newTabuSamples <- unlist(lapply(keptSols,function(x)x$samples))

        fullTabuFeatures <- union(fitnessArgs$tabu$features,newTabuFeatures)
        fullTabuSamples <- union(fitnessArgs$tabu$samples,newTabuSamples)

        fitnessArgs$tabu$features <- fullTabuFeatures
        fitnessArgs$tabu$samples <- fullTabuSamples
        cat("updating arguments to re-run biclustering loop\n")
        environment(featureSelFun) <- environment()
        environment(fitnessFun) <- environment()

        loopCount <- loopCount + 1

        }

        allSols

}

