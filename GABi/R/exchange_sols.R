exchangeSols <- function(demes,fitnesses,fittestonly,proximity){

    cat("exchanging solutions\n")
    nsubpops <- length(demes)
    if (fittestonly){
        chrs <- list(rep(NA,nsubpops))
        swapindices <- rep(NA,nsubpops)
        for(i in 1:nsubpops){
            fittest <- sort(fitnesses[[i]],index.return=TRUE,decreasing=TRUE)$ix[1];
            chrs[[i]] <- demes[[i]][fittest,];
            swapindices[i] <- fittest
        }
    }
    else{
        chrs <- list(rep(NA,nsubpops))
        swapindices <- rep(NA,nsubpops)
        for(i in 1:nsubpops){
            swapchoice <- runif(1,min=1,max=dim(demes[[i]])[1])
            chrs[[i]] <- demes[[i]][swapchoice,];
            swapindices[i] <- swapchoice
        }
    }
    if (proximity){
    # enforces exchange of solutions only to adjacent 'islands'
        swapindex <- swapindices[1]
        demes[[1]][swapindex,] <- chrs[[nsubpops]]
        for (i in 2:nsubpops){
            swapindex <- swapindices[i]
            demes[[i]][swapindex,] <- chrs[[i-1]]
        }
    }
    else {
        listrefs <- sample(c(1:nsubpops))
        # indirect for random permutation of list of solutions to be exchanged
        for (i in 1:nsubpops){
            swapindex <- swapindices[i]
            demes[[i]][swapindex,] <- chrs[[listrefs[i]]]
        }
    }
        demes
}

