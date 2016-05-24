fps <- function(population,fitnesses,elitism){

        popsize <- nrow(population)
        goodSols <- which(fitnesses>0)
        probShare <- fitnesses[goodSols]/sum(fitnesses[goodSols])
        cutoffs <- cumsum(probShare)
        intpop <- array(dim=dim(population))
        if (elitism){
                fittest <- sort(fitnesses,index.return=TRUE,decreasing=TRUE)$ix[1]
                intpop[1,] <- population[fittest,]
                selectionPoints <- runif(nrow(intpop)-1)
                selectedSols <- goodSols[unlist(lapply(selectionPoints,function(x,cutoffs)min(which(cutoffs>x)),cutoffs=cutoffs))]
                intpop[c(2:popsize),] <- population[selectedSols,]
        }
  
        else {
                selectionPoints <- runif(nrow(intpop))
                selectedSols <- goodSols[unlist(lapply(selectionPoints,function(x,cutoffs)min(which(cutoffs>x)),cutoffs=cutoffs))]
                intpop <- population[selectedSols,]
        }
        intpop
        
}
