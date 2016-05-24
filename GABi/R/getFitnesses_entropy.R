getFitnesses.entropy <- function(chr){
        score <- 0
        if(sum(chr)>1){
                features <- which(rowMeans(x[diffRows,chr==1])>fitnessArgs$consistency)
                score <- -length(setdiff(which(chr==1),fitnessArgs$tabu$samples))*sum(log(fitnessArgs$featureWeights[diffRows[setdiff(features,fitnessArgs$tabu$features)]]))
        }
        score
}

