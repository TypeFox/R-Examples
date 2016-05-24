featureSelection.basic <- function(cols){
        features <- 0
        if(length(cols)>1){
                features <- diffRows[which(rowMeans(x[diffRows,cols])>fitnessArgs$consistency)]
        }
        features
}

