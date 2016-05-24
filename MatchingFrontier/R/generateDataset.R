generateDataset <-
function(frontier.object, N){
    X <- nrow(frontier.object$dataset) - frontier.object$frontier$Xs
    ind <- which(abs(X - N) == min(abs(X - N)))[1]
    this.dat.inds <- unlist(frontier.object$frontier$drop.order[ind:length(frontier.object$frontier$drop.order)])
    dataset <- frontier.object$dataset[this.dat.inds,]

    # Add weights
    if(frontier.object$ratio == 'variable'){    
        w <- makeWeights(dataset, frontier.object$treatment)
        dataset$weights <- w
    }
    
    return(dataset)
}
