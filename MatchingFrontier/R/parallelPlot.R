parallelPlot <-
function(frontier.object, N, variables, treated.col = 'grey', control.col = 'black'){

    dataset <- generateDataset(frontier.object, N)
    
    col <- rep(NA, nrow(dataset))
    col[dataset[[frontier.object$treatment]] == 1] <- treated.col
    col[dataset[[frontier.object$treatment]] == 0] <- control.col

    dataset <- dataset[,variables]
    parcoord(dataset, col = col)
}
