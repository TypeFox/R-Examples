columnmean <-
function(x){
    colvec <- colMeans(x)
    nrow <- dim(x)[1]
    ones <- rep(1,nrow)
    ones %*% t(colvec)
 }
