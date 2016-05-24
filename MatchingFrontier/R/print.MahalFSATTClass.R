print.MahalFSATTClass <-
function(x, ...){
    msg <- paste('An imbalance frontier with', as.character(length(x$frontier$Xs)), 'points.\n', sep = ' ')
    cat(msg)
}
