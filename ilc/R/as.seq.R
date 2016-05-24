as.seq <-
function(op){
    if (!is.logical(op)) op <- bool(op)
    seq(length(op))[op]
}
