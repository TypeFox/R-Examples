ent <-
function(y){
    p <- table(y)/length(y)
    e <- -sum(p*mylog(p))
    return(e)
}
