getConfs <- function(nvertex, ncolor){
    if(nvertex < 2) stop("There are at least 2 vertices!")
    if(ncolor < 2) stop("There are at least two possible choices for each vertex!")

    subconfs <- function(subpartial){
        lapply(1:ncolor, function(x) list(c(subpartial, x)))
    }

    confs <- function(nvertex){
        if(nvertex > 2){
            partial <- matrix(unlist(confs(nvertex-1)), nrow=nvertex-1)
            lapply(1:ncol(partial), function(x) subconfs(partial[,x]))
        }
        else{
            lapply(1:ncolor, function(x) subconfs(x))
        }

    }

    matrix(unlist(confs(nvertex)), nrow=nvertex)
}
