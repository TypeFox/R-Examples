dist2list <-
function(dist){
    if(!class(dist) == "dist"){
        stop("the input data must be a dist object.")
    }
    dat <- as.data.frame(as.matrix(dist))
    if(is.null(names(dat))){
        rownames(dat) <- paste(1:nrow(dat))
    }
    value <- stack(dat)$values
    rnames <- rownames(dat)
    namecol <- expand.grid(rnames,rnames)
    colnames(namecol) <- c("col", "row")
    res <- data.frame(namecol, value)
    return(res)
}

