#' @export
print.loq <- function(x, ... ) {
    if (!inherits(x,"loq")) {
        stop("not loq object")
    }
    
    lloq <- unlist(lapply(x, function(y) y$lloq))
    uloq <- unlist(lapply(x, function(y) y$uloq))
    method <- unlist(lapply(x, function(y) y$method))  
    ans <- data.frame(analyte = names(x), lloq=lloq, uloq=uloq, method=method)
    rownames(ans) <- 1:nrow(ans)
    print(ans)
}

