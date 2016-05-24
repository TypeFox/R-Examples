#' @export
getnext.comb <- function(I, d=1L, drop=TRUE){
    if (d * I$r > .Machine$integer.max) stop("d is too large.")
    if (I$status > 0){
        I$status <- -1L
        return(NULL)
    }
    if (I$replace){
        if (I$status == -1L) {
            I$currInd <- rep(0L, I$r)
        }
        C <- next_combinations_replace(I, d)
    }else{
        if (I$status == -1L) {
            if (I$is.multiset){
                # add 0L to blame lazy evaluation
                I$currInd <- I$multiset[1:I$r] + 0L
            }else{
                I$currInd <- (1:I$r) - 1L
            }
        }
        if (I$is.multiset){
            C <- next_multiset_combinations(I, d)
        }else{
            C <- next_combinations(I, d)
        }
    }
    if (is.null(C)){
        I$status <- -1L
        return(NULL)
    }else if (I$status > 0){
        C <- C[1:I$status, , drop = FALSE]
        d <- I$status
    }
    if (is.null(I$labels)){
        if (drop || d > 1){
            return(C)
        }else{
            return(matrix(C, nrow = 1))
        }
    }else{
        if (I$r == 1){
            return(matrix(I$labels[C], ncol = 1))
        }else if (drop && d == 1) {
            return(I$labels[C])
        }else if (d == 1){
            return(matrix(I$labels[C], nrow = 1))
        }else{
            return(matrix(I$labels[C], ncol = I$r))
        }
    }
}

#' @export
#' @method getlength comb
getlength.comb <- function(I){
    if (I$replace){
        return(choose(I$unique_n + I$r - 1, I$r))
    }else{
        if (I$is.multiset){
            return(nc_multiset(I$f, I$r))
        }else{
            return(choose(I$n, I$r))
        }
    }
}
