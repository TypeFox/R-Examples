weighted.residuals.manylm <- function(obj, drop0 = TRUE){
    w <- weights(obj)
    r <- as.matrix(residuals(obj, type = "deviance"))
    if (drop0 && !is.null(w)) 
        r[w != 0,,drop=FALSE]
    else r

}

#setMethod("weighted.residuals", "manylm", weighted.residuals.manylm)


