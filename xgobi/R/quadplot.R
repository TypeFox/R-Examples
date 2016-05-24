quadplot <-
    function(mat4,
             pointlabs = rownames(mat4),
             vertexlabs = paste(1:4),
             normalize = median(abs(c(mat4))) > 1)
{
    mat4 <- if(is.data.frame(mat4)) data.matrix(mat4) else as.matrix(mat4)
    n <- nrow(mat4)
    m <- ncol(mat4)
    if(m != 4)
        stop("Data matrix `mat4' must have four columns.")
    if(normalize)
        mat4 <- mat4 / c(mat4 %*% c(1,1,1,1))
    ## == sweep(mat4, 1, apply(mat4,1,sum), "/")

    rt3  <- 1/sqrt(3)
    rt14 <- 1/sqrt(14)
    projct <- cbind(c(0, -1,-4, 5)*(rt3*rt14),
                    c(3, -1,-1,-1)* rt3/2,
                    c(0, -3, 2, 1)* rt14)
    tetralines <- cbind(c(1,1,1,2,2,3),
                        c(2,3,4,3,4,4))

    mat3 <- (rbind(diag(4), mat4) - 1/4) %*% projct

    if(is.null(pointlabs)) pointlabs <- as.character(1:n)
    xgobi(mat3, lines = tetralines, rowlab = c(vertexlabs, pointlabs),
	  resources = c("*showLines: True", "*showAxes: False"))
    invisible(mat3) # or invisible()
}
