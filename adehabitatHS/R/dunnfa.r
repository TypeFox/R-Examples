dunnfa <- function(dudi, pr, scannf = TRUE, nf = 2) #(Z, pr)
{
    ## Verifications
    if (!inherits(dudi, "pca"))
        stop("Currently only implemented for objects of class \"pca\"")
    call <- match.call()
    call <- match.call()
    if (any(is.na(dudi$tab)))
        stop("na entries in table")
    if (!is.vector(pr))
        stop("pr should be a vector")

    ## Bases for the analysis
    prb <- pr
    pr <- pr/sum(pr)
    Z <- as.matrix(dudi$tab)
    n <- nrow(Z)


    ## centering and scaling ("used" weighting)
    f1 <- function(v) sum(v * pr)
    center <- apply(Z, 2, f1)
    Zu1 <- sweep(Z, 2, center)
    f2 <- function(v) sum((v^2) * pr)
    sdu <- apply(Zu1, 2, f2)
    Zu <- sweep(Zu1, 2, sqrt(sdu), "/")

    ## centering for "available" weighting
    center <- apply(Zu, 2, mean)
    Za <- sweep(Zu, 2, center)

    ## correlation matrices
    Su <- t(apply(Zu, 2, function(x) x*pr))%*%Zu
    Sa <- t(Za)%*%Za/nrow(Za)

    ## inverse of Su
    Sumo <- solve(Su)

    ## Cholesky factorization
    Th <- chol(Sa)

    ## The core of the analysis:
    mat <- Th%*%Sumo%*%t(Th)
    res <- eigen(mat)


    ## Number of eigenvalues
    if (scannf) {
        barplot(res$value)
        cat("Select the number of axes: ")
        nf <- as.integer(readLines(n = 1))
    }
    if (nf <= 0 | nf > ncol(Zu))
        nf <- 1


    ## The results
    B <- res$vec
    A <- solve(Th)%*%B

    ## scale the vectors so that their length is 1
    A <- apply(A,2, function(x) x/sqrt(sum(x^2)))

    ## The results:
    sor <- list()
    sor$eig <- res$val
    sor$pr <- prb
    sor$co <- A[,1:nf]
    sor$liU <- as.data.frame(Zu %*% A)[,1:nf]
    sor$liA <- as.data.frame(Za %*% A)[,1:nf]
    sor$corA <- cor(as.data.frame(Za), as.data.frame(sor$liA))
    if (length(dim(sor$liU))>1) {
        sor$mahasu <- apply(sor$liU, 1, function(x) sum(x^2/sor$eig[1:nf]))
    } else {
        sor$mahasu <- sor$liU^2/sor$eig[1]
    }
    sor$call <- call
    sor$tab <- as.data.frame(Za)
    sor$nf <- nf
    class(sor) <- "dunnfa"
    return(sor)
}



print.dunnfa <- function (x, ...)
{
    if (!inherits(x, "dunnfa"))
        stop("Object of class 'dunnfa' expected")
    cat("Factorial analysis of the specialization of James Dunn")
    cat("\n$call: ")
    print(x$call)
    cat("\neigen values: ")
    l0 <- length(x$eig)
    cat(signif(x$eig, 4)[1:(min(5, l0))])
    if (l0 > 5)
        cat(" ...")
    cat("\n$nf:", x$nf, "axes saved")
    cat("\n")
    cat("\n")
    sumry <- array("", c(3, 4), list(1:3, c("vector", "length",
                                            "mode", "content")))
    sumry[1, ] <- c("$pr", length(x$pr), mode(x$pr), "vector of presence")
    sumry[2, ] <- c("$mahasu", length(x$mahasu), mode(x$mahasu), "squared Mahalanobis distances")
    sumry[3, ] <- c("$eig", length(x$eig), mode(x$eig), "eigen values")
    class(sumry) <- "table"
    print(sumry)
    cat("\n")
    sumry <- array("", c(5, 4), list(1:5, c("data.frame", "nrow",
                                            "ncol", "content")))
    sumry[1, ] <- c("$tab", nrow(x$tab), ncol(x$tab), "modified array (centering on available points)")
    sumry[2, ] <- c("$liA", nrow(x$liA), ncol(x$liA), "row coordinates (centering on available points)")
    sumry[3, ] <- c("$liU", nrow(x$liU), ncol(x$liU), "row coordinates (centering on used points")
    sumry[4, ] <- c("$co", nrow(x$co), ncol(x$co), "column coordinates")
    sumry[5, ] <- c("$cor", nrow(x$cor), ncol(x$cor), "cor(habitat var., scores) for available points")
    class(sumry) <- "table"
    print(sumry)
    if (length(names(x)) > 15) {    cat("\nother elements: ")
                                    cat(names(x)[16:(length(x))], "\n")
                                }
}
