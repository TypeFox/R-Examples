engen2008II <- function(us, av, id, nsim=500, nsimra=500)
{
    if (is.data.frame(us))
        us <- as.matrix(us)
    if (is.data.frame(av))
        av <- as.matrix(av)
    if (ncol(us)!=ncol(av))
        stop("us and av should have the same number of column")
    if (any(colnames(us)!=colnames(av)))
        stop("us and av should have the same column names")
    if (length(id)!=nrow(us))
        stop("id should have same length as the number of rows in us")
    ma <- apply(av,2, function(x) max(table(x)))
    ma2 <- apply(us,2, function(x) max(table(x)))
    ma <- max(c(ma,ma2))
    if (max(ma)==1) {
        nsimra <- 1
        warning("there were no ties in the data,\n no randomizations were carried out")
    }
    if (!is.factor(id))
        id <- factor(id)
    if (min(table(id))<2)
        stop("at least two used units are required")

    toto <- .C("engen2008r", as.double(t(av)), as.double(t(us)),
               as.integer(nrow(av)), as.integer(nrow(us)),
               as.integer(ncol(av)), as.integer(as.numeric(id)),
               as.integer(nlevels(id)), as.integer(nsim),
               double(ncol(av)*5*nsimra), as.integer(nsimra),
               PACKAGE="adehabitatHS")

    to <- as.data.frame(matrix(toto[[9]], nrow=nsimra, byrow=TRUE))
    res <- list()
    res$raw <- lapply(1:ncol(av), function(i) {
        oo <- as.data.frame(to[,c(((i-1)*5+1):((i-1)*5+5))])
        names(oo) <- c("sigma2", "tau", "rho", "mu", "sd.mu")
        return(oo)
    })
    names(res$raw) <- colnames(av)

    res$results <- as.data.frame(do.call("rbind", lapply(res$raw,
                                           function(x) apply(x,2,mean))))
    res$results <- res$results[,-c(1:2)]

    res$results[,3] <- sqrt(res$results[,3])

    class(res) <- "engenetalII"
    return(res)
}


print.engenetalII <- function(x, ...)
{
    if (!inherits(x, "engenetalII"))
        stop("x should be of class \"engenetalII\"")
    cat("*************************************\n")
    cat("** Method of Engen et al. (2008)\n\n")
    cat("Preferences and correlation:\n")
    print(x$results)
    cat("\nThis data frame is stored in the component $result of this list.\n\n")
}





engen2008I <- function(us, av, nsimra=500)
{
    if (is.data.frame(us))
        us <- as.matrix(us)
    if (is.data.frame(av))
        av <- as.matrix(av)
    if (ncol(us)!=ncol(av))
        stop("us and av should have the same number of column")
    if (any(colnames(us)!=colnames(av)))
        stop("us and av should have the same column names")
    ma <- apply(av,2, function(x) max(table(x)))
    ma2 <- apply(us,2, function(x) max(table(x)))
    ma <- max(c(ma,ma2))
    if (max(ma)==1) {
        nsimra <- 1
        warning("there were no ties in the data,\n no randomizations were carried out")
    }

    toto <- .C("engen2008Ir", as.double(t(av)), as.double(t(us)),
               as.integer(nrow(av)), as.integer(nrow(us)),
               as.integer(ncol(av)), double(ncol(av)*2*nsimra),
               as.integer(nsimra), PACKAGE="adehabitatHS")

    to <- as.data.frame(matrix(toto[[6]], nrow=nsimra, byrow=TRUE))
    res <- list()
    res$raw <- lapply(1:ncol(av), function(i) {
        oo <- as.data.frame(to[,c(((i-1)*2+1):((i-1)*2+2))])
        names(oo) <- c("mu", "sd.mu")
        return(oo)
    })
    names(res$raw) <- colnames(av)

    res$results <- as.data.frame(do.call("rbind", lapply(res$raw,
                                                         function(x) apply(x,2,mean))))
    res$results <- res$results[,c(1:2)]

    res$results[,2] <- sqrt(res$results[,2])

    class(res) <- "engenetalI"
    return(res)
}


print.engenetalI <- function(x, ...)
{
    if (!inherits(x, "engenetalI"))
        stop("x should be of class \"engenetalI\"")
    cat("*************************************\n")
    cat("** Method of Engen et al. (2008)\n\n")
    cat("Preferences:\n")
    print(x$results)
    cat("\nThis data frame is stored in the component $result of this list.\n\n")
}
