partial.mk.test <-
function(x, z ){
    # Tests
    na.fail(x)
    na.fail(z)
    n1 <- length(x)
    n2 <- length(z)
    if (n1 != n2) {
        stop("Error: objects x and y must be of same length")
    }

    # Ausgabe fuer Klasse 'htest'
    res <- list(statistic=NULL,
 #               parameter = NULL,
                estimate = NULL,
                p.value =NULL,
                statistic =NULL,
                alternative = "true trend exists in series",
                method = "Partial Mann-Kendall trend test",
                data.name = NULL)


    DNAMEX <- deparse(substitute(x))
    DNAMEY <- deparse(substitute(z))

    inp <- data.frame(x,z)

    out <- multivar.MK.test(inp, method = "Partial")

    Z <- out$Z
    names(Z) <- "Z"
    res$statistic <- Z
    p.value <- out$pvalue
    names(p.value) <- "p-value"
    res$p.value <- p.value
    res$data.name <- paste("t AND ",DNAMEX," . ",DNAMEY)
    S <- out$Stot
    names(S) <- "S"
    varS <- out$Varianz
    names(varS) <- "varS"
    res$estimate <- c(S, varS)
    class(res) <- "htest"
    return(res)
}

