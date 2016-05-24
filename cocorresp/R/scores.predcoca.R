`scores.predcoca` <- function(x, choices = c(1,2),
                              display = c("sites", "species"), ...) {
    if (!inherits(x, "predcoca"))
        stop("x must be of class \"predcoca\"")
    scoreOpts <- c("species", "sites")
    names(scoreOpts) <- c("species", "sites")
    take <- scoreOpts[display]
    retval <- list()
    retval$species <- if ("species" %in% take)
        list(Y = x$scores$species$Y[, choices, drop = FALSE],
             X = x$scores$species$X[, choices, drop = FALSE])
    retval$sites <- if ("sites" %in% take)
        list(Y = x$scores$site$Y[, choices, drop = FALSE],
             X = x$scores$site$X[, choices, drop = FALSE])
    ##class(retval) <- "scores.predcoca"
    retval
}
