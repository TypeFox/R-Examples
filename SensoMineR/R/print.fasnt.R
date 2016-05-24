print.fasnt=function (x, file = NULL, sep = ";", ...){
    res.fasnt <- x
    if (!inherits(res.fasnt, "fasnt"))
        stop("non convenient data")

    cat("**Results for the sorted napping data**\n\n")
    cat("The analysis was done on", nrow(res.fasnt$ind$coord), "objects, described by",
        nrow(res.fasnt$group[[2]]), "subjects\n\n")
    cat("*The results are available in the following objects:\n\n")
    res <- array("", c(18, 2), list(1:18, c("name", "description")))
    res[1, ] <- c("$eig", "eigenvalues")
    res[2, ] <- c("$ind", "results for the objects")
    res[3, ] <- c("$ind$coord", "coordinates for the objects")
    res[4, ] <- c("$ind$cos2", "cos2 for the objects")
    res[5, ] <- c("$ind$contrib", "contributions for the objects")
    res[6, ] <- c("$ind$partial", "results for the partial points associated with the objects")
    res[7, ] <- c("$quali.var", "results for the categories of the qualitative variables")
    res[8, ] <- c("$quali.var$coord", "coordinates for the categories of the qualitative variables")
    res[9, ] <- c("$quali.var$contrib", "contributions for the categories of the qualitative variables")
    res[10, ] <- c("$quanti.var", "results for the quantitative variables")
    res[11, ] <- c("$quanti.var$coord", "coordinates for the quantitative variables")
    res[12, ] <- c("$quanti.var$cos2", "cos2 for the quantitative variables")
    res[13, ] <- c("$quanti.var$contrib", "contributions for the quantitative variables")
    res[14, ]  <- c("$group", "results for the subjects")
    res[15, ]  <- c("$indicator", "indicators associated with Napping and sorting task")
    res[16, ]  <- c("$textual", "textual analysis")
    res[17, ]  <- c("$validity", "results for the elements of validity")
    res[18, ]  <- c("$call", "some statistics")

    print(res[1:18, ])
    if (!is.null(file)) {
        write.infile(res.fasnt, file = file, sep = sep)
        print(paste("All the results are in the file", file))
    }
}
