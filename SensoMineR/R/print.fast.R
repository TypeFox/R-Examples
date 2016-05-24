print.fast=function (x, file = NULL, sep = ";", ...)
{
    res.catego <- x
    if (!inherits(res.catego, "fast"))
        stop("non convenient data")
    cat("**Results for the categorization**\n\n")
    cat("The analysis was done on ", nrow(res.catego$ind$coord),
        "products, described by", nrow(res.catego$group$coord),
        "consumers\n\n")
    cat("*The results are available in the following objects:\n\n")
    res <- array("", c(24, 2), list(1:24, c("name", "description")))
    res[1, ] <- c("$eig", "eigenvalues")
    res[2, ] <- c("$var", "results for the levels")
    res[3, ] <- c("$var$coord", "coordinates of the levels")
    res[4, ] <- c("$var$cos2", "cos2 for the levels")
    res[5, ] <- c("$var$contrib", "contributions for the levels")
    res[6, ] <- c("$var$vtest", "vtest for the levels")
    res[7, ] <- c("$ind", "results for the products")
    res[8, ] <- c("$ind$coord", "coord. for the products")
    res[9, ] <- c("$ind$cos2", "cos2 for the products")
    res[10, ] <- c("$ind$contrib", "contributions for the products")
    res[11, ] <- c("$group", "results for the consumers")
    res[12, ] <- c("$group$coord", "coord. for the consumers")
    res[13, ] <- c("$group$cos2", "cos2 for the consumers")
    res[14, ] <- c("$group$contrib", "contributions for the consumers")
    res[15, ] <- c("$reord", "reordered data set")
    res[16, ] <- c("$cooccur", "cooccurrence matrix")
    res[17, ] <- c("$cramer", "Cramer's coefficient between consumers")
    res[18, ] <- c("$textuel", "Characterization of products by the words")
    print(res[1:18, ])
    if (!is.null(file)) {
        write.infile(res.catego, file = file, sep = sep)
        print(paste("All the results are in the file", file))
    }
}
