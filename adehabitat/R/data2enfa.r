"data2enfa" <- function (kasc, pts)
{
    ## Verifications
    if (!inherits(kasc, "kasc"))
        stop("should be an object of class \"kasc\"")
    if (ncol(pts) != 2)
        stop("pts should have 2 columns")

    ## prepares the output
    attr <- storemapattr(kasc)
    tab <- kasc2df(kasc)
    index <- tab$index
    tab <- tab$tab
    pr <- as.vector(count.points(pts, kasc))[index]
    dataenfa <- list(tab = data.frame(tab), pr = pr,
        index = index, attr = attr)
    class(dataenfa) <- "dataenfa"

    ## output
    return(invisible(dataenfa))
}

