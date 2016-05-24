renamefields <-
function (dat, ID = "ID", x = "x", y = "y", Species = "Species") 
{
    cn <- names(dat)
    m0 <- match(ID, cn, nomatch = 0)
    m1 <- match(x, cn, nomatch = 0)
    m2 <- match(y, cn, nomatch = 0)
    m3 <- match(Species, cn, nomatch = 0)
    cn[m0] <- "ID"
    cn[m1] <- "x"
    cn[m2] <- "y"
    cn[m3] <- "Species"
    names(dat) <- cn
    return(dat)
}
