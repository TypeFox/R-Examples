genZ2freq <-
function (genZ) 
{
    s <- apply(genZ, 2, sum)
    if (ncol(genZ) == 3) {
        return(s/(sum(s)))
    }
    else {
        d <- NULL
        n <- log(ncol(genZ))/log(3)
        for (l in 1:n) {
            d <- c(d, sum(s[(3 * l - 2):(3 * l)]))
        }
        return(s/(rep(d, each = 3)))
    }
}
