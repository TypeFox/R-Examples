prime.factor <-
function (n) 
{
    data(primlist)
    dummy <- 2
    end <- 0
    faclist <- 0
    test <- is.prim(n)
    if (test == TRUE) {
## The two following lines were modified:
        msg <- cat(n, "is a prime!\n")  ## Initially \r instead of \n, why?
        return(n)                       ## initially return(msg)
## End of modifications
    }
    while (end == 0) {
        prim <- primlist[dummy]
        if (prim > n) {
            end <- 1
        }
        else {
            test <- n/prim
            test2 <- is.whole(test)
            if (test2 == TRUE) {
                faclist <- c(faclist, prim)
                test3 <- is.prim(test)
                if (test3 == TRUE) {
                  end <- 1
                  faclist <- c(faclist, test)
                }
                else {
                  n <- test
                  dummy <- 1
                }
            }
            dummy <- dummy + 1
        }
    }
    faclist <- faclist[-1]
    return(faclist)
}

