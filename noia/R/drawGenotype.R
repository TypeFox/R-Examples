drawGenotype <-
function (nloc = 1, type = "F2", freqmat = NULL) 
{
    ans <- ""
    for (l in 1:nloc) {
        r <- runif(1)
        g <- ""
        if (type == "F2") {
            if (r < 1/4) {
                g <- "1"
            }
            else if (r < 3/4) {
                g <- "2"
            }
            else {
                g <- "3"
            }
        }
        else if (type == "Finf") {
            if (r < 1/2) {
                g <- "1"
            }
            else {
                g <- "3"
            }
        }
        else if (type == "F1") {
            g <- "2"
        }
        else if (type == "UWR") {
            if (r < 1/3) {
                g <- "1"
            }
            else if (r < 2/3) {
                g <- "2"
            }
            else {
                g <- "3"
            }
        }
        else if (type == "G2A") {
            if (is.numeric(freqmat) && length(freqmat) == nloc) {
                if (r < freqmat[l]^2) {
                  g <- "1"
                }
                else if (r < freqmat[l]^2 + 2 * freqmat[l] * 
                  (1 - freqmat[l])) {
                  g <- "2"
                }
                else {
                  g <- "3"
                }
            }
            else {
                stop("freqmat argument does not fit population type \"G2A\", must be numeric and have length()==nloc")
            }
        }
        else if (type == "noia") {
            if (is.numeric(freqmat) && all(dim(freqmat) == c(nloc, 
                3))) {
                if (r < freqmat[l, 1]) {
                  g <- "1"
                }
                else if (r < freqmat[l, 1] + freqmat[l, 2]) {
                  g <- "2"
                }
                else {
                  g <- "3"
                }
            }
            else {
                stop("freqmat argument does not fit population type \"noia\", must be numeric and have dim()==c(nloc,3)")
            }
        }
        else {
            stop("Population type ", type, " unknown.")
        }
        ans <- paste(ans, g, sep = "")
    }
    return(ans)
}
