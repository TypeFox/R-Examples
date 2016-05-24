symclos <-
function (q) 
{
    for (i in 1:nrow(q)) {
        for (j in 1:ncol(q)) {
            if (isTRUE(q[i, j] != q[j, i]) == TRUE) {
                if (isTRUE(q[i, j] == "o") == TRUE) {
                  q[i, j] <- q[j, i]
                }
                else if (isTRUE(q[j, i] == "o") == TRUE) {
                  q[j, i] <- q[i, j]
                }
            }
        }
        rm(j)
    }
    rm(i)
    for (i in 1:nrow(q)) {
        for (j in 1:ncol(q)) {
            if (isTRUE(q[i, j] != q[j, i]) == TRUE) {
                if (isTRUE(q[i, j] == "p") == TRUE) {
                  q[j, i] <- q[i, j]
                }
                else if (isTRUE(q[j, i] == "p") == TRUE) {
                  q[i, j] <- q[j, i]
                }
                else if (isTRUE(q[i, j] == "q") == TRUE || isTRUE(q[j, 
                  i] == "q") == TRUE) {
                  q[i, j] <- q[j, i] <- "p"
                }
                else {
                  q[i, j] <- "a"
                }
            }
        }
        rm(j)
    }
    rm(i)
    return(q)
}
