# This file creates a hyperparameter object.
# It is designed to be called by file apprex_1d.R


"hpa.fun.1d" <-
function (x) 
{
    if (length(x) != 9) {
        stop("x must have 11 elements")
    }
    "pdm.maker" <- function(x) {
        jj <- diag(x,nrow=1)
        rownames(jj) <- LETTERS[1]
        colnames(jj) <- LETTERS[1]
        return(jj)
    }
    sigma_squareds <- x[1:3]
    names(sigma_squareds) <- paste("level", 1:3, sep = "")
    B <- list()
    B[[1]] <- pdm.maker(x[5])
    B[[2]] <- pdm.maker(x[6])
    B[[3]] <- pdm.maker(x[7])
    rhos <- x[8:9]
    names(rhos) <- paste("level", 1:2, sep = "")
    return(list(sigma_squareds = sigma_squareds, B = B, rhos = rhos))
}

