con.regmat <- function (regmat, neighbors, count = FALSE) 
{
    nart <- ncol(regmat)
    nreg <- nrow(regmat)
    out <- c()
    for (i in 1:nart) {
        nspec <- sum(regmat[, i])
        if (nspec>0){
          if (count) 
            cat("Species ", i, " size ", nspec, "\n")
          regions <- (1:nreg)[as.logical(regmat[, i])]
          comat <- matrix(FALSE, ncol = nspec, nrow = nspec)
          for (j in 1:nspec)
            for (k in neighbors[[regions[j]]])
              comat[j, (1:nspec)[regions == k]] <- TRUE
          out[i] <- max(con.comp(comat))
        }
        else
          out[i] <- NA
    }
    out
}



