
## One curve only
"simDR" <- function(mpar, sigma, fct, noSim = 1000, conc, edVec = c(10, 50), seedVal = 20070723)
{
    set.seed(seedVal)

    ## Calculating the true ED values
    lened <- length(edVec)
    edTRUE <- rep(0, lened)
    for (i in 1:lened)
    {
        edTRUE[i] <- fct$edfct(mpar, edVec[i], type="relative")[[1]]
    }

    ## Run simulations
    edMat1 <- array(NA, c(length(conc)-4, 6, lened))
#    edMat2 <- array(NA, c(length(conc)-4, 6, 3))
#    edMat3 <- array(NA, c(length(conc)-4, 6, 3))    
    tempMat <- matrix(NA, noSim, lened)    
    for (i in 5:length(conc))
    {
        cVec1 <- sort(conc[1:i])
        for (j in 1:6)
        {
            cVec2 <- rep(cVec1, rep(j, i))
            sim1 <- rdrm(noSim, LL.4(), mpar, cVec2, ypar = sigma)

            for (k in 1:noSim)
            {
                tempFit <- try(drm(sim1$y[k, ]~sim1$x[k, ], fct = fct), silent = TRUE)
                if (!inherits(tempFit, "try-error"))
                {
                    edVal <- ED(tempFit, edVec, display = FALSE, multcomp = TRUE)[["EDdisplay"]]
                    tempMat[k, ] <- edVal[, 1] - edTRUE
                }
            }
            edMat1[i - 4, j, ] <- apply(tempMat, 2, sd, na.rm = TRUE)
#            edMat2[i - 4, j, ] <- apply(tempMat, 2, mean, na.rm = TRUE)
#            edMat3[i - 4, j, ] <- apply(tempMat, 2, function(x) {mean(x^2, na.rm = TRUE)})
        }
    } 

#print(edMat1)
    cat("Concentrations used:", conc, "\n\n")
    for (i in 1:lened)
    {
        tempMat <- edMat1[, , i]
        colnames(tempMat) <- 1:6
        rownames(tempMat) <- 5:9
        
        cat("ED value considered:", edVec[i], "\n")
        cat("Conc. no.\\Replicates:", "\n")
        print(tempMat)
        cat("\n\n")
    } 
    
    invisible(list(se=edMat1))  # , bias=edMat2, mse=edMat3))
}
