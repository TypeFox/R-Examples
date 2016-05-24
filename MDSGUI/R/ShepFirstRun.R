ShepFirstRun <-
function (ShepComp, tShepx, Shepx) 
{
    f.DistFunc <- as.vector(0)
    ShepPointindex <- matrix(nrow = length(tShepx), ncol = 2)
    rownames(ShepPointindex) <- names(tShepx)
    ShepPointindex[, 2] <- seq(1:length(tShepx))
    for (i in 1:length(tShepx)) {
        for (j in 1:length(Shepx)) {
            if (names(tShepx)[i] == names(Shepx)[j]) {
                ShepPointindex[i, 1] <- j
                f.DistFunc[j] <- ShepComp$yf[i]
                names(f.DistFunc)[j] <- names(Shepx)[j]
            }
        }
    }
    return(list(ShepPointindex = ShepPointindex, f.DistFunc = f.DistFunc))
}
