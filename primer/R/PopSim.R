`PopSim` <-
function (Rs, N0, years = 50, sims = 10) 
{
    sim.RM = matrix(sample(Rs, size = sims * years, replace = TRUE), 
        nrow = years, ncol = sims)
    output <- numeric(years + 1)
    output[1] <- N0
    outmat <- sapply(1:sims, function(i) {
        for (t in 1:years) output[t + 1] <- round(output[t] * 
            sim.RM[t, i], 0)
        output
    })
    return(outmat)
}
