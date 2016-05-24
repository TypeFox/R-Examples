heavy.control <-
function(maxIter = 500, tolerance = 1e-6, fix.shape = FALSE, ndraws = 500, algorithm = c("EM", "NEM"), ncycles = 5)
{
    algorithm <- match.arg(algorithm)
    choice <- switch(algorithm, "EM" = 0, "NEM" = 1)
    list(maxIter = maxIter, tolerance = tolerance, fix.shape = fix.shape, ndraws = ndraws,
         algorithm = choice, ncycles = ncycles)
}
