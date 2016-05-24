bootstrapsamples <-
function(data, numSamples)
{
    boot <- boot(data, function(d, f) f, numSamples)
    dimension <- dim(data)
    samples <- list()
    noSamples <- list()
    array(boot.array(boot, indices = FALSE), c(numSamples, dimension[1])) -> howMany
    array(boot.array(boot, indices = TRUE), c(numSamples, dimension[1])) -> who
    result <- list(who, howMany)
    names(result) <- c("who", "howMany")
    return(result)
}