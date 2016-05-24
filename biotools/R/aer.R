aer <-
function(obs, predict)
{
   obs <- as.vector(obs)
   predict <- as.vector(predict)
   if (length(obs) != length(predict))
      stop("incompatible dimensions!")
   n <- length(obs)
   dif <- function(x, y) mean(x != y)
   rate <- dif(obs, predict)
   return(rate)
}
