"boa.geweke" <-
function(link, p.first, p.last)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
{
   iter <- unique(boa.iter(link))
   n <- length(iter)
   link.first <- boa.getiter(link, iter[1:round(p.first * n)])
   link.last <- boa.getiter(link, iter[(n - round(p.last * n) + 1):n])
   result <- (colMeans(link.first) - colMeans(link.last)) /
                       sqrt(boa.gewekePwr(link.first) / nrow(link.first) +
                            boa.gewekePwr(link.last) / nrow(link.last))
   result <- cbind(result, 2 * (1 - pnorm(abs(result))))
   dimnames(result)[[2]] <- c("Z-Score", "p-value")

   return(result)
}
