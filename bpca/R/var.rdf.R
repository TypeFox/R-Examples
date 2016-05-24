var.rdf <- function(x,
                    var.rb,
                    limit)
{
  dif <- 100 * abs(var.rb - cor(x))

  big <- dif > limit
  leq <- dif <= limit

  dif[big] <- '*'
  dif[leq] <- ''
  diag(dif) <- '-'

  var.rd <- dif
  dimnames(var.rd) <- list(dimnames(x)[[2]],
                           dimnames(x)[[2]])

  return(as.data.frame(var.rd))
}
