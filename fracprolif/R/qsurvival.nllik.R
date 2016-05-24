qsurvival.nllik <- function(dist, complete.lifespans, censored.lifespans, Q, ...)
{
  -sum(log(1-Q) +
       do.call(paste("d",dist, sep=''), list(x=complete.lifespans, log=TRUE, ...))
       ) -
  sum(log(Q+(1-Q)*
             do.call(paste("p",dist, sep=''), list(q=censored.lifespans, lower.tail=FALSE, ...))
      ))
}
