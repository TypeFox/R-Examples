plot.dose.radir <-
function(x, ci=FALSE, cr=NA, distr=FALSE, prob=NA, col.ci="grey", col.pr="grey", ...)
{
  if (class(x)!="dose.radir") stop("Wrong object")
  if (ci==TRUE & distr==TRUE) distr <- FALSE
  if (!is.na(prob)[1] & (length(prob) != 2 | prob[1] > prob[2] | prob[1] < 0 | prob[2] < 0)) stop("Wrong probability interval")
  if (ci==TRUE & !is.na(prob)[1]) ci <- FALSE
  if (distr==TRUE & !is.na(prob)[1]) distr <- FALSE
  
  if (distr==FALSE) 
  {
    plot(x[[2]], x[[1]], type="l", ylab="Probability Density", xlab="Dose, x, Gy", ...)
    if (ci==TRUE)
    {
      if (is.na(cr)) cr <- 0.95
      polygon(x=seq(ci.dose.radir(x, cr)[1], ci.dose.radir(x, cr)[2], length.out=length(c(min(x[[1]]),x[[1]][which(x[[2]] > ci.dose.radir(x, cr)[1] & x[[2]] < ci.dose.radir(x, cr)[2])],min(x[[1]])))),
              y=c(min(x[[1]]),x[[1]][which(x[[2]] > ci.dose.radir(x, cr)[1] & x[[2]] < ci.dose.radir(x, cr)[2])],min(x[[1]])), col=col.ci)
    }
    if (!is.na(prob)[1])
    {
      polygon(x=c(prob[1], x[[2]][which(x[[2]] > prob[1] & x[[2]] < prob[2])], prob[2]),
              y=c(min(x[[1]]), x[[1]][which(x[[2]] > prob[1] & x[[2]] < prob[2])], min(x[[1]])), col=col.pr)
    }
  }
  if (distr==TRUE) 
  {
    probs <- sapply(1:length(x[[2]]), FUN=function(i){pr.dose.radir(x, 0, x[[2]][i])})
    plot(x[[2]], probs, type="l", ylab="Probability", xlab="Dose, x, Gy")
  }
}
