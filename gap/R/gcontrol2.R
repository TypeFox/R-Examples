gcontrol2 <- function(p,col=palette()[4],lcol=palette()[2],...)
{
   p <- p[!is.na(p)]
   n <- length(p)
   x2obs <- qchisq(p,1,lower.tail=FALSE)
   x2exp <- qchisq(1:n/n,1,lower.tail=FALSE)
   lambda <- median(x2obs)/median(x2exp)
   qqunif(p,col=col,lcol=lcol,...)
   invisible(list(x=x2exp,y=x2obs,lambda=lambda))
}
