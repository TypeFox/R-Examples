kulife.colors <- function (n=4, KUgrey=FALSE)
{

  n <- as.integer(n[1L])

  if (n<1)
    n <- 1
  
  if (KUgrey)
    res <- palette(c("black", "#666666", "#541800", paste("#541800", round(seq(100, 20, length=n)), sep="")[-1]))
  else
    res <- palette(c("black", "#541800", paste("#541800", round(seq(100, 20, length=n)), sep="")[-1]))

  return(res)
}
