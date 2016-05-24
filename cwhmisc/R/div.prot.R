div.prot <- function(num,den) {
  res <- num/den
  res[is.infinite(res)] <- .Machine$double.xmax^(3/4)
  return(res)
}
