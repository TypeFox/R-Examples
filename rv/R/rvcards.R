
rvcarddeck <- function (n.jokers=0L) {
  faces.df <- expand.grid(y=c("A", 2:10, "J", "Q", "K"), x=c("H", "S", "D", "C"))
  faces <- do.call(paste, c(list(sep=""), as.list(faces.df)))
  if (n.jokers > 0L) {
    faces <- c(faces, paste("JOKER", 1:n.jokers, sep=""))
  }
  i <- seq_along(faces)
  deck0 <- rvpermut(i)
  deck <- rvfactor(deck0, levels=faces)
  return(deck)
}
