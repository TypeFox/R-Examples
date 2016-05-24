rocs <- function(Y, pv, cex = 1){
  L <- dim(pv)[2]
  par(mai = c(0.1, 0.1, 0.1, 0.1), cex = cex)
  plot(c(0, L + 0.5, L + 0.5, 0, 0), c(0, 0, L + 0.5, L + 0.5, 0),
       type = "l", axes = FALSE, ann = FALSE,col="blue")
  lines(c(0, 0.5), c(L + 0.5, L), col = 'blue')
  text(0.35, L + 0.35, bquote(theta), cex = cex)
  text(0.15, L + 0.15, bquote(b), cex = cex)
  w <- nchar(attr(Y, 'levels'))
  if(any(w > 7)) {
    w <- nchar(attr(Y, 'levels'))
  }
  for(b in seq_len(L)) {
    lines(c(0, L + 0.5), c(b, b), col = "blue")
    lines(c(b - 0.5, b - 0.5), c(0, L + 0.5), col="blue")
    text(b, L + 0.25, bquote(.(attr(Y, 'levels')[b])), cex = cex)
    text(0.25, b - 0.5, bquote(.(substr(attr(Y, 'levels'), 1, 7)[L - b + 1])), cex = cex)
    Nb = sum(Y == b)
    for(a in seq_len(L)) {
      lines(c(a - 0.5, a + 0.5), c(L - b, L - b + 1), col = "green")
      tmpx <- sort(c(0, pv[Y == b, a], pv[Y == b, a], 1));
      tmpy <- matrix(0, 2 * Nb + 2, 1)
      tmpy[seq(3, 2 * Nb + 2, 2)] = seq(1 / Nb, 1, by = 1 / Nb)
      tmpy[seq(4, 2 * Nb + 2, 2)] = seq(1 / Nb, 1, by = 1 / Nb)
      points(a - 0.5 + tmpx, L - b + tmpy, type = "s", lwd = 3)
    }
  }
}
  
