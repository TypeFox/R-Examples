nappeplot <- function (donnee, numr = 2, numc = 2, color = "blue", lim = c(60, 40)) {
  nbjuge <- ncol(donnee)/2
  mult <- nbjuge%/%(numr * numc)
  if (nbjuge == (nbjuge%/%(numr * numc)) * (numr * numc))   mult = mult - 1
  for (m in 0:mult) {
    par(mfrow = c(numr, numc)) ###
    for (nbd in 1:(numr * numc)) {
      nb <- (m * (numr * numc) + nbd)
      if (nb <= nbjuge) {
        plot(donnee[, (2 * (nb - 1) + 1):(2 * nb)], col = color,
             xlim = c(0, lim[1]), ylim = c(0, lim[2]), xlab = "", 
             ylab = "", main = paste(colnames(donnee)[2 * nb]), type = "n", asp = 1)
        points(donnee[, (2 * (nb - 1) + 1):(2 * nb)], col = color, cex = 0.8, pch = 20)
        text(donnee[, (2 * (nb - 1) + 1):(2 * nb)], label = rownames(donnee), col = color, cex = 0.8, pos = 4, offset = 0.2)
      }
    }
    if (m < mult) dev.new()
  }
  par(las = 0)
  par(mfrow = c(1, 1))
}
