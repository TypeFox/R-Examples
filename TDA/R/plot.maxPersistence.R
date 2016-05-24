plot.maxPersistence <-
function(x, features = "dimension", colorBand = "pink", colorBorder = NA, ...) {
  
  parameter <- x[["parameters"]]
  Pers <- x[["Persistence"]]
  eps <- x[["bands"]]
  Kseq <- length(parameter)
  maxPers <- 0
  for (k in seq_len(Kseq)) {
    maxPers <- max(maxPers, Pers[[k]][, 2]) 
  }

  if (x[["bandFUN"]] == "bootstrapBand" || features == "all")
  {
    dimension <- NULL
  } else {
    dimension <- x[["dimension"]]
  }

  graphics::plot(parameter, rep(maxPers, Kseq), type = "n",
      ylim = c(0, 1.18 * maxPers), xlab = "", ylab = "", axes=FALSE, ...)
  graphics::title(xlab = "parameter", ylab = "Persistence", line = 2.2)
  graphics::axis(1)
  graphics::axis(2)
  
  eps <- pmin(eps, rep(1.1 * maxPers / 2, length(eps)))
  graphics::polygon(c(parameter, parameter[Kseq:1]), c(2 * eps, rep(0, Kseq)),
      col = colorBand, border = colorBorder, lwd = 1.5)
  
  for (i in seq_len(Kseq)) {  
    
    if (x[["bandFUN"]] == "bootstrapBand" || features == "all")
    {
      selectDimension <- seq_len(nrow(Pers[[i]]))
    } else {
      selectDimension <- which(Pers[[i]][, 1] == dimension)
    }

    symb <- Pers[[i]][, 1]
    for (j in seq(along = symb)){
      if (symb[j] == 0) {
        symb[j] <- 16
      } else if (symb[j] == 1) {
        symb[j] <- 2
      } else if (symb[j] == 2) {
        symb[j] <- 5
      }
    }

    col <- Pers[[i]][, 1] + 1  # betti0 black, betti1 red
    for (j in seq(along = symb)) {
      if (col[j] == 3) {
        col[j] <- 4            # betti2 blue
      }
      if (symb[j] == 3) {
        col[j] <- 3            # betti3 green
      }
    }

    graphics::points(rep(parameter[i], length(selectDimension)),
        Pers[[i]][selectDimension, 2], col = col[selectDimension],
        pch = symb[selectDimension], lwd = 2) 
    if (x[["bandFUN"]] == "bootstrapBand" || x[["dimension"]] == 0) {
        graphics::points(parameter[i], 1.2 * maxPers, col = 1, pch = 16,
            lwd = 2)
    } 
  }
  
  
  if (x[["bandFUN"]] == "bootstrapBand" || x[["dimension"]] == 0) {
    graphics::abline(h = 1.18 * maxPers, lty = 2)
  }
  
  # legend(parameter[floor(0.6*Kseq)], maxPers, c("components", "loops", "voids"), pch=c(16,2,5), pt.lwd=2, col=c(1,2,4))
    
}
