"hist.enfa" <-
function (x, scores = TRUE, type = c("h", "l"),
                       adjust = 1, Acol, Ucol, 
                       Aborder, Uborder, Alwd = 1, Ulwd = 1, ...)
{
  type <- match.arg(type)
  if (!inherits(x, "enfa")) 
    stop("Object of class 'enfa' expected")
  if (scores)
    tab <- x$li
  else
    tab <- x$tab
  pr <- x$pr
  if (missing(Acol)) {
    Acol <- NULL
    Acolf <- "white"
    Acold <- "black"
  }
  else {
    Acold <- Acol
    Acolf <- Acol
  }
  if (missing(Aborder)) 
    Aborder <- "black"
  if (missing(Ucol)) {
    Ucol <- gray(0.8)
    Ucold <- gray(0.8)
  }
  else
    Ucold <- Ucol
  if (missing(Uborder)) 
    Uborder <- gray(0.8)
  clas <- rep("", ncol(tab))
  for (j in 1:ncol(tab)) {
    w1 <- "q"
    if (is.factor(tab[, j])) 
      w1 <- "f"
    clas[j] <- w1
  }
  if (any(clas == "f") & type == "l") 
    warning("Type = 'l' is not possible for factors, type = 'h' used instead.\n")
  old.par <- par(no.readonly = TRUE)
  on.exit(par(old.par))
  par(mar = c(0.5, 0.5, 2, 0.5))
  par(mfrow = rev(n2mfrow(ncol(tab))))
  f1 <- function(j) {
    tmpU <- rep(tab[, j], pr)
    tmpA <- tab[, j]
    name <- names(tab)[j]
    if (clas[j] == "f") {
      par(mar = c(3, 0.5, 2, 0.5))
      mat <- t(cbind(table(tmpA), table(tmpU)))
      mat <- lapply(1:2, function(i) mat[i, ]/sum(mat[i, 
                                                      ]))
      mat <- rbind(mat[[1]], mat[[2]])
      max <- max(mat)
      max <- max + max/20
      ylim <- c(0, max)
      barplot(mat, col = c(Acolf, Ucol), border = c(Aborder, Uborder), 
              ylim = ylim, main = name, ylab = NULL, axes = FALSE, 
              beside = TRUE, ...)
      par(mar = c(0.5, 0.5, 2, 0.5))
    }
    else {
      if (type == "h") {
        xrange <- range(tmpA)
        H <- hist(tmpU, plot = FALSE, br = seq(min(xrange), 
					max(xrange), length = 15))
        G <- hist(tmpA, plot = FALSE, br = seq(min(xrange), 
					max(xrange), length = 15))
        yrange <- c(0, max(H$density, G$density))
        plot(H, freq = FALSE, col = Ucol, border = Uborder, 
             xlim = xrange, ylim = yrange, main = name, xlab = NULL, 
             ylab = "Density", axes = FALSE, ...)
        plot(G, freq = FALSE, col = Acol, border = Aborder, add = TRUE)
      }
      if (type == "l") {
        densA <- density(tmpA, adjust = adjust)
        densU <- density(tmpU, adjust = adjust, from = min(densA$x), 
                         to = max(densA$x))
        max <- max(densU$y, densA$y)
        max <- max + max/20
        ylim <- c(0, max)
        plot(densU, col = Ucol, ylim = ylim, type = "l", 
             lwd = Ulwd, main = name, xlab = NULL, ylab = "Density", 
             axes = FALSE, ...)
        lines(rep(mean(tmpU), 2), c(0, densU$y[512 - sum(densU$x > 
                                                         mean(tmpU))]),
              col = Ucol, lty = 2, lwd = Ulwd)
        lines(densA, col = Acold, lwd = Alwd)
        lines(rep(mean(tmpA), 2), c(0, densA$y[512 - sum(densA$x > 
                                                         mean(tmpA))]),
              col = Acold, lty = 2, lwd = Alwd)
      }
    }
    box()
  }
  lapply(1:ncol(tab), f1)
  return(invisible(NULL))
}

