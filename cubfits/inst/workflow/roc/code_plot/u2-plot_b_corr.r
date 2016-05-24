### For plotting b corr.

my.range <- function(x){
  lim <- range(x)
  lim + c(-1, 1) * (lim[2] - lim[1]) * 0.20
} # End of my.range().

plot.b.corr <- function(x, y, label, x.ci = NULL, y.ci = NULL,
    xlim = NULL, ylim = NULL, xlab = "", ylab = "", main = "",
    workflow.name = NULL, add.lm = FALSE, add.ci = FALSE){
  if(!is.null(x.ci)){
    x.ci <- matrix(x.ci, ncol = 2)
  }
  if(!is.null(y.ci)){
    y.ci <- matrix(y.ci, ncol = 2)
  }
  if(is.null(xlim)){
    xlim <- my.range(x)
  }
  if(is.null(ylim)){
    ylim <- my.range(y)
  }

  ### main part.
  plot(NULL, NULL, xlim = xlim, ylim = ylim,
       xlab = xlab, ylab = ylab, main = main)
  mtext(workflow.name, line = 3, cex = 0.6)

  ### Add lm.
  if(add.lm){
    m.1 <- try(lm(y ~ x), silent = TRUE)
    if(class(m.1) != "try-error"){
      a <- m.1$coef[1]
      b <- m.1$coef[2]
      R2 <- summary(m.1)$r.squared
      abline(a = a, b = b, col = 2)

      rho <- ifelse(b > 0, sqrt(R2), -sqrt(R2))
      
      width <- xlim[2] - xlim[1]
      height <- ylim[2] - ylim[1]
      text(xlim[1] - width * 0.04, ylim[2] - height * 0.05,
           parse(text = paste("y == ", sprintf("%.4f", a),
                              " + ", sprintf("%.4f", b), " * x", sep = "")),
           pos = 4, cex = 0.6)

      text(xlim[1] - width * 0.04, ylim[2] - height * 0.10,
           parse(text = paste("rho == ",
                              sprintf("%.4f", rho), sep = "")),
           pos = 4, cex = 0.6)

      ### Add 95% CI.
      if(add.ci){
        text(xlim[1] - width * 0.04, ylim[2] - height * 0.15,
             "95% CI ",
             pos = 4, cex = 0.6)

        tmp <- summary(m.1)

        ### For intercept.
        CI.95 <- tmp$coefficients[1, 1] +
                 qt(c(0.025, 0.975), tmp$df[2]) * tmp$coefficients[1, 2]
        text(xlim[1] - width * 0.04, ylim[2] - height * 0.20,
             paste("intercept: ",
                   sprintf("(%.4f, %.4f)", CI.95[1], CI.95[2]), sep = ""),
             pos = 4, cex = 0.6)

        ### For slop.
        CI.95 <- tmp$coefficients[2, 1] +
                 qt(c(0.025, 0.975), tmp$df[2]) * tmp$coefficients[2, 2]
        text(xlim[1] - width * 0.04, ylim[2] - height * 0.25,
             paste("slope: ",
                   sprintf("(%.4f, %.4f)", CI.95[1], CI.95[2]), sep = ""),
             pos = 4, cex = 0.6)
      }
    }
  }

  ### Add one-to-one.
  abline(a = 0, b = 1, col = 4, lty = 2)

  ### Add main lines.
  for(i in 1:length(x)){
    if(!is.null(x.ci)){
      lines(x = x.ci[i,], y = rep(y[i], 2), col = 1)
    }
    if(!is.null(y.ci)){
      lines(x = rep(x[i], 2), y = y.ci[i,], col = 1)
    }
  }

  ### Build index by codons.
  tmp.id <- order(x)
  tmp.x <- x[tmp.id]
  tmp.y <- y[tmp.id]
  tmp.label <- label[tmp.id]
  tmp.pch <- rep(1, length(tmp.x))
  tmp.col <- rep(1, length(tmp.x))
  ### 2-codons amino acids: C, D, E, F, H, K, N, Q, Y, Z.
  id.aa <- c("C", "D", "E", "F", "H", "K", "N", "Q", "Y", "Z")
  id.pch <- c(0, 1, 2, 5, 6, 21:25, 15:18)
  for(i.id.aa in 1:length(id.aa)){
    id.codon <- grep(paste("^", id.aa[i.id.aa], "\\.", sep = ""), tmp.label)
    tmp.pch[id.codon] <- id.pch[i.id.aa]
  }
  ### 3-codons amino acids: I. 
  ### 4-codons amino acids: A, G, P, T, V.
  ### 6-codons amino acids: L, R.
  id.aa <- c("I", "A", "G", "P", "T", "V", "L", "R")
  id.pch <- c(15:18, 7:9, 11)
  id.col <- c(2, 3, 4, 5, 6, 8, 2, 3, 4)
  for(i.id.aa in 1:length(id.aa)){
    id.codon <- grep(paste("^", id.aa[i.id.aa], "\\.", sep = ""), tmp.label)
    tmp.pch[id.codon] <- id.pch[i.id.aa]
    tmp.col[id.codon] <- id.col[i.id.aa]
  }

  ### Add points.
  if(is.null(x.ci) || is.null(y.ci)){
    # points(tmp.x, tmp.y, pch = 20, cex = 0.6, col = 2)
    points(tmp.x, tmp.y, pch = tmp.pch, cex = 0.8, col = tmp.col)
  }

  ### Add label.
  # x.split <- xlim[1] + (xlim[2] - xlim[1]) / 2
  # x.adj <- (xlim[2] - xlim[1]) * 0.1 *
  #          rep(c(1, -1), c(sum(tmp.x <= x.split), sum(tmp.x > x.split))) *
  #          ((1:length(tmp.x) - 1) %% 3 + 1)
  x.adj <- (xlim[2] - xlim[1]) * 0.1 *
           rep(c(1.2, -1.4, 1.4, -1.6),
               as.integer(length(x) / 4) + 1)[1:length(x)]
  tmp.x <- tmp.x + x.adj
  text(tmp.x, tmp.y, labels = tmp.label, cex = 0.6, col = tmp.col)
} # End of plot.b.corr().

