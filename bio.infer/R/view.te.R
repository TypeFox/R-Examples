"view.te" <-
function(coef, plotform = "pdf") {

  done0 <- tkmessageBox(message="View all taxon-environment relationships?",
                        icon="question",type="yesno",default="yes")
  if (as.character(done0) != "yes") {
    tnames0 <- tklist.modal("Select taxa.", coef$tnames, selectmode = "multiple")
  }
  else {
    tnames0 <- coef$tnames
  }

  if (! is.null(coef$raw.data)) {
    nlist <- as.list(rep(NA, times = length(coef$raw.data)))
    for (i in 1:length(coef$raw.data)) {
      nlist[[i]] <- names(coef$raw.data[[i]])
    }
  }

  if (plotform == "pdf") {
    pdf(file = "taxon.env.pdf", width = 9, height = 6.5, pointsize = 10)
    par(mfrow = c(2,3), pty = "s", mar = c(4,4,3,2))

  }
  else {
    if (plotform == "windows") {
#      windows(width = 6, height = 3, pointsize = 10)
      par(mfrow = c(1,2), pty = "m", mar = c(4,4,2,0))
    }
  }
  np <- 100

  xnew <- seq(from = 0, to = 1, length = np)
  ynew <- seq(from = 0, to = 1, length = np)
  df1 <- expand.grid(x = xnew, y = ynew)
  if (length(tnames0) > 0) {
    for (i in 1:length(tnames0)) {
      if ((plotform == "windows") & (i > 2) & (floor((i-1)/2) == (i-1)/2)) {
          dev.new()
        par(mfrow = c(1,2), pty = "m", mar = c(4,4,2,0))
      }
      isel <- match(tnames0[i], coef$tnames)
      cc <- coef$csave[isel,]

      if (length(coef$xvar) == 2) {
        z <- cc[1] + cc[2]*df1$x + cc[3]*(df1$x^2) + cc[4]*df1$y +
          cc[5]*(df1$y^2) + cc[6]*df1$x*df1$y
        z2 <- exp(z)/(1+exp(z))
        dim(z2) <- c(np,np)

        contour(xnew, ynew, z2, axes = FALSE)
        at0 <- seq(from = 0, to = 1, length = 5)
        xlims <- coef$xlims
        lab1 <- round(at0*diff(xlims[[1]]) + xlims[[1]][1], digits = 0)
        lab2 <- round(at0*diff(xlims[[2]]) + xlims[[2]][1], digits = 0)
        box()
        axis(1, at = at0, labels = lab1)
        axis(2, at = at0, labels = lab2)
        mtext(toupper(coef$xvar[1]), side = 1, line=2.3)
        mtext(toupper(coef$xvar[2]), side = 2, line = 2.3)
        mtext(tnames0[i], side = 3, line = 0.5)
      }
      else {
        if (length(coef$xvar) == 1) {
          z <- cc[1] + cc[2]*xnew + cc[3]*xnew^2
          z2 <- exp(z)/(1+exp(z))
          ylim0 <- range(z2)
          if (! is.null(coef$raw.data)) {
            imatch <- NA
            j <- 0
            while ((is.na(imatch)) & (j < length(coef$raw.data))) {
              j <- j + 1
              imatch <- match(tnames0[i], names(coef$raw.data[[j]]))
            }

            cutp <- quantile(coef$raw.data[[j]][, coef$xvar[1]],
                             probs = seq(from = 0, to = 1,
                               length = floor(nrow(coef$raw.data[[j]])/40)))
            cutm <- 0.5*(cutp[-1] + cutp[-length(cutp)])
            cutf <- cut(coef$raw.data[[j]][, coef$xvar[1]], cutp,
                        include.lowest = TRUE)
            vals <- tapply(coef$raw.data[[j]][, imatch], cutf,
                           function(x) mean(as.numeric(x>0)))
            ylim0 <- range(c(ylim0, vals), na.rm = TRUE)
          }
          plot(xnew, z2, axes = FALSE, type = "l", xlab = "", ylab = "",
               ylim = ylim0)
          if (! is.null(coef$raw.data)) {
            points(cutm, vals)
          }
          at0 <- seq(from = 0, to = 1, length = 5)
          xlims <- coef$xlims
          lab1 <- round(at0*diff(xlims[[1]]) + xlims[[1]][1], digits = 1)
          box(bty = "l")
          axis(1, at = at0, labels = lab1)
          axis(2)
          mtext(toupper(coef$xvar[1]), side = 1, line = 2.3)
          mtext("Capture Probability", side = 2, line = 2.3)
          mtext(tnames0[i], side = 3, line = 0.5)
        }
        else {
          print("Coefficient file not suitable.")
        }
      }
    }

  }
  if (plotform == "pdf") {
    dev.off()
  }

}

