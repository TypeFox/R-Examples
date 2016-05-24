overlap.plot <-
function(over.stat, nbreaks = 50, equal.axis = FALSE, species.names, col,
                         mean.cred = TRUE, mean.cred.col = "green", xlab) {
  if(length(dim(over.stat)) != 3 || dim(over.stat)[1] != dim(over.stat)[2])
    stop("Incorrect specification of over.stat.")
  nspec <- dim(over.stat)[1]
  if(missing(species.names)) species.names <- dimnames(over.stat)[[1]]
  # histograms
  over.hist <- apply(over.stat, 1:2,
                     function(x) {
                       if(any(is.na(x))) return(NULL)
                       else {
                         tmp <- hist(x*100, breaks = nbreaks, plot = FALSE)
                         tmp$density <- tmp$density/max(tmp$density)
                         tmp$counts <- tmp$counts/max(tmp$counts)
                       }
                       tmp
                     })
  # x-axis limits
  xlim <- lapply(over.hist,
                 function(h) {
                   if(is.null(h)) tmp <- matrix(NA, 2, 2)
                   else {
                     tmp <- cbind(range(h$breaks), range(h$density))
                   }
                   tmp
                 })
  xlim <- matrix(xlim, nspec, nspec)
  if(equal.axis) {
    for(ii in 1:nspec) {
      xlim[,ii] <- rep(list(cbind(range(sapply(xlim[,ii],
                                               function(ll) ll[,1]), na.rm = TRUE),
                                  range(sapply(xlim[,ii],
                                               function(ll) ll[,2]), na.rm = TRUE))),
                       nspec)
    }
  }
  # mean and credible intervals
  if(mean.cred) {
    over.mean <- apply(over.stat, 1:2, mean)*100
    over.cred <- apply(over.stat, 1:2, quantile, prob = c(.025, .975), na.rm = TRUE)*100
    over.cred <- array(over.cred, dim = c(2, nspec, nspec))
  }
  # plot
  par(mfcol = c(nspec,nspec), mar = c(1.5,rep(.5, 3)), oma = rep(4,4))
  for(ci in 1:nspec) {
    for(ri in 1:nspec) {
      # initialize plot
      plot.new()
      if (ri != ci) {
        # off diagonals: overlap histograms
        plot.window(xlim[ri,ci][[1]][,1], xlim[ri,ci][[1]][,2])
        if(equal.axis) {
          # recalculate histograms with new limits
          tmp <- hist(over.stat[ri,ci,]*100,
                      breaks = seq(xlim[ri,ci][[1]][1,1],
                        xlim[ri,ci][[1]][2,1], len = nbreaks+1),
                      plot = FALSE)
          tmp$density <- tmp$density/max(tmp$density)
          over.hist[[ri,ci]] <- tmp
        }
        plot(over.hist[[ri,ci]], add = TRUE, freq = FALSE, col = col[ci],
             border = "white")
        if(mean.cred) {
          # mean and 95% credible intervals
          abline(v = c(over.mean[ri,ci], over.cred[,ri,ci]),
                 col = mean.cred.col, lty = c(1,2,2), lwd = 2)
        }
      } else {
        # diagonals: empty plots
        plot.window(xlim = c(0,1), ylim = c(0,1))
      }
      if(ri == 1 && ci == 1) {
        text(x = .5, y = .5,
             labels = expression("Niche Overlap: "*bgroup("(",
                                                          atop("Row", "Column"), ")")),
             adj = c(.5, .5), cex = 1)
      }
      box()
      if(ci != ri) axis(side = 1)
      if(ci > 1) axis(side = 2, labels = FALSE)
      if(ci < nspec) axis(side = 4, labels = FALSE)
      if(ri == 1) mtext(text = species.names[ci], side = 3, line = 1, col = col[ci])
      if(ci == 1) mtext(text = species.names[ri], side = 2, line = 1)
      if(mean.cred && ri == nspec && ci == nspec) {
        legend(x = "center", legend = c("Mean", "95% Credible Interval"),
               lty = c(1, 2), lwd = 2, col = mean.cred.col)
      }
    }
  }
  if(!missing(xlab)) {
    mtext(text = xlab, side = 1, line = 1.5, cex = 1, outer = TRUE)
  }
}
