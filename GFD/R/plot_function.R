plotting <- function(plot.object, descr.object, factor, col, pch, legendpos, ...){

  nf <- plot.object$nf
  color <- col
  # one-way 
  if(nf == 1){
    xmax <- length(plot.object$levels[[1]])
    y <- descr.object[, 3]
    li <- descr.object[, 5]
    ui <- descr.object[, 6]
    
    plotrix::plotCI(x = 1:xmax, y, li = li,
                    ui = ui, xlim = c(0.8, length(plot.object$levels[[1]]) + 0.3), xaxt = "n", col = color[1], pch = pch[1], ...)
    axis(side = 1, at = 1:1:length(plot.object$levels[[1]]), labels = plot.object$levels[[1]], ...)
  } else {
    
    TYPE <- plot.object$Type    
    Faktor <- factor
    nadat2 <- plot.object$nadat2
    levels <- plot.object$levels
    mu <- plot.object$mu
    lower <- plot.object$lower
    upper <- plot.object$upper
    lev_names <- plot.object$lev_names
    mean_out <- descr.object$Means
    CI <- matrix(c(descr.object$Lower, descr.object$Upper), ncol = 2)
    colnames(CI) <- c("CIl", "CIu")
    fac_names <- plot.object$fac_names
    fl <- plot.object$fl
    fac_names_original <- plot.object$fac_names_original

    # plots for interactions
    if (TYPE == "nested") {
      # main effect
      if (Faktor == nadat2[1]) {
        plotrix::plotCI(x = 1:length(levels[[1]]), mu[[1]],
                        li = lower[[1]], ui = upper[[1]], xlim = c(0.8, length(levels[[1]]) + 0.3),
                        col = color[1], pch = pch[1], xaxt = "n", ...)
        axis(side = 1, at = 1:1:length(levels[[1]]), labels = levels[[1]], ...)
      } else if (Faktor == fac_names[2] && nf == 2) {
        plotrix::plotCI(x = 1:length(levels[[2]]), mean_out,
                        li = CI[, 1], ui = CI[, 2], xlim = c(0.8, length(levels[[2]]) + 0.3),
                        ylim = c(min(CI) - 1, max(CI) + 1), col = color[1], pch = pch[1], xaxt = "n", ...)
        axis(side = 1, at = 1:1:length(levels[[2]]), labels = lev_names[, 2], ...)
        aa <- length(levels[[2]]) / fl[1] ^ 2
        bb <- length(levels[[2]]) / fl[1]
        cc <- length(levels[[2]]) + aa
        ss <- seq(from = - aa, to = cc, by = bb)
        axis(side = 3, at = ss[2:(length(ss) - 1)], labels = levels[[1]], ...)
      } else if(Faktor %in% fac_names && nf == 3 && !(Faktor == nadat2[1])) {
        error <- "For three-way nested design, only the main effect
    can be plotted."
      }}
    
    
    if (TYPE == "crossed") {
      # plot of main effects
      for (i in 1:nf) {
        if (Faktor == nadat2[i]) {
          plotrix::plotCI(x = 1:length(levels[[i]]), mu[[i]],
                          li = lower[[i]], ui = upper[[i]], xlim = c(0.8, length(levels[[i]]) + 0.3),
                          col = color[1], pch = pch[1], xaxt = "n", ...)
          axis(side = 1, at = 1:1:length(levels[[i]]), labels = levels[[i]], ...)
        }}
      
      # two-fold interactions for three- and higher-way layout
      fac_names_twofold <- plot.object$fac_names_original[ - (1:nf)]
      fac_names_twofold <- fac_names_twofold[1:choose(nf, 2)]
      
      if (Faktor %in% fac_names_twofold) {
        nmu <- list()
        nsigma <- list()
        nn_groups <- list()
        nupper <- list()
        nlower <- list()
        new_levels <- list()
        dat2 <- plot.object$dat2
        fl <- plot.object$fl
        alpha <- plot.object$alpha
        
        counter <- 1
        for (i in 2:nf) {
          for (j in (i + 1):(nf + 1)) {
            nmu[[counter]] <- matrix(by(dat2[, 1], dat2[, c(i, j)], mean),
                                     nrow = fl[i - 1])
            nsigma[[counter]] <- matrix(by(dat2[, 1], dat2[, c(i, j)], var),
                                        nrow = fl[i - 1])
            nn_groups[[counter]] <- matrix(by(dat2[, 1], dat2[, c(i, j)],
                                              length), nrow = fl[i - 1])
            nlower[[counter]] <- nmu[[counter]] -
              sqrt(nsigma[[counter]] / nn_groups[[counter]]) *
              qt(1 - alpha / 2, df = nn_groups[[counter]])
            nupper[[counter]] <- nmu[[counter]] +
              sqrt(nsigma[[counter]] / nn_groups[[counter]]) *
              qt(1 - alpha / 2, df = nn_groups[[counter]])
            new_levels[[counter]] <- list(levels[[i - 1]], levels[[j - 1]])
            counter <- counter + 1
          }
        }
        names(nmu) <- fac_names_twofold
        names(nupper) <- fac_names_twofold
        names(nlower) <- fac_names_twofold
        place <- which(Faktor == fac_names_twofold)
        xxx <- rep(NA, nf)
        for (ii in 1:nf) {
          xxx[ii] <- grepl(plot.object$fac_names_original[ii], Faktor, fixed = TRUE)
        }
        plotrix::plotCI(x = 1:length(new_levels[[place]][[2]]),
                        nmu[[Faktor]][1, ],
                        li = nlower[[Faktor]][1, ],
                        ui = nupper[[Faktor]][1, ], xlim = c(0.8, length(new_levels[[place]][[2]]) + 0.3),
                        ylim = c(min(nlower[[Faktor]]) - 1, max(nupper[[Faktor]]) + 1),
                        col = color[1], pch = pch[1], xaxt = "n", ...)
        axis(side = 1, at = 1:1:length(new_levels[[place]][[2]]),
             labels = new_levels[[place]][[2]], ...)
        for (i in 2:length(new_levels[[place]][[1]])) {
            plotrix::plotCI(x = ((1:length(new_levels[[place]][[2]])) + 0.07 * i),
                          nmu[[Faktor]][i, ], li = nlower[[Faktor]][i, ],
                          ui = nupper[[Faktor]][i, ], add = TRUE, col = color[i], pch = pch[1], ...)
            
        }
        legend(legendpos, 
               legend = new_levels[[place]][[1]],
               col = color[1:length(new_levels[[place]][[1]])],
               seg.len = 0.5, 
               lty = rep(1, length(new_levels[[place]][[1]])))
      } else if (nf == 3 && Faktor == fac_names[length(fac_names)]) {
        # three-way
        Var1 = Var2 = Var3 = CIl = CIu = NULL
        for_plots <- cbind(lev_names, mean_out, CI)
        group <- list()
        for (i in 1:length(levels[[1]])) {
          group[[i]] <- subset(for_plots, Var1 == as.factor(levels[[1]])[i],
                               select = c(Var2, Var3, mean_out, CIl, CIu))
        }
        next_group <- list()
        new_group <- list()
        for (j in 1:length(group)) {
          for (l in 1:length(levels[[2]])) {
            next_group[[l]] <- subset(group[[j]],
                                      Var2 == as.factor(levels[[2]])[l],
                                      select = c(mean_out, CIl, CIu))
          }
          new_group[[j]] <- next_group
        }
        counter <- 1
        delta <- seq(from = 0, by = 0.05,
                     length = length(levels[[1]]) * length(levels[[2]]) + 1)
        plotrix::plotCI(x = 1:length(levels[[3]]), new_group[[1]][[1]][, 1],
                        li = new_group[[1]][[1]][, 2],
                        ui = new_group[[1]][[1]][, 3], ylim = c(min(CI) - 1, max(CI) + 1),
                        xlim = c(0.8, length(levels[[3]]) + 0.3),
                        col = color[1], pch = pch[1], xaxt = "n", ...)
        axis(side = 1, at = 1:1:length(levels[[3]]), labels = levels[[3]], ...)
        for (j in 1:length(levels[[1]])) {
          for (i in 1:length(levels[[2]])) {
            plotrix::plotCI(x = ((1:length(levels[[3]])) + delta[counter]),
                            new_group[[j]][[i]][, 1], li = new_group[[j]][[i]][, 2],
                            ui = new_group[[j]][[i]][, 3], add = TRUE, col = color[i], pch = pch[j], ...)
            counter <- counter + 1
          }}
        legend(legendpos, 
               legend = c(fac_names_original[[1]], levels[[1]], fac_names_original[[2]], levels[[2]]), box.lty = 0,
               col = c(0, rep(1, length(levels[[1]])), 0, color[1:length(levels[[2]])]),
               pch = c(NA, pch[1:length(levels[[1]])], NA, rep(NA, length(levels[[2]]))),
               seg.len = 0.5,
               lty = c(NA, rep(NA, length(levels[[1]])), NA, rep(1, length(levels[[2]]))))
      } else if (Faktor %in% fac_names && nf >= 4) {
        error <- "Higher-way interactions cannot be plotted!"
      }
    }   
  }
}
