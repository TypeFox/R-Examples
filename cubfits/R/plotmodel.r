
### This plots one amino acid for ROC or NSE model.
plotmodel <- function(ret.model, main = NULL,
    xlab = "Production Rate (log10)", ylab = "Proportion",
    xlim = NULL, lty = 1, x.log10 = TRUE, ...){
  if(x.log10){
    ret.model$center <- log10(ret.model$center)
  }

  ### Observed dots and whiskers.
  if(is.null(xlim)){
    x.lim <- range(ret.model$center)
  } else{
    x.lim <- xlim
  }
  y.lim <- c(0, 1) + c(-0.05, 0.05)

  u.codon <- sort(colnames(ret.model)[-ncol(ret.model)])
  color <- get.color(u.codon)

  ### Reorder R for better legend.
  if(all(u.codon %in% .CF.GV$synonymous.codon$R)){
    u.codon <- u.codon[c(3:6, 1:2)]
    color <- color[c(3:6, 1:2)]
  }

  ### Make an empty plot
  plot(NULL, NULL, xlim = x.lim, ylim = y.lim,
       main = main, xlab = xlab, ylab = ylab, ...)

  ### Add modeled lines.
  plotaddmodel(ret.model, lty, u.codon, color, x.log10 = FALSE)

  ### Add focal codon.
  u.codon.star <- attr(ret.model, "u.codon.star")
  if(!is.null(u.codon.star)){
    if(all(u.codon %in% .CF.GV$synonymous.codon$R)){
      u.codon.star <- u.codon.star[c(3:6, 1:2)]
    }
  } else{
    u.codon.star <- u.codon
  }

  ### Add legends.
  legend(x.lim[1], y.lim[2], u.codon.star, col = color,
         box.lty = 0, lty = 1, pch = 19, cex = 0.8)
} # End of plotmodel().

plotaddmodel <- function(ret.model, lty, u.codon = NULL, color = NULL,
    x.log10 = TRUE){
  if(x.log10){
    ret.model$center <- log10(ret.model$center)
  }

  if(is.null(u.codon)){
    u.codon <- sort(colnames(ret.model)[-ncol(ret.model)])
  }
  if(is.null(color)){
    color <- get.color(u.codon)
  }
  
  if(is.list(ret.model)){
    linetype <- c(1, 2, 4, 3)
    for(i.pos in 1:length(ret.model)){
      for(i.codon in 1:length(u.codon)){
        lines(list(x = ret.model[[i.pos]]$center,
                   y = ret.model[[i.pos]][, colnames(ret.model[[i.pos]]) == u.codon[i.codon]]),
              col = color[i.codon],
              lty = linetype[i.pos], lwd = 1.2)
      }
    }
  }else{
    for(i.codon in 1:length(u.codon)){
      lines(list(x = ret.model$center,
                 y = ret.model[, colnames(ret.model) == u.codon[i.codon]]),
            col = color[i.codon],
            lty = lty, lwd = 1.2)
    }
  }
} # End of plotaddmodel().

