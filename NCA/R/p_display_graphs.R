p_display_graphs <-
function (results.nca, loop.data, title="NCA Plot", use.title=TRUE, pdf=FALSE, prefix="out") {
  create_outline <- function () {
    abline(h=loop.data$y.low,  lty=2, col="grey")
    abline(h=loop.data$y.high, lty=2, col="grey")
    abline(v=loop.data$x.low,  lty=2, col="grey")
    abline(v=loop.data$x.high, lty=2, col="grey")
  }

  lines <- list(ab=c("ols", "cols", "qr", "sfa", "cr_vrs", "cr_fdh"),
                lh=c("lh"),
                dea=c("ce_vrs", "ce_fdh"))
  
  # Open new window or PDF
  if (pdf) {
    p_new_pdf(prefix, loop.data$id.x, loop.data$id.y)
  } else {
    window.title <- paste(title, p_generate_title(loop.data), sep=" : ")
    p_new_window(title=window.title)
    par(family="")
    par(mfrow=c(1, 1))
  }
  
  # Plot the data points
  plot (loop.data$x, loop.data$y, 
        xlab=loop.data$names[1], ylab=loop.data$names[2], col="blue",
        xlim=c(loop.data$x.low, loop.data$x.high), 
        ylim=c(loop.data$y.low, loop.data$y.high))

  # Plot the scope outline
  create_outline()
  
  # Plot the lines, collect data for the legend
  legendParams = list()
  for (typeName in names(results.nca)) {
    type      <- results.nca[[typeName]]
    if (is.null(type$line)) {
      next
    }

    lineColor <- p_LINE_COLORS[[typeName]]
    lineType  <- p_LINE_TYPES[[typeName]]
    if (typeName %in% lines$ab) {
      abline(type$line, lty=lineType, col=lineColor, lwd=1.5)
    } else if (typeName %in% lines$lh) {
      lines (type$line[[1]], type$line[[2]], type="l",
             lty=lineType, col=lineColor, lwd=1.5)
    } else if (typeName %in% lines$dea) {
      dea.plot(type$line[[1]], type$line[[2]],
               RTS=unlist(strsplit(typeName, "_"))[2],
               ORIENTATION="graph", add=TRUE,
               lty=lineType, col=lineColor, lwd=1.5)
    }
    
    legendParams$names  = append(legendParams$names,  p_pretty_name(typeName))
    legendParams$types  = append(legendParams$types,  lineType)
    legendParams$colors = append(legendParams$colors, lineColor)
  }

  # Plot the legend
  if (length(legendParams) > 0) {
    legend("topleft", cex=0.7, legendParams$names,
           lty=legendParams$types, col=legendParams$colors)
  }

  # Plot the title
  if (use.title) {
    title <- paste0(title, " : ", p_generate_title(loop.data))
    title(title, cex.main=1)
  }
}
