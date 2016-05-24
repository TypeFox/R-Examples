plotPanels <- function(kitname, data, xlim = NULL, cex = 0.75, alpha = 0.5){
  df <- data[[kitname]]
  df$marker <- as.character(df$marker)
  df[df$marker == "AMEL", "marker"] <- "A"
  dcoln <- unique(as.character(df$dye.col))
  ncol <- length(dcoln)
  bmin <- min(df$min.bp)
  bmax <- max(df$max.bp)
  if(is.null(xlim)) xlim <- c(bmin, bmax)
  plot.new()
  plot.window(xlim = xlim, ylim = c(0, ncol))
  yscale <- (ncol-1):0
  names(yscale) <- dcoln
  for(i in 1:nrow(df)){
    col <- as.character(df[i, "dye.col"])
    colalpha <- col2alpha(col, alpha)
    rect(df[i, "min.bp"], yscale[col] + 0.25, df[i, "max.bp"], yscale[col] + 0.5, col = colalpha)
    text(df[i, "min.bp"], yscale[col]+0.75, df[i, "marker"], pos = 4, cex = cex)
  }
  title(main = kitname, xlab = "Amplicon Size Ranges [bp]")
  axis(1)
}
