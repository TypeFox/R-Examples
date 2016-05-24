p_display_bottleneck <-
function (bottlenecks, names, y.index, prefix, use.title, pdf, bottleneck.x, bottleneck.y, nx, steps) {
  # Put all tables in 1 file
  if (pdf) {
    p_new_pdf(prefix, "bottleneck", y.index, paper="A4r")
  }

  for (name in names(bottlenecks)) {
    p_display_table(bottlenecks[[name]], names, name, y.index, prefix, use.title, pdf, bottleneck.x, bottleneck.y, nx, steps)
  }
  
  if (pdf) {
    dev.off()
  }
}

p_display_table <-
function (mp, names, name, y.index, prefix, use.title, pdf, bottleneck.x, bottleneck.y, nx, steps) {
  # Set precision
  precision.x <- if (p_bottleneck_id(bottleneck.x) %in% c(1, 2, 4)) 1 else 3
  precision.y <- 3
  if (p_bottleneck_id(bottleneck.y) %in% c(1, 2, 4)) {
    precision.y <- if ((100 / steps) %% 1 == 0) 0 else 1
  }

  tmp <- matrix(sapply(mp, p_pretty_number, "NN", precision.x), nrow=steps+1)[,-1]
  if (nx == 1) {
    tmp <- matrix(tmp, ncol=1)
  }

  if (length(tmp) == 0) {
    return()
  }

  # Use numbered rows, legend as subtitle to the plot
  colnames(tmp) <- as.character(c(1:length(names[1:nx])))
  rownames(tmp) <- p_pretty_number(mp[,1], "", precision.y)
  
  start <- 1
  while (start < (steps+1)) {
    end <- min(start + 30, steps + 1)
    
    # A new window for each part
    if (!pdf) {
      dims <- p_calc_dims(names, nx, end-start)
      title <- paste("NCA Bottleneck table:", names[nx + y.index], p_pretty_name(name))
      p_new_window(title=title, width=dims[1], height=dims[2])
      par(family="")
    }

    textplot(if (nx == 1) tmp else tmp[start:end,],
             cex=1, halign="left", valign="top", mar=c(0, 0, 3, 0))

    title <- ""
    if (use.title) {
      title <- paste("Bottleneck", p_pretty_name(name), ":", names[nx + y.index])
    }
    legend <- paste0(bottleneck.y, " / ", bottleneck.x, "\n")
    for (i in seq(length(names[1:nx]))) {
      legend <- paste0(legend, i, " ", names[1:nx][i], "\n")
    }
    title(title, cex.main=1, sub=legend)
    
    start <- end
  }
}

p_calc_dims <-
function (names, nx, rows) {
  height <- 0.5 + 0.21 * (rows + length(names))
  
  width <- 1.5
  for (name in names[1:nx]) {
    width <- width + 0.8
  }
  
  return( c(width, height) )
}
