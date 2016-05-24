RAM.pal <- function(cols.needed=20) {
  # color palette total of 40
  col.pal <- c(RColorBrewer::brewer.pal(12, "Set3"), RColorBrewer::brewer.pal(8, "Accent"), .ram.pal(20))
    
  if (cols.needed <= 40  ) {
    # need to manually construct a palette to use
    col.func <- function(n){col.pal[1:n]}
    values=col.func(cols.needed)
  } else {
    col.func  <-  grDevices::colorRampPalette(col.pal)
    values=col.func(cols.needed)
  }
  return(values)
}
