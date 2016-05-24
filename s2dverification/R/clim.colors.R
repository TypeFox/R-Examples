clim.colors <- function(n) {
  colorbar <- colorRampPalette(c("dodgerblue4", "dodgerblue1", "forestgreen", 
                                 "yellowgreen", "white", "white", "yellow", 
                                 "orange", "red", "saddlebrown"))
  colorbar(n)
}
