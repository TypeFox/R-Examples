kristen.colors <- function(n=64){

  ## Set up palette.
  col.ramp <- colorRampPalette(c(grey(.9),"darkorchid4", "blue","darkgreen","yellow","orange","red", "brown"))

  return(col.ramp(n))
}
