setColors <- function(colors){
  palette(colors)
  defaultColor <- colors[1]
  par(col = defaultColor,
      col.axis = defaultColor,
      col.lab = defaultColor,
      col.main = defaultColor,
      col.sub = defaultColor,
      fg = defaultColor)
  return(palette())
}

frColors <- function(dark = FALSE){
  if(dark)
    c("white", "plum1", "firebrick1", "sienna1", "yellow", "green",
      "cyan", "royalblue1", "sienna", "pink1", "gray69")
  else 
    c("black", "magenta2", "red2", "darkorange2", "gold2",
      "green3", "cyan2", "blue", "brown", "pink1", "grey85")
}
