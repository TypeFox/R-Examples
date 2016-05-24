munsell_to_rgb <-
function(color, name) {
  color <- as.character(color)
  separated <- unlist(strsplit(color, ' '))
  hue <- separated[1]
  remaining <- separated[2]
  separated2 <- unlist(strsplit(remaining, '/'))
  value <- round(as.numeric(separated2[1]))
  chroma <- round(as.numeric(separated2[2]))
  col <- suppressWarnings(munsell2rgb(hue, value, chroma))
  if (length(col)!=1) {
        col <- '#39302CFF'
    warning(paste('munsell color in',name,'missing or not valid, replaced with 5YR 2/1'))
  }
  invisible(col)
}
