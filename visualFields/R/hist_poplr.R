hist_poplr <- function( scomb_obs, pcomb_obs, scomb, txtfont = "mono", pointsize = 7 ) {

  limExtention <- 2
  sizeBreaks   <- 4
  breaks       <- c( 0:ceil( max( scomb[scomb != Inf] ) / sizeBreaks ) * sizeBreaks )
  hcomb <- hist( scomb, breaks = breaks, plot = FALSE )
  hcomb$counts <- hcomb$counts / max( hcomb$counts )
  ops     <- par()$ps
  ofamily <- par()$family
  oplt    <- par()$plt
  par( ps     = pointsize )
  par( family = txtfont )
  par( plt    = c( 0, 1, 0.3, 1 ) )
  plot( hcomb, main = "", xlab = NULL, ylab = NULL, xlim = c(0, limExtention * max( scomb[scomb != Inf] ) ), ylim = c( 0, 1.1 ), border = rgb( 0.7, 0.7, 0.7 ), col = rgb( 0.9, 0.9, 0.9 ), axes = FALSE )
  axis( 1, las = 1, tcl = -.3, lwd = 0.5, lwd.ticks = 0.5 )
  title( xlab = "S statistic", mgp = c( 2, 1, 0 ) )
  lines( c( scomb_obs, scomb_obs ), c( 0, 0.60 ), col = "red" )
  points( c( scomb_obs ), c( 0.60 ), pch = 22, col = "red", bg = "red" )
  scomb_txt <- as.character( round( scomb_obs ) )
  pcomb_txt <- as.character( round( pcomb_obs * 1000 ) / 1000 )
  text( c( scomb_obs ), c( 0.90 ), labels = paste( "S =", scomb_txt ) )
  text( c( scomb_obs ), c( 0.75 ), labels = paste( "p =", pcomb_txt ) )
  par( ps     = ops )
  par( family = ofamily )
  par( plt    = oplt )
}