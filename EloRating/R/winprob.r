winprob <- function(elo1, elo2) {
  elo_diff <- (as.numeric(as.character(elo1)) - as.numeric(as.character(elo2)))
  # z score based on fixed SD=200 (see Elo 1978)   
  z_score <- elo_diff/(200*sqrt(2))  
  # calculates the winning probabilty
  p_win <- pnorm(z_score)
  return(p_win)
}
