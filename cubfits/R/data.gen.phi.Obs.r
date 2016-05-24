### Convert phi.df to phi.Obs

gen.phi.Obs <- function(phi.df){
  phi.Obs <- phi.df[, 2]
  names(phi.Obs) <- phi.df[, 1]
  phi.Obs
} # End of gen.phi.Obs().
