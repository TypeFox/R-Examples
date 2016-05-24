.fonctiontsdMCMC <- function(data, x) {

# .tsd_fit <- function(par, males, N, temperatures, equation)
  
return(.tsd_fit(par=x, males=data$males, N=data$N, temperatures=data$temperatures, equation=data$equation))

}
