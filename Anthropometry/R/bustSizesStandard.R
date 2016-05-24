bustSizesStandard <- function(bustCirc_4, bustCirc_6){
  bustCirc <- c(bustCirc_4, bustCirc_6)
  nsizes <- length(bustCirc)
  return(list(bustCirc=bustCirc, nsizes=nsizes))
}