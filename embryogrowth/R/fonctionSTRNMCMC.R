.fonctionSTRNMCMC <- function(data, x) {

# .STRN_fit <- function(par, EmbryoGrowthTRN, tsd, Sexed, Males, Temperatures) {
  
return(.STRN_fit(par=x, EmbryoGrowthTRN=data$EmbryoGrowthTRN, tsd=data$tsd, 
                Sexed=data$Sexed, Males=data$Males, 
                Temperatures=data$Temperatures))

}
