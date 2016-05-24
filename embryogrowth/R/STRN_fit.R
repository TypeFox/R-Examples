.STRN_fit <- function(par, EmbryoGrowthTRN, tsd, Sexed, Males, Temperatures) {

  infoall <- info.nests(NestsResult=EmbryoGrowthTRN, SexualisationTRN=par, out="summary", replicate.CI = 1)

    temp_TSD <- infoall[, Temperatures]
    
#    par_TSD <<- par
    
  -sum(dbinom(prob=predict(tsd, temp_TSD, replicates=1)$sexratio, size=Sexed, x=Males, log=TRUE))
    
}
