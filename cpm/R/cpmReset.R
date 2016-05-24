cpmReset <- function(cpm) {
    #Lepage is handed separately since both its component CPMS must
    #be reset
    if (class(cpm) == "ChangePointModelLepage") {
        return(cpmResetLepage(cpm))
    }
    
    cpm@n <- 0
    cpm@changeDetected <- FALSE
    for (i in 1:length(cpm@windowStatistic)) {
        cpm@windowStatistic[[i]] <- numeric()
    }
    
  return(cpm)
}