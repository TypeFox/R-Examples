coding.baseline.allele <-
  function(geno, coding = c("minor", "major")) {
    
    if (is.na(match(coding , c("minor", "major"))))
    { 
      stop( "coding should be minor or major.") 
    }
    tmptab <- rep(0,4)
    names(tmptab) <- c("0", "1", "2", "3")
    if (coding == "minor") {
      for (i in 1:dim(geno)[2]) {
        tt <- table(geno[,i])
        tmptab[names(tt)] <- as.numeric(tt)
        if (tmptab["1"] > tmptab["2"]) {
          # change "1" to "2"
          geno[, i] <- replace(geno[, i], geno[, i] == 1, 6 )
          geno[, i] <- replace(geno[, i], geno[, i] == 2, 1 )
          geno[, i] <- replace(geno[, i], geno[, i] == 6, 2 )
        }
        tmptab[1:4] <- rep(0, 4)
      }
    } # coding "minor"
    if (coding == "major") {
      for (i in 1:dim(geno)[2]) {
        tt <- table(geno[, i])
        tmptab [names(tt)] <- as.numeric(tt)
        if (tmptab["1"] < tmptab["2"]) {
          # change "1" to "2"
          geno[, i] <- replace(geno[, i] , geno[, i] == 1, 6)
          geno[, i] <- replace(geno[, i] , geno[, i] == 2, 1)
          geno[, i] <- replace(geno[, i] , geno[, i] == 6, 2)
        }
        tmptab[1:4] <- rep(0, 4)
      }
    } # coding "major"
    return ( geno )
  }
