is.error <- function(Zt, Zp) {
  error <- array(data = 0, dim = length(Zp))

  if (length(Zt) != 0) {
    zt <- as.numeric(levels(Zt))
    zp <- as.numeric(levels(Zp))

    for (i in 1:length(zt)) {
      Zti <- Zt[which(Zt == zt[i])]
      Zpi <- Zp[which(Zt == zt[i])]

      j <- as.numeric(names(sort(table(Zpi), decreasing = TRUE)))

      k <- which(j %in% zp == TRUE)

      if (length(k) == 0) {
        error[which(Zt == zt[i])] = 1
      }
      else {
        j <- j[k[1]]
    
        error[which(Zt == zt[i] & Zp != j)] = 1

        zp <- zp[which(zp != j)] 
      }
    }
  }
  
  error
} ## is.error
