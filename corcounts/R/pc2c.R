`pc2c` <-
function(Theta) {
  n <- dim(Theta)[2]
  Theta.alt <- Theta
  for (k in 2:(n-1)) {      # Anfangszeile
    for (i in (n-1):k) {    # Zeile (beginne unten)
      for (j in n:(i+1)) {  # Spalte (beginne rechts)
        Theta[i,j] <- Theta[i,j]*sqrt((1-Theta.alt[i-k+1,i]^2)*(1-Theta.alt[i-k+1,j]^2)) +
                      Theta.alt[i-k+1,i]*Theta.alt[i-k+1,j]
      }
    }
  }

  C <- diag(rep(1,n))
  for (i in 2:n) {
    for (j in 1:(i-1)) {
      C[i,j] <- Theta[j,i]
      C[j,i] <- Theta[j,i]
    }
  }
  return(C)
}

