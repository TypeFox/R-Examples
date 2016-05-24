gd <-
function(samp, libr, nuc, key) {
  n <- length(samp)
  if (n==1) {
    return(0)
  } else {
    gd <- matrix(0,n,n)
    for (i in 2:n) {
      vi <- which(key==samp[i])
      for (j in 1:(i-1)) {
        vj <- which(key==samp[j])
        if (samp[i]!=samp[j]) {
          if (!is.na(libr[[ vi ]][1])) {
            for (k in 1:length(libr[[ vi ]])) {
              if (!libr[[ vi ]][k] %in% libr[[ vj ]]) {
                gd[i,j] <- gd[i,j]+1
              } else if (libr[[ vi ]][k] %in% libr[[ vj ]]) {
                t <- which(libr[[ vj ]]==libr[[ vi ]][k])
                if (nuc[[ vi ]][k] != nuc[[ vj ]][t]) {
                  gd[i,j] <- gd[i,j]+1
                }
              }
            }
          }
          if (!is.na(libr[[ vj ]][1])) {
            for (k in 1:length(libr[[ vj ]])) {
              if (!libr[[ vj ]][k] %in% libr[[ vi ]]) {
                gd[i,j] <- gd[i,j]+1
              } else if (libr[[ vj ]][k] %in% libr[[ vi ]]) {
                t <- which(libr[[ vi ]]==libr[[ vj ]][k])
                if (nuc[[ vj ]][k] != nuc[[ vi ]][t]) {
                  gd[i,j] <- gd[i,j]+1
                }
              }
            }
          }
        }
        gd[j,i] <- gd[i,j]
      }
    }
    return(gd)
  }
}
