meansnps <-
function(strain.log, freq.log, libr, nuc, key) {
  types <- length(strain.log) # number of unique types
  exp <- 0
  if (types!=1) {
    for (i in 2:types) {
      str1 <- which(key==strain.log[[i]])
      for (j in 1:(i-1)) {
        str2 <- which(key==strain.log[[j]])
        exp <- exp + (freq.log[i]/sum(freq.log))*(freq.log[j]/sum(freq.log))*
                      (sum(!libr[[str1]] %in% libr[[str2]]) + 
                         sum(!libr[[str2]] %in% libr[[str1]]))
      }
    }
  }
  return(exp)
}
