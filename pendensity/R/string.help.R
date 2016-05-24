string.help <- function(string, star = " ") {
  str.len <- nchar(star)
  str.loc <- NULL
  for (i in 1:nchar(string)) {
    if (substring(string,i,(i+str.len-1)) == star) {
      str.loc <- c(str.loc, i)
    }
  }
  no <- length(str.loc)+1
  output <- rep(0,no)
  if(no == 1)
    output[1] <- string
  else {
    output[1] <- substring(string,1,(str.loc[1] - 1))
    if (no > 2) {
      for(i in 2:(no-1)) {
        output[i] <- substring(string,(str.loc[i-1]+str.len), (str.loc[i]-1))
        }
      }
  output[no] <- substring(string, (str.loc[no - 1] +str.len), nchar(string))
  }
  return(output)
}
