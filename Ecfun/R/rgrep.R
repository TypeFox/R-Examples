rgrep <- function(pattern, x, ignore.case = FALSE, 
      perl = FALSE, value = FALSE, fixed = FALSE, 
      useBytes = FALSE, invert = FALSE){
##
## 1, 2:  np & g.
##
  np <- length(pattern)
  g. <- rep(NA, np)
##
## 3.  for each pattern
##
  for(i in seq(length=np)){
    gi <- grep(pattern[i], x, ignore.case = ignore.case, 
        perl = perl, value = value, fixed = fixed, 
               useBytes = useBytes, invert = invert)
    g.[i] <- (length(gi)>0) 
  }
##
## 3.  done 
##
  if(value){
    return(pattern[g.])
  } else return(which(g.))
}
