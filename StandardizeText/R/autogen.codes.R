autogen.codes <-
function(n) {
  if (n > 456976) {
    cat("\nYou have requested more than 26^4 codes. Please specify a manual coding scheme.\n")
    return(NULL)
  }
  generated.codes <- expand.grid(LETTERS[-24],LETTERS[-24])
  generated.codes <- paste(generated.codes[,2],generated.codes[,1],sep="")
  while (n > length(generated.codes)) {
    generated.codes <- expand.grid(LETTERS[-24],generated.codes)
    generated.codes <- paste(generated.codes[,2],generated.codes[,1],sep="")
  }
  return(generated.codes[seq(n)])
}
