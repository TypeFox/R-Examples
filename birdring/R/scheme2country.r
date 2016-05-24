scheme2country <- function(x){
    # x: variable with scheme codes as given by EURING
  #filename <- system.file("data", "codes_scheme.txt", package = "birdring")
  #data(schemes, envir = environment())
  #codesdat <- read.table(filename, header=TRUE)
  schemes$Code <- gsub(" ", "", schemes$Code)
  namesx <- schemes$Country[match(as.character(x), as.character(schemes$Code))][drop=TRUE]
  return(namesx)
}