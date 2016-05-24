run1.moss <- function (filename, alpha = 1, c = 0.1, cPrime = 0.0001, q = 0.1, replicates = 5, maxVars = 3, dimens=NULL, k = 2) {


  if(missing(filename)) stop("Name of input file must be provided.")

  data <- read.table(filename, header=FALSE, sep="\t", stringsAsFactors=FALSE)
  if(ncol(data) <= 1)
  	data <- read.table(filename, header=FALSE, sep=" ", stringsAsFactors=FALSE)

  # If user did not provide the 'dimens' vector, we build it
  # assuming that all columns have the same number of values, 
  # except for the last column, which is representing CASE/Control.
  if(is.null(dimens)) {
    datacol <- dim(data)[2]
    dimens <- array(2, datacol)
    paramCounter <- apply(data, 2, unique)
    i <- 1
    #maxlen <- 2
    while(i<datacol) {
      if(length(unlist(paramCounter[i])) > 2)
        dimens[i] <- length(unlist(paramCounter[i]))
      i <- i + 1
    }
    #dimens <- array(maxlen, datacol)
    #dimens[datacol] <- 2  # last variable CASE/CONTROL
#print("**************************************")
#    print(dimens)
#print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    #print(maxlen)


  }
 
  return(MOSS.GWAS(alpha = alpha, c = c, cPrime = cPrime, q = q, replicates = replicates, maxVars = maxVars, data=data, dimens=dimens, k = k))



}
