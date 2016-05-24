## everything for binaryRep except the function itself
setClass("binaryRep",
         representation( original = "numeric",
                        sign = "integer",
               exponent = "integer",
               bits = "raw"))

setMethod("show", "binaryRep",
          function(object){
              cat("Object of class \"", class(object), "\"\n", sep="")
              sign <- ifelse(object@sign < 0, "-", " ")
              pasteBits <- function(bits)
                paste(ifelse(bits>0, "1","0"),collapse="")
              origSep <- "\n      " # could be " " if options()$width large
              bits <- matrix(object@bits, nrow = length(sign))
              bits <- apply(bits, 1, pasteBits)
              lines <- paste(format(1:length(sign)),": ",sign,".",
                             bits, " * 2^", object@exponent, origSep,
                             "(", object@original, ")", sep="")
              cat(lines, sep="\n")
          })

binaryRepA <- function(m)
    matrix(rep(c(1, rep(0, m-1), -1), length=m*m), m, m)

binaryRepPowers <- function(n, m)
    (matrix(2 ^ (1:m), 1, m, byrow=TRUE))[rep(1,n),]

binaryRepBits <- function(data) {
    data <- as.vector(data) # throw away matrix structure
    if(any((data%%1)^2 > .Machine$double.eps))
        stop("The \"bits\" have fractional parts.")
    data <- round(data)
    if(!identical(sort(unique(data)), c(0,1)))
        stop("The \"bits\" must be 0 or 1.")
    as.raw(data)
}

