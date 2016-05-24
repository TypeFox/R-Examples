fround <- function (x, digits) {
    format (round (x, digits), nsmall=digits)
}
  
pfround <- function (x, digits) {
    print (fround (x, digits), quote=FALSE)
}
 
