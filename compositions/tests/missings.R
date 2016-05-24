#library(compositions,lib.loc="/home/boogaart/R/compositions.Rcheck")


testvalues <- structure(c(1, 0, NA, NaN, -Inf, Inf, 1, -2), .Dim = c(2, 4), class = "acomp")

# testvalues <- matrix(c(TRUE,0.0,NA,NaN,-Inf,Inf,1,-2),nrow=2)
# class(testvalues) <-"acomp"

classi <- function(x) {
  ifelse( is.NMV(x) , "NMV" , ifelse( is.BDL(x) , "BDL" , ifelse( is.MAR(x), "MAR", ifelse( is.SZ(x) , "SZ", ifelse(is.MNAR(x),"MNAR","ERR")))))
}

missingset <- matrix(
c( 1, 1 , 1 ,1,
  NA, 1 , 1 ,1,
  NaN,1 , 1 ,1,
  -Inf,1,1 ,1,
  0.0 ,1,1 ,1,
  NA,NaN,1,1,
  NA,-Inf,1,1,
  NA,0.0,1,1,
  NaN,-Inf,1,1,
  NaN,0.0,1,1,
  -Inf,0.0,1,1,
  NA,-Inf,0.0,1,
  NA,NaN,0.0,1,
  NaN,NA,0.0,1,
  NaN,NA,0.0,Inf
  ),ncol=4,byrow=TRUE)

