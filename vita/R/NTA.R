NTA = function(PerVarImp)  UseMethod("NTA")

NTA.default = function (PerVarImp){

  # M1 = {VI_j|VI_j<0; j=0,...,p}
  M1 = PerVarImp[which(PerVarImp<0)]
  # M2 = {VI_j|VI_j=0; j=0,...,p}
  M2 = PerVarImp[which(PerVarImp==0)]
  # M3 = {-VI_j|VI_j<0; j=0,...,p}
  M3 = -M1
  # M = M1 u M2 u M3
  if(length(M2)==0){
      M = c(M1,M3)
  } else {
    M = c(M1,M2,M3)
  }
  # The empirical cumulative null distribution function Fn0 of M
  Fn0 = ecdf(M)
  out = list ( PerVarImp = PerVarImp,
              M = M,
              pvalue = matrix(1-Fn0(PerVarImp),ncol = 1,dimnames = list(dimnames(PerVarImp)[[1]] ,"p-value")  ),
              call = match.call())
  class(out) = "NTA"
  return(out)

}

print.NTA = function(x, ...){
  op = options()
  options(width = 90, digits = 2)
  if (!inherits(x, "NTA")) stop(" is not of class NTA")
  cat("Call:\n \n")
  print(x$call)
  cat("------------------------------------------------------------------------------------\n")
  cat("The non-positive variable importance values with the mirrored values:\n")
  cat("------------------------------------------------------------------------------------\n")
  print(format(t(x$M),digits = 2, nsmall = 1,justify = "left"),quote = F)
  cat("------------------------------------------------------------------------------------\n")
  cat("p-value :\n")
  cat("------------------------------------------------------------------------------------\n")
  print(format(t(x$pvalue),digits = 2, nsmall = 1,justify = "left"),quote = F)
  options(op)
}


