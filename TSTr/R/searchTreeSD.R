searchTreeSD <- function(keeper, string) {
  result <- as.logical(0)
  nodo <- 1
  i <- 1
  while (!is.na(nodo)) {
    if (str_sub(string,i,i) < keeper$ch[nodo]) {
      nodo <- keeper$L[nodo]
    }
    else if (str_sub(string,i,i) > keeper$ch[nodo]) {
      nodo <- keeper$R[nodo]
    }
    else{
      i <- i + 1
      if (i - 1 == nchar(string)) {
        result <- keeper$word[[nodo]]
      }
      else{
        nodo <- keeper$C[nodo]
      }
    }
  }
  result
}
