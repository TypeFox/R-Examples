# Internal function of MRtest
ProcTest <- function(Teste) {
  if (is.matrix(Teste) == FALSE) {
    g  <- diag(length(Teste))
    Teste <- cbind(Teste, g)
  }
  return(Teste)
}
