atCfn <-
function(aM,
                  M,
                  At,
                  at)
{
  #Estimate atC
  atC <- mapply(function(a, b){c(t(t(aM) + M %*% inv(a + M) %*% (b - t(aM))))},
                At, at, SIMPLIFY = FALSE)
  #Input names
  atC <- lapply(atC, function(a) {names(a) <- colnames(aM); a})
  #Output
  atC
}
