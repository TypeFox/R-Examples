`klin.ls` <-
function(A,b) {
  if (class(A)!="klin.prepls")
    A <- klin.preparels(A)
  klin.solve(A$lhs,klin.eval(A$rhs,b,transpose=TRUE))
}

