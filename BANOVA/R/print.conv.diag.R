print.conv.diag <-
function(x, ...){
  cat('Geweke Diag.\n')
  print(as.table(x$sol_geweke))
  cat('\n')
  cat("Heidelberger and Welch's Diag.\n")
  print(as.table(x$sol_heidel[,1:3]))
  if(x$pass_ind){
    cat('\n')
    cat("The Chain has converged.\n")
  }else{
    warning("The Chain may not have converged. Consider a longer burn-in, speeding up the convergence by setting conv_speedup = T, or modification of the model.\n", call. = F, immediate. = T)
  }
}
