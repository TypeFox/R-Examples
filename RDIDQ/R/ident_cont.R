ident_cont <-
function(test){
  cols <- 1:ncol(test)
  cont_var <- sapply(cols,  function(i) is.numeric(test[, i]))
  cont_var_test <- test[, cont_var]
  return(cont_var_test)
}
