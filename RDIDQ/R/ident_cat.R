ident_cat <-
function(test){
  cols <- 1:ncol(test)
categ_var <- sapply(cols,  function(i) is.character(test[, i]))
categ_var_test <- test[, categ_var]
return(categ_var_test)
}
