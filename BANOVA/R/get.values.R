get.values <-
function(data, mf, id_index = array(dim = 0), intercept_names = character(0)){

  var_names <- c(intercept_names, attr(mf,'names'))
  results <- list()
  if (length(var_names) > 0)
    for (i in 1: length(var_names)){
      if (length(id_index) == 0)
        eval(parse(text = paste('results[[',i,']] <- data$',var_names[i],sep='')))
      else
        eval(parse(text = paste('results[[',i,']] <- data$',var_names[i], '[id_index]', sep='')))
    }
  return(results)
}
