prob.table <-
function(mod){
  if("omitted" %in% names(mod)){
    warning("The co-occurrence model was run using 'thresh = TRUE.' The probability table may not include all species pairs")
  }
  mod$results
}
