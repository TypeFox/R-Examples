"applyWeightingModel" <- 
function (model) 
{
  if (length(model@weightpar) > 0 || length(model@weightList) > 0) {
    w <- weightPsi(model)
    model@weight <- TRUE
    model@psi.weight <- w$psi.weight
    model@weightM <- w$weight
  }
  else model@psi.weight <- model@psi.df
  
  model
}
