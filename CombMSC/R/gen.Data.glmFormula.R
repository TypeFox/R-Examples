`gen.Data.glmFormula` <-
function(fmla, pars, var.Frame, data.Sigma=1, data.Size, train.Frame, test.Size, ...){
  if(dim(train.Frame)[1] > dim(var.Frame)[1]) stop ("Screwed up training frame / var.Frame")
  Terms <- delete.response(terms(fmla))
  m <- model.frame(Terms, var.Frame)
  X <- model.matrix(Terms, m)
  temp <- X %*% pars
  prob.Vector <- exp(temp)/(1+exp(temp))
  temp <- rbinom(n = length(prob.Vector), size = 1, prob = prob.Vector)
  train.Size <- dim(train.Frame)[1]
  train.Inds = 1:train.Size
  test.Inds = (train.Size + 1):dim(var.Frame)[1]
  if(dim(train.Frame)[1] == dim(var.Frame)[1]) test.Inds = NULL
  list(data.Vector = temp[train.Inds], test.Vector = temp[test.Inds])
}

