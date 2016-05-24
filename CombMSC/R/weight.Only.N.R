`weight.Only.N` <-
function(model.List, n){
  temp <- sapply(model.List, num.Terms)
  ifelse(temp==n, 1, 0)
}

