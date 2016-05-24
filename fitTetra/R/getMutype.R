getMutype <-
function(model) {
  rep(c(1,2,5,6),2)[((model-1)%%8)+1]
}
