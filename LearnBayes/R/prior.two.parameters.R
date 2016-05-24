prior.two.parameters = function(parameter1, parameter2) {
  prior = matrix(1, length(parameter1), length(parameter2))
  prior = prior/sum(prior)
  dimnames(prior)[[1]] = parameter1
  dimnames(prior)[[2]] = parameter2
  prior
 }