dmatrix(x, model) %::% a : WignerModel : a
dmatrix(x, model) %as% {
  sqrt(4 - x^2) / (2 * pi)
}

dmatrix(x, model) %::% a : WishartModel : a
dmatrix(x, model) %as% {
  var <- model$sd^2
  bounds <- domain(model)
  b.neg <- bounds[1]
  b.pos <- bounds[2]
  ind <- ifelse(b.neg < x & x < b.pos, 1, 0)
  ind * model$Q / (2*pi*var*x) * sqrt((x - b.neg) * (b.pos - x))
}

dmatrix(x, model) %::% a : JacobiModel : a
dmatrix(x, model) %as% {
  flog.warn("This function is incomplete")
  c1 <- model$n / model$m1
  c2 <- model$n / model$m2
  c0 <- c1*x + x^3*c1 - 2*c1*x^2 - c2*x^3 + c2*x^2

  b0 <- c1*x - c2*x - c1 + 2
  b1 <- -2*c2*x^2 + 2*x - 3*c1*x + c1 + c2*x - 1 + 2*c1*x^2
  b2 <- c0

  #num <- (1 - c1)^2
  #b.neg <- num / (c1^2 - c1 + 2 + c2 - c1 * c2 + 2 * sqrt(c1 + c2 - c1 * c2))
  #b.pos <- num / (c1^2 - c1 + 2 + c2 - c1 * c2 - 2 * sqrt(c1 + c2 - c1 * c2))

  #ind <- ifelse(b.neg < x & x < b.pos, 1, 0)
  #sqrt(ind * (x - b.neg) * (x + b.neg)) / (2 * pi * c0)
  sqrt(4*b2*b0 - b1^2) / (2*pi*b2)
}

# Get the bounds of the eigenvalues for the given model
domain(model) %::% WishartModel : numeric
domain(model) %when% {
  model %hasa% sd
  model %hasa% Q
} %as% {
  domain.min <- model$sd^2 * (1 - sqrt(1/model$Q))^2
  domain.max <- model$sd^2 * (1 + sqrt(1/model$Q))^2
  c(domain.min, domain.max)
}


