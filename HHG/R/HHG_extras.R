hhg.example.datagen = function(n, example) {
  if (example == '') {
  } else if (example == '4indclouds') {
    .datagen4indclouds(n)
  } else if (example == '2Parabolas') {
    .datagen2Parabolas(n)
  } else if (example == 'W') {
    .datagenW(n)
  } else if (example == 'Parabola') {
    .datagenParabola(n)
  } else if (example == 'Diamond') {
    .datagenDiamond(n)
  } else if (example == 'Circle') {
    .datagenCircle(n)
  } else if (example == 'TwoClassUniv') {
    .datagenTwoClassUniv(n)
  } else if (example == 'FourClassUniv') {
    .datagenFourClassUniv(n)
  } else if (example == 'TwoClassMultiv') {
    .datagenTwoClassMultiv(n)
  } else {
    stop('Unexpected example specified. Please consult the documentation.')
  }
}

.datagen4indclouds = function(n) {
  dx = rnorm(n) / 3
  dy = rnorm(n) / 3
  cx = sample(c(-1, 1), size = n, replace = T)
  cy = sample(c(-1, 1), size = n, replace = T)
  u = cx + dx
  v = cy + dy 
  return (rbind(u, v))
}

.datagen2Parabolas = function(n) {
  x = seq(-1, 1, length = n)
  y = (x ^ 2 + runif(n) / 2) * (sample(c(-1, 1), size = n, replace = T))
  return (rbind(x, y))
}

.datagenW = function(n) {
  x = seq(-1, 1, length = n)
  u = x + runif(n)/3
  v =  4*( ( x^2 - 1/2 )^2 + runif(n)/500 )
  return (rbind(u,v))
}

.datagenParabola = function(n) {
  x = seq(-1, 1, length = n)
  y = (x ^ 2 + runif(n)) / 2
  return (rbind(x,y))
}

.datagenDiamond = function(n) {
  x = runif(n, min = -1, max = 1)
  y = runif(n, min = -1, max = 1)

  theta = -pi / 4
  rr = rbind(c(cos(theta), -sin(theta)),
             c(sin(theta),  cos(theta)))
  tmp = cbind(x, y) %*% rr
  u = tmp[,1]
  v =  tmp[,2]
  return (rbind(u, v))
}

.datagenCircle = function(n) {
  x = seq(-1, 1, length = n)
  u = sin(x * pi) + rnorm(n) / 8
  v = cos(x * pi) + rnorm(n) / 8
  return (rbind(u, v))
}

.datagenTwoClassUniv = function(n) {
  y = as.double(runif(n) < 0.5)
  x = y * rnorm(n, mean = -0.2) + (1 - y) * rnorm(n, mean = 0.2)
  return (list(x = x, y = y))
}

.datagenFourClassUniv = function(n) {
  y = as.double(sample(x = 0:3, size = n, replace = T))
  x = (y == 1) * rnorm(n, mean = -0.4) + 
      (y == 2) * rnorm(n, mean = -0.2) +
      (y == 3) * rnorm(n, mean =  0.2) +
      (y == 4) * rnorm(n, mean =  0.4)
  return (list(x = x, y = y))
}

.datagenTwoClassMultiv = function(n) {
  m = 10
  x = matrix(as.double((runif(n * m) < 0.4) + (runif(n * m) < 0.4)), ncol = m)
  y = as.double(xor(rowSums(x[, 1:5] > 0) > 2, rowSums(x[, 6:10] > 0) > 2))
  return (list(x = x, y = y))
}








