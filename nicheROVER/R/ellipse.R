ellipse <-
function(mu, V, alpha = .95, n = 100) {
  tmp <- eigen(V)
  hlen <- sqrt(qchisq(alpha, df = 2)*tmp$val)
  theta <- atan2(tmp$vec[2,1], tmp$vec[1,1])
  t <- seq(0, 2*pi, len = n+1)
  x <- hlen[1] * cos(t)
  y <- hlen[2] * sin(t)
  alpha <- atan2(y, x)
  rad <- sqrt(x^2 + y^2)
  cbind(x = rad * cos(alpha + theta) + mu[1],
        y = rad * sin(alpha + theta) + mu[2])
}
