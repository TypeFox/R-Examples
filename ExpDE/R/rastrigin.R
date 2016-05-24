rastrigin <- function(Pop,...){
  apply(X = Pop,
        FUN = rastfun,
        MARGIN = 1)
}

rastfun <- function(x){
  n <- length(x);
  f <- 0;
  r <- 50;
  s <- 50;

  for (k in 1:n) {
    f <- f + (x[k] ^ 2 - 10 * cos(2 * pi * x[k]));
  }

  f1 <- 0;
  f2 <- 0;

  for (p in 1:n) {
    f1 <- f1 + max(0, (sin(2 * pi * x[p]) + 0.5)); #menor igual
    f2 <- f2 + abs(cos(2 * pi * x[p]) + 0.5); #igualdade
  }
  out <- f + r * f1 + s * f2;
  return (out);
}