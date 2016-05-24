propCI_wilson <-
function(x, n, l){

  p <- x/n
  c <- l[2] - l[1]

  if (x == 0){
    z <- qnorm(c)
    z_sq <- z ^ 2
    ci <- c(0,
           (p +
              z_sq / (2 * n) +
              z * sqrt(p * (1 - p) / n + z_sq / (4 * n^2))) /
            (1 + z_sq / n))
  } else if (x == n){
    z <- qnorm(c)
    z_sq <- z ^ 2
    ci <- c((p +
               z_sq / (2 * n) -
               z * sqrt(p * (1 - p) / n + z_sq / (4 * n^2))) /
            (1 + z_sq / n),
            1)
  } else{
    z <- qnorm(l[1])
    z_sq <- z ^ 2
    ci <- (p +
             z_sq / (2 * n) +
             c(1, -1) * z * sqrt(p * (1 - p) / n + z_sq / (4 * n^2))) /
          (1 + z_sq / n)
  }

  return(ci)
}
