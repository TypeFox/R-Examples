avedistTwo <-
function(xb1, xb2) {

  distances1 <- apply(xb1, 1, function(y) sqrt(rowSums(sweep(xb2, 2, y, "-")^2)))
  distances2 <- apply(xb2, 1, function(y) sqrt(rowSums(sweep(xb1, 2, y, "-")^2)))

  mean(c(apply(distances1, 2, min), apply(distances2, 2, min)))

}
