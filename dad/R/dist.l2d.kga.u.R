dist.l2d.kga.u <-
function(x1, x2, check = FALSE)  {
  # x1, x2 :       samples
  return(sqrt(l2d.kga.u(x1, x1, check = check) + l2d.kga.u(x2, x2, check = check) - 2*l2d.kga.u(x1, x2, check = check)))
}
