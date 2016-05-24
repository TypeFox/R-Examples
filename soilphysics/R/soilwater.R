soilwater <-
function(x, theta_R, theta_S, alpha, n, m = 1-1/n,
   saturation.index = FALSE)
{
   sat.index <- (1 + (alpha * x) ^ n) ^ (-m)
   if (saturation.index) {
      out <- sat.index
      } else {
      out <- theta_R + (theta_S - theta_R) * sat.index
   }
   return(out)
}
