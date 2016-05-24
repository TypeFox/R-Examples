risk.colors <- function(n)
{
    j <- n%/%4
    c(rgb(0.86*(1:j)/j, 0, 0), 
      rainbow(n-2*j, start=1/20,end=1/7),
      hsv(h = 1/6,s = seq(from = 1 - 1/(2 * j), to = 1/(2 * j), length = j), v = 1))[n:1]
}

