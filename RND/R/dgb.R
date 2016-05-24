dgb <-
function(x,a,b,v,w)
{
  x.new = (x/b)^a / (1 + (x/b)^a)
  out   = dbeta(x.new, shape1 = v, shape2 = w, ncp = 0, log = FALSE) * (a * b^a) * (x^(a-1)) / (x^a + b^a)^2
  out
}
