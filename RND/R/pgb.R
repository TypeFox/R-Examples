pgb <-
function(x,a,b,v,w)
{
  x.new = (x/b)^a / (1 + (x/b)^a)
  out = pbeta(q = x.new, shape1 = v, shape2 = w, ncp = 0, lower.tail = TRUE)
  out
}
