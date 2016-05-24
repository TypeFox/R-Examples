approximate.max <-
function(x,y, k = 5)
{
  out = ( 1/(1 + exp(-k*(x-y))) )*x + ( 1 - 1/(1 + exp(-k*(x-y))) ) * y
  out
}
