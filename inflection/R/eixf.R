eixf <-
function(x,y,f,i)
{
  p=0.5*(x[i + 1] - x[i]) * (y[i] - f(x[i]) + y[i + 1] - f(x[i + 1])) 
  p
}
