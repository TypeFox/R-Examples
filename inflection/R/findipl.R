findipl <-
function(x,y,j)
{
  n=dim(x)[1]
  x1=x[1]
  y1=y[1]
  x2=x[j]
  y2=y[j]
  flin1=function(x)
  {
    lin2(x1,y1,x2,y2,x)
  }
  sl=0
  for (i in 1:(j-1))
  {
    sl=sl+eixf(x,y,flin1,i)
  }
  x1=x[j]
  y1=y[j]
  x2=x[n]
  y2=y[n]
  flin2=function(x)
  {
    lin2(x1,y1,x2,y2,x)
  }
  sr=0
  for (k in j:(n-1))
  {
    sr=sr+eixf(x,y,flin2,k)
  }
  c(j,x[j],sl,sr)
}
