lin2 <-
function(x1,y1,x2,y2,x)
{
  y1 + (y2 - y1) * (x - x1) / (x2 - x1)
}
