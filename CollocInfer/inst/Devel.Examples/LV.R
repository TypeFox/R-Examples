make.LV = function()
{
LV = function(t,y,p,more=NULL)
{
  x = y
  x[,1] = p[1]*y[,1]-p[2]*y[,1]*y[,2]
  x[,2] = -p[3]*y[,2]+p[4]*y[,1]*y[,2]
  return(x)
}

LVdx =  function(t,y,p,more=NULL)
{
  x = array(0,c(dim(y),dim(y)[2]))
  x[,1,1] = p[1] - p[2]*y[,2]
  x[,1,2] = -p[2]*y[,1]
  x[,2,1] = p[4]*y[,2]
  x[,2,2] = -p[3]+p[4]*y[,2]
  return(x)
}

LVdp =  function(t,y,p,more=NULL)
{
  x = array(0,c(dim(y),length(p)))
  x[,1,1] = y[,1]
  x[,1,2] = -y[,1]*y[,2]
  x[,2,3] = -y[,2]
  x[,2,4] = y[,1]*y[,2]
  return(x)
}

LVdx2 = function(t,y,p,more=NULL)
{
  x = array(0,c(dim(y),dim(y)[2],dim(y)[2]))
  x[,1,1,2] = -p[2]
  x[,1,2,1] = -p[2]
  x[,2,1,2] = p[4]
  x[,2,2,2] = p[4]
  return(x)
}


LVdxdp = function(t,y,p,more=NULL)
{
  x = array(0,c(dim(y),dim(y)[2],length(p)))
  x[,1,1,1] = 1
  x[,1,1,2] = -y[,2]
  x[,1,2,2] = -y[,1]
  x[,2,1,4] = y[,2]
  x[,2,2,3] = -1
  x[,2,2,4] = y[,2]
  return(x)
}

return(list(fn = LV, dfdx = LVdx, dfdp = LVdp, d2fdx2 = LVdx2, d2fdxdp = LVdxdp))
}