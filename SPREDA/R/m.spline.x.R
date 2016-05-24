m.spline.x <-
function(x,tt,i,k)
{
  ti=tt[i]
  tik=tt[i+k]
  if(x<ti|x>=tik)#x<ti | x>tik |ti==tik)
  {
    res=0
    return(res)
  }else{
    if(k==1)
    {
      res=1/(tik-ti)
      return(res)
    }else{
      a1=m.spline.x(x,tt,i,k-1)
      a2=m.spline.x(x,tt,i+1,k-1)
      res=(k*((x-ti)*a1+(tik-x)*a2))/((k-1)*(tik-ti))
      return(res)
    }
  }
}
