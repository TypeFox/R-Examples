SMfzero <-
function(f,x1,x2,num=1000,eps=1e-5,eps1=1e-5)
  {
    i=0
    while((abs(x1-x2)>eps) & (i<num))
    {
      c=x2-f(x2)*(x2-x1)/(f(x2)-f(x1))
      x1=x2
      x2=c
      i=i+1
    }
    print(x2)
    print(f(x2))
    if(abs(f(x2))<eps1)
    {
      print("finding root is successful")
    }
    else
      print("finding root is fail")
  }
