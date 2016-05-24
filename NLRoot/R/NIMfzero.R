NIMfzero <-
function(f,f1,x0=0,num=100,eps=1e-5,eps1=1e-5)
  {
    a=x0;b=a-f(a)/f1(a)
    i=0
    while((abs(b-a)>eps)&(i<num))
    {
      a=b
      b=a-f(a)/f1(a)
      i=i+1
    }
    print(b)
    print(f(b))
    if(abs(f(b))<eps1)
    {
      print("finding root is successful")
    }
    else
      print("finding root is fail")
  }
