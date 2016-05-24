simpson <- function (f, a, b, n, ...)
{
  if (n<1 || a>b) 
  { 
  	stop("Invalid argument")
  }

  h = (b-a)/(2*n)

  s <- f(a,...)
   for (i in 1:n)
  { 
  	s <- s + 4*f(a+(2*i-1)*h,...)
  }
  
  if (n>1)
  { 
    for (i in 1:(n-1))
    { 
    	s <- s + 2*f(a+2*i*h,...)
    }
  }
  s <- s + f(b,...)

  s * h/3
}
