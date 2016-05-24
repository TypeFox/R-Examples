sintegral<-function(x,fx,n.pts=256){
  ## numerically integrates fx over x using Simpsons rule
  ## x - a sequence of x values
  ## fx - the value of the function to be integrated at x
  ## n.pts - the number of points to be used in the integration


  n.x<-length(x)

  if(n.x!=length(fx))
    stop("Unequal input vector lengths")

  if(n.pts<64)
    n.pts<-64

  if(length(x)>n.pts)
      n.pts<-length(x)

  ## use linear approximation to get equally spaced x values


  ap<-approx(x,fx,n=2*n.pts+1)

  h<-diff(ap$x)[1]

  integral<-h*(ap$y[2*(1:n.pts)-1]
               +4*ap$y[2*(1:n.pts)]
               +ap$y[2*(1:n.pts)+1])/3

  (list(x=ap$x[2*(1:n.pts)],y=cumsum(integral),int = sum(integral)))
}
