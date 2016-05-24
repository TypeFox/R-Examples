lipper<-function(a,b,c,d=0,f=0,g=-1)
  {
    #######  given the cartesian parameters for an ellipse
    ###  return the parametric form, a, b, phi
    ###  
    if(missing(d)) d = 0
    if(missing(f)) f = 0
    if(missing(g)) g = -1

    acot<-function(z)
      {
        return( atan(1/z) )
      }
    
    
    M = matrix(c(a,b,d, b,c,f, d, f, g), ncol=3, byrow=TRUE )
    
    DEL = det(M)
    J = det(M[1:2, 1:2])
    I = a +c

    if(DEL==0 | J<= 0 | DEL/I >=0)
      {
        print("NOT an Ellipse")
        return(NULL)
      }
    
    ##  if a = c the ellipse is a circle (do something else?)

    BBAC = b*b - a*c
    
    x0 = (c*d - b*f)/BBAC 
    y0 = (a*f - b*d)/BBAC 


    num1 = 2*(a*f^2 +c*d^2 + g*b^2 -2*b*d*f -a*c*g)
    num2 = sqrt((a-c)^2 +4 * b^2)

    
    ap = sqrt(num1 / (BBAC*(num2 - I)))
    bp = sqrt(num1 / (BBAC*(-1*num2 - I)))

    if(b==0 & a<c)
      {
        phi = 0

      }
    if(b==0 & a>c)
      {
        phi = pi/2
      }
    if(b!=0 & a<c)
      {
        phi = 0.5*acot( (a-c)/(2*b))
      }
    if(b!=0 & a>c)
      {
        phi = pi/2 + 0.5*acot( (a-c)/(2*b))
      }
    
    
    return( c(ap, bp, phi) )
  }
