omega <-
function(p,mu,theta0)
{
omega <- double(1)
x0 <- double(1)
g  <- numeric(1)
d  <- numeric(1)
alpha <- double(1)
t <- double(1)
a <- numeric(1)
b <- numeric(1)

image <- double(1)


if (p>1 && p!=2) 
{

  ## calcul de alpha
  alpha <- (p-2)/(p-1)

    
  
  ## initialisation des bornes g et d de l'intervalle

	  g <- -(theta0 - 1)
	  d <- 0


  ## initialisation de la solution x0
  x0 <- (g+d)/2
	
  ## r\'{e}duction de l'intervalle ou se trouve la solution de l'\'{e}quation imfx0(x)=0

  while ((d-g)>(10^(-14))) 
    {
	   image <-imfx0(x0,p,mu,theta0)
	   if (image>0) { d <- x0}
	   else { g <- x0}
	   x0 <- (g+d)/2
	   ecart <- d - g
    }  
  ## changement de variable 
  omega <- log(x0)
  
}
if (p==1) 
{


  ## posons les contantes a et b pour simplifier les calculs
  a <- (theta0-1)+log(mu)
  

 
  ## initialisation des bornes g et d de l'intervalle
   if (a < -1)  
      {
        g <- 0
        d <- -a
      }
   else
      {
        g <- -a-1
        d <- -a
      } 

  ## initialisation de la solution x0
  x0 <- (g+d)/2
	
  ## r\'{e}duction de l'intervalle ou se trouve la solution de l'\'{e}quation imfx0(x)=0	
  while ((d-g)>(10^(-14))) 
    {
	    image <-imfx0(x0,p,mu,theta0)
		if (image>0) { d <- x0}
	    else { g <- x0}
	    x0 <- (g+d)/2
      
    }  
  ## changement de variable 
  omega <- x0

}
if  (p==2)
{
  omega <- log((mu * (1 - theta0))/(1 + mu))

}
  omega
}

