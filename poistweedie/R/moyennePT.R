moyennePT <-function(p,omega,theta0)
{
## fonction qui permet  de calculer mu=k'(w)
mu <- double(1)
alpha <- double(1)

## calcul de alpha
alpha <- (p-2)/(p-1)

## calcul de la moyenne mu

if (p>1 && p!=2) 
  {
    mu <- exp(omega)*(((1-p)*(exp(omega)-1+theta0))^(alpha-1))
  }
  
if (p==1) 
  {
    mu <- exp(omega)*exp((exp(omega)-1+theta0))
  }
if (p==2)
  {
    mu <- -exp(omega)/(exp(omega)-1+theta0)
  }     
mu
}

