imfx0 <-function(x0,p,mu,theta0)
{
## fonction qui permet  de calculer k'(w)
imfx <- double(1)
a <- double(1)
b <- double(1)


if (p>1 && p!=2)  
   {
      ## posons les contantes a et b pour simplifier les calculs
      a <- (1-p)*((mu)^(p-1))
      b <- (theta0 - 1)* a
      ## calcul de l'image de x0
      imfx <- exp((p-1)*log(x0))- (a*(x0))- b 
   }
if (p==1)
   {
      ## posons la contante a  pour simplifier les calculs
      a <- (theta0-1)-log(mu)
      ## calcul de l'image de x0
      imfx <- exp(x0)+ (x0)+ a
   }    
imfx
}

