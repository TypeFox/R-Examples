dist.l2d.gp <-
function(mean1, var1, mean2, var2, check=FALSE)  
{
  # mean1, var1 :   les (vecteur) moyenne et (matrice de) variance du premier jeu de données.
  # mean2, var2 :   les (vecteur) moyenne et (matrice de) variance du deuxième jeu de données.
   if(check)
    {if(abs(det(var1))<.Machine$double.eps | abs(det(var2))<.Machine$double.eps )
      {stop("One of the sample variances is degenerate")
      }
    }  
  return(sqrt(l2d.gp(mean1, var1, mean1, var1) + l2d.gp(mean2, var2, mean2, var2) - 2*l2d.gp(mean1, var1, mean2, var2)))
}

