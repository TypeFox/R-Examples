knots.transform <-  function( d , alpha = 0, symmetric = TRUE)
  {
    knots <-  c(0,1)
    if(d>0) {
      if (symmetric == FALSE)
        {
          for ( k in 1:d)
            {
              knots.new <-  knots
              for (j in 1: length(knots))
                {
                  knots.new <-  c( knots.new, knots[j-1]+ (exp(alpha)/(1+exp(alpha))) * (knots[j] - knots[j-1]))
                }
              knots <-  sort(knots.new)
            }
        }
      else {
        knots <- knots.transform( d-1,alpha = alpha, symmetric=FALSE)
        knots <-  sort(unique(c(knots *0.5 , 1- 0.5 * knots)))
      }
    }
    return(knots)
  }
