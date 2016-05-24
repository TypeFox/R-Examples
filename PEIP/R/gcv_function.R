gcv_function <-
function(alpha,gamma2,beta)
  {

    f = (alpha^2)/(gamma2 + alpha^2);
   ### length(f);
   ### length(beta);
    if(length(f)>length(beta))
      {
        f=f[1:length(beta)];
      }
    else
      {
        if(length(beta)>length(f))
          {
            iend = length(beta)
            beta=beta[(iend-length(f)+1):iend];
          }
      }
    
    g = (Mnorm(f*beta)^2)/(sum(f))^2;

    return(g)


  }
