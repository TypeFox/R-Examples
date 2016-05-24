gausskernel <-
function(X=NULL,sigma=NULL)
      {
       return(exp(-1*as.matrix(dist(X)^2)/sigma))
       }

