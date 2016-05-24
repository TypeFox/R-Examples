rnk <-
function(G, tol=1e-14)
  {
    ###  tries to duplicate the matlab function rank
    ###   which returns the number of non-zero singular values
    SV=svd(G);
    return(length(which(SV$d>tol)))
  }
