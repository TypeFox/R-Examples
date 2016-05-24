chi2inv <-
function(x, n)
  {
###########  Rick Aster chi2nv function is simple in R
    ### use qchisq
    return(qchisq(x, n))

  }
