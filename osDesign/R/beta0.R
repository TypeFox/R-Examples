beta0 <-
function(betaX,
									X,
									N,
									rhoY,
									expandX="all")
{
	##
  value <- optimize(beta0eval,
  									interval=logit(rhoY) + c(-10,10),
                    betaX=betaX,
                    X=X,
                    N=N,
                    rhoY=rhoY,
                    expandX=expandX,
                    tol=.Machine$double.eps^0.5)$minimum
  ##
  return(value)
}
