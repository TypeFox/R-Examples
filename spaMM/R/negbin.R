negbin <- function(theta=1,link="log") {
  fam <- MASS::negative.binomial(theta=theta,link=link)
  fam$family <- structure("negbin",
                          withArgs=quote(paste("Negative Binomial(",.Theta,")",sep="")))
  ## keep theta as environment(negbin(1)$aic)$.Theta
  return(fam) 
}
#nb <- negbin(1)
#nb$family
#eval(attr(nb$family,"withArgs"),envir=environment(nb$aic))
#assign(".Theta",2,envir = environment(nb$aic))
#eval(attr(nb$family,"withArgs"),envir=environment(nb$aic))
