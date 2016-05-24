modele1 <- function(eps,thetavec,xvec) {

  n <- length(eps)

  p <- length(thetavec) # y compris l'intercept
  
  xvec <- matrix(xvec,nrow=n,ncol=p,byrow=FALSE) # on le remet sous forme de matrice. !! Doit déjà contenir la colonne de 1's

  yvec <- xvec %*% as.matrix(thetavec) + eps # on fabrique des y_i's à partir des epsilon_i's

  thetavec <- solve(t(xvec) %*% xvec) %*% t(xvec) %*% as.matrix(yvec) # thetachap   On estime les paramètres
  eps <- yvec - xvec %*% as.matrix(thetavec) # epschap  On calcule les résidus

# equivalent de:
  # eps <- residuals(lm(yvec~xvec))

 
  return(eps)

}
