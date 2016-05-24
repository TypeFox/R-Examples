.init.exp.model.variogram   = function(dist,gamma){
  ###############################################################################################
  # cette fonction estime par moindres carres les parametres de  lissage d'un semi variogram
  # Le modele choisi est fixe (modele Exp)
  # Ces parametres sont utilises comme valeurs initiales du programme d'ajustement fit.variogram appele  
  #       par la fonction krige (quand on est en mode auto)
  ###############################################################################################
  gammaFunc = function(h, A, a){
    A*(1 -exp(-3*h/a)) 
    A*(1 -exp(-h/a)) 
  }
    objectif = function(par,dist,gamma){
    A = par[1]
    a = par[2]
    obj = sum((gammaFunc(dist,A,a) - gamma)^2)
    obj
  }
  par = c(mean(gamma), max(dist)/2)
  opt = optim(par,objectif,dist =  dist, gamma =gamma,lower = c(-Inf,0.02),method = "L-BFGS-B")
  if(opt$convergence!=0) return(NULL)
  c(psill = opt$par[1],range = opt$par[2])
}

  

