cor.arma = function(mod)
{
# calcul des corrélations des estimateurs dans une modélisation arima
# input un modele, résultat d'un appel à arima()
# output la matrice des correlations des estimateurs
aa= mod$var.coef
bb = diag(diag(aa)^(-.5))
cc = bb%*% aa%*%bb
dimnames(cc) = dimnames(aa)
cc
}
