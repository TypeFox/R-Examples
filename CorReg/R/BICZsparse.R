# ' calcul du BIC en format creux
# ' @param Z_Zi numero des lignes des 1 dans la matrice Z
# ' @param vecteur du nombre de 1 dans chaque colonne
# ' @param Bic_vide_vect BIC vectoriel de la matrice nulle
# ' @param BicOld BIC vectoriel de la matrice precedente
# ' @param methode 1:utilisation de la fonction householderQr, 2:utilisation de la fonction colPivHouseholderQr
# ' @param Zold_Sj nombre de 1 dans la colonne de la matrice precedente
# ' @param nrow nombre de ligne de la matrice x
BICZsparse<-function(X=X, Z_Zi=Z_Zi, Z_Sj=Z_Sj, Bic_vide_vect=Bic_vide_vect, BicOld=BicOld, methode=methode, Zold_Sj=Zold_Sj, nrow=nrow){
  res=.Call( "BICZsparse", X, Z_Zi, Z_Sj, Bic_vide_vect, BicOld, methode, Zold_Sj, nrow, PACKAGE = "CorReg")
  return (res$BIC)
}