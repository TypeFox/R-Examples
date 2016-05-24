.isInputOk.prevR = function(...){
  ###############################################################################################
  # cette fonction sert a tester la validite des parametres d'entree N R et var
  # ... est une liste contenant N R (et var)
  # De nombreuses fonctions ont pour parametres d'entree optionnels N R et parfois var
  # N et R correspondent aux parametres des differents rings calcules et var contient la variable sur laquelle on veut travailler
  # Par exemple on veut faire un krigeage de la variable r.prev pour les prametres N = 100 R = Inf
  # Par exemple on veut faire un lissage par noyau pour les prametres N = 100 R = Inf  
  # L'utilisateur peut fournir plusieurs N plusieurs R et plusieurs var cependant nous avons defini des contraintes d'utilisation
  # L'idee premiere est que l'utilisateur doit fournir tous les couples (ou triplet) N R  (et var )
  #  Exemple N = c(100, 100, 200, 400) R = c(Inf, Inf, Inf, Inf) var = c("r.prev", "r.wprev", "r.prev", "r.wprev")
  #         Danns ce cas on traite 4 cas
  #            N = 100 R = Inf var = "r.prev"
  #            N = 100 R = Inf var = "r.wprev"
  #            N = 200 R = Inf var = "r.prev"
  #            N = 400 R = Inf var = "r.wprev"
  # Pour donner de la souplesse a l'utilisateur on autorise a chacune des variables (N R  et var) d'etre
  #   soit de longueur 1
  #   soit d'etre de longueur egal a la longueur de la variable de plus grande longueur ainsi
  #   N = c(100, 100, 200, 400) R = Inf var = c("r.prev", "r.wprev", "r.prev", "r.wprev")  est correct 
  #         (On considere que R contient 4 fois Inf)  
  #   N = c(100, 200, 300, 400) R = Inf var = "r.prev"  est correct 
  #         (On considere que R contient 4 fois Inf et var 4 fois r.prev)
  #   N = c(100, 200, 300, 400) R = Inf var = c("r.prev","r.wprev")  est incorrect 
  # cette fonction teste aussi si les variables N et R sont numeriques
  # Cette fonction arrete la simulation si une erreur est detecte
  ###############################################################################################
  N= list(...)$N
  if(!is.null(N) && !is.numeric(N)){
    stop("N is not numeric.", call.=F)
  }
  R= list(...)$R
  if(!is.null(R) && !is.numeric(R)){
    stop("R is not numeric.", call.=F)
  }
  nc = sapply(list(...),length)
  if((sum(nc==1) + sum(nc == max(nc))) == length(nc)) return(NULL)
  if(all(nc==nc[[1]])) return(NULL)
  if(sum(nc!=1)!=1) {
    nom.param = paste(names(nc),sep=", ")
    stop(gettextf("the input parameters are not valid: %s must have the same length or a length equal to 1.",nom.param,domain="R-prevR"), call.=F)
  }
  NULL
}

