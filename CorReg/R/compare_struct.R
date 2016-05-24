#' To compare sub-regression structures
#' @export
#' @description Compares two sub-regression structures, considering one of them as the "true one".
#' @param trueZ first structure (binary adjacency matrix)
#' @param Zalgo second structure (binary adjacency matrix)
#' @param all (boolean) Also compute the ratio for each stat.
#' @param mode how to modify the structures before comparison. mode=c("NULL","hybrid","clique","sym")
#' It allows to compare groups instead of exact sub-regressions. Does nothing by default.
#' @return
#' \item{true1}{ Number of links that exist in both matrices }
#' \item{false1}{ Number of links that exist only in Zalgo }
#' \item{false0}{ Number of links that exist only in trueZ }
#' \item{deltadr}{Number of sub-regressions in trueZ -  Number of sub-regressions in Zalgo (i.e. : negative if too much sub-regressions in Zalgo)}
#' \item{true_left}{ Number of variables redundant in both matrices }
#' \item{false_left}{ Number of variables redundant in Zalgo but not in trueZ  }
#' \item{ratio_true1}{ ratio of links in trueZ that exist also in Zalgo  }
#' \item{ratio_true0}{ ratio of links not in trueZ that don't exist in  Zalgo.}
# ' \item{ratio_false1}{ ratio of links that exist only in Zalgo  }
# ' \item{ratio_false0}{ ratio of links that exist only in trueZ  }
compare_struct<-function(trueZ=trueZ,Zalgo=Zalgo,all=TRUE,mode="NULL"){
  if(mode=="hybrid"){
    trueZ=trueZ+t(trueZ)+t(trueZ)%*%trueZ#attention ? bien multiplier par la transpos?e ? gauche
    trueZ[trueZ>1]=1     
    Zalgo=Zalgo+t(Zalgo)+t(Zalgo)%*%Zalgo#attention ? bien multiplier par la transpos?e ? gauche
    Zalgo[Zalgo>1]=1      
  }else if(mode=="clique"){#on fait des cliques
    trueZ=cliquefaction(trueZ)
    Zalgo=cliquefaction(Zalgo)
  }else if (mode=="sym"){#on fait juste la d?sorientation
    trueZ=trueZ+t(trueZ)
    Zalgo=Zalgo+t(Zalgo)
  }else{
    #on garde les structures d'origine
  }
  trueZ=as.matrix(trueZ)
  Zalgo=as.matrix(Zalgo)
  res=as.matrix(trueZ-Zalgo)
  nbbon1=sum(Zalgo*trueZ)#produit de hadamard
  nbtrop=-sum(res[res==-1])#faux 1
  nbmank=sum(res[res==1]) #faux 0
  nbbon0=ncol(trueZ)^2-nbbon1-nbtrop-nbmank #vrai 0
  #dist=nbtrop+nbmank
  #on passe aux pourcentages
  taux_bon1=nbbon1/sum(trueZ)#taux de v?rit? d?couverte
  taux_bon0=nbbon0/(ncol(trueZ)^2-sum(trueZ))
  taux_faux1=nbtrop/sum(Zalgo) #taux d'ajout dans ce qui est dit doit tendre vers 0 (max si on ne dit que des aneries)
  taux_faux0=1-taux_bon0 #taux d'oublis doit tendre vers 0 (max si on n'a rien dit de vrai)
  vraissreg=which(colSums(trueZ)>0) 
  ssregalgo=which(colSums(Zalgo)>0)
  deltap2=length(vraissreg)-length(ssregalgo)
  bon_gauche=sum(duplicated(c(vraissreg,ssregalgo)))
  faux_gauche=length(ssregalgo)-bon_gauche
  if(all){
     #return(list(ratio_true1=taux_bon1,ratio_true0=taux_bon0,ratio_false1=taux_faux1,ratio_false0=taux_faux0,true1=nbbon1,false1=nbtrop,false0=nbmank,deltadr=deltap2,true_left=bon_gauche,false_left=faux_gauche))
     return(list(ratio_true1=taux_bon1,ratio_true0=taux_bon0,true1=nbbon1,false1=nbtrop,false0=nbmank,deltadr=deltap2,true_left=bon_gauche,false_left=faux_gauche))
  }else{
    return(list(ratio_true1=taux_bon1,ratio_false1=taux_faux1))
  }
}