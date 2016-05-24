# ' Clean Z's columns based on p-values (coefficients or global)
cleanZtest<-function (Z = Z, X = X, pvalmin = 0.05, methode = 1, adj = T,global=F,bonferroni=F) 
{
  if(global){
    if(bonferroni){
      pvalmin=pvalmin/length(quicol)
    }
    quicol = which(colSums(Z) != 0)
    for (i in quicol) {
      qui = which(Z[, i] != 0)
      Xloc = X[, qui]
      Yloc = X[, i]
      lmloc=lm(Yloc~.,data=data.frame(Xloc))
      summar=summary(lmloc)
      global_pval=pf(summar$fstatistic[1],summar$fstatistic[2],summar$fstatistic[3],lower.tail=FALSE)   
      if (global_pval>pvalmin) {#equation pourrie, on la supprime
        Z[, i] = 0
      }
    }
     }else{#on regarde chaque coef donc on boucle jusqu'a stabilite
    pvalminini=pvalmin
    change=(colSums(Z) != 0)#colonnes a regarder
    quicol=which(change)
    while(length(quicol)>0 ){
      if(bonferroni){pvalmin=pvalminini/sum(Z[,quicol])}
      for (i in quicol) {#nettoyage
        qui = which(Z[, i] != 0)
        Xloc = X[, qui]
        Yloc = X[, i]
        lmloc=lm(Yloc~.,data=data.frame(Xloc))
        summar=summary(lmloc)
        coefs_pval=coef(summar)[,4]#p-values des coefficients
        #on elague juste les coefs pourris
        if(length(Z[qui,i][coefs_pval[-1]>pvalmin])>0){#on a des choses a changer
          Z[qui,i][coefs_pval[-1]>pvalmin]=0#on supprime
        }else{#rien n'a change donc on marque la colonne
          change[i]=F
        }         
      }
      quicol=which(change & (colSums(Z) != 0))#on doit reverifier les colonnes modifiees encore non nulles
    }
  }
  
  return(Z)
}