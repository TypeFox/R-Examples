`as.PhViD` <-
function(DATA.FRAME,MARGIN.THRES=1){
## Fonction qui va permettre de calculer le nombre de notifications pour chaque "médicament" et "effet indésirable".

## FILE :         Emplacement + nom du fichier
##		            Les données doivent contenir 3 colonnes 
##			           - le label ATC
##			           - le label meddra
##			           - le nombre de notifications		

## MARGIN.THRES :  Seuil sur les marges. on peut considérer uniquement les lignes et colonnes qui dépassent
##                un certain effectif

## FILE="ATC7-pt.csv"
## FILE="ATC5-hlt.csv"
## rm(list=ls())

#############################################################################
##			 
##			  Event  
##   			 ______________
## 		Drug    |  n11  |  n10 | <- Drug.marge (n1.)
##			|_______|______|
## 			|  n01  |  n00 |
##			|_______|______|
##          		   |         	N (n..)
##			Event.marge (n.1)
##############################################################################                  

data <- DATA.FRAME
data[,1] <- as.factor(DATA.FRAME[,1])
data[,2] <- as.factor(DATA.FRAME[,2])
data[,3] <- as.double(DATA.FRAME[,3])
#data <- read.table(file=FILE,header=TRUE,sep=";",colClasses=c("factor","factor","double")) 

coln <- names(data) # On stocke le nom des colonnes
names(data)[3] <-"n11"
data_cont <- xtabs(n11 ~ .,data=data) # "matrice" de contingence
n1._mat <- apply(data_cont,1,sum)     # marges des lignes
n.1_mat <- apply(data_cont,2,sum)     # marges des colonnes


## Nettoyage en fonction du seuil sur les marges
if (MARGIN.THRES > 1){
  while(sum(n1._mat<MARGIN.THRES)>0 | sum(n.1_mat<MARGIN.THRES)>0){             
    data_cont <- data_cont[n1._mat  >= MARGIN.THRES,] # on supprime les lignes dont la marge n'est pas supérieure au seuil
    data_cont <- data_cont[,n.1_mat >= MARGIN.THRES] # on supprime les colonnes dont la marge n'est pas supérieure au seuil
    n1._mat <- apply(data_cont,1,sum) # on recalcule les marges des lignes...
    n.1_mat <- apply(data_cont,2,sum) # ...et des colonnes
  }
}
##---------------------------------------
coord <- which(data_cont != 0, arr.ind=TRUE) # on récupère les coordonnées des notifications (traduit par un nombre >=1 dans la matrice data_cont)
coord <- coord[order(coord[,1]),] # la fonction which les a ordonnées par colonnes. On les veut ordonnées par lignes

#data_net <- cbind(rownames(data_cont)[coord[,1]],colnames(data_cont)[coord[,2]],data_cont[coord]) # on compile les
#libellés des médicaments, des effets indésirables et des notifications après le nettoyage sur les marges. (si pas de
#nettoyage, cela revient à prendre data)

# on récupère les nouvelles dimensions de la matrice
Nb_n1. <- length(n1._mat) # nombre de lignes
Nb_n.1 <- length(n.1_mat) # nombre de colonnes

libel.medoc <- rownames(data_cont)[coord[,1]] # on conserve les libellés des médicaments qui restent
libel.effet <- colnames(data_cont)[coord[,2]] # on conserve les libellés des effets indésirables qui restent
n11 <- data_cont[coord] 
N <- sum(n11) # le nombre total de notifications
n1. <- n1._mat[coord[,1]] # on affecte à chaque notification sa marge "ligne"...
n.1 <- n.1_mat[coord[,2]] # ... et sa marge "colonne"


RES <- vector(mode="list")
#RES$PARAM <- data.frame(FILE,MARGIN.THRES) # input parameters # FILE, not necessary
RES$L <- data.frame(libel.medoc,libel.effet)
colnames(RES$L) <- coln[1:2]
RES$data <- cbind(n11,n1.,n.1) # matrice "équivalente" à la matrice de départ avec en plus les marges et effectifs attendus
rownames(RES$data) <- paste(libel.medoc,libel.effet)
RES$N <- N # nb de notifications total
RES
}

