"valchisq" <-
function(matr){

# Giraudoux 08.2003
# renvoie les valeurs du Xi2 d'un tableau de contingence
    
    expect<-chisq.test(matr)$expected
    tabchisq<-expect
    
    for (i in 1:length(matr[,1])){ # de la premiere  la derniere ligne
        for (j in 1:length(matr[1,])){ # de la premiere a la derniere colonne
            tabchisq[i,j]<-(matr[i,j]-expect[i,j])^2/expect[i,j]     
        }
    }
    
    tabchisq
    
}

