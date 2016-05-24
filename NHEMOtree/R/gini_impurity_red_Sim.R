gini_impurity_red_Sim <-
function(Daten, index_set, vari_index, Varis, cutoff){
  
  # Betrachtete Variable
    variable<- Varis[which(Varis[,1]==vari_index), 2]

  # Linke und Rechte Teilmenge der Beobachtungen berechnen
    left_index_set <- intersect(index_set, which(Daten[, variable]<=cutoff))
    right_index_set<- intersect(index_set, which(Daten[, variable]>cutoff))
    
  # Tochterknoten
    class      <- Daten[index_set,1]
    class_left <- Daten[left_index_set,1]
    class_right<- Daten[right_index_set,1]

  # impurity reduction
    N_class      <- length(class)
    n_classL     <- length(class_left)/N_class
    n_classR     <- length(class_right)/N_class
    imp_reduction<- gini_index(class) - n_classL*gini_index(class_left) - n_classR*gini_index(class_right)
    
  # Ausgabe:
    return(imp_reduction)
}
