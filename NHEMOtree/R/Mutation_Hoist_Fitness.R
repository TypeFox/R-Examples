Mutation_Hoist_Fitness <-
function(Tree, Datensatz, CV_Laeufe){

  # Falls Tree Wurzelbaum ist, keine Mutation:
    if (Tree$Knoten==1) return(Tree)

  # Zusammenhangstabelle (ohne Wurzelknoten)
    Tree   <- Werte_Reduktion_ZHT2(Tree)
    FKR    <- FKR_CV(Tree=Tree, Daten=Datensatz, CV_Laeufe=CV_Laeufe)
    ZHT_oWK<- Tree$Zusammenhangstabelle[3:nrow(Tree$Zusammenhangstabelle),]

  # Erstellung aller moeglichen Subbaeume aus Tree
    VK     <- ZHT_oWK[2*(1:(nrow(ZHT_oWK)/2)), 1] # moegliche Subbaum-Wurzelknoten
    Subbaum<- list()
    Sub_FKR<- c()
    for (i in 1:length(VK)){
         Subbaum[[i]]<- Mutation_Hoist(Tree=Tree, KC=VK[i])
         Sub_FKR[i]  <- FKR_CV(Tree=Subbaum[[i]], Daten=Datensatz, CV_Laeufe=CV_Laeufe)

    }

  # Falls FKR eines Subbaums < FKR Ausgangsbaum, Auswahl des Subbaums
    if (min(Sub_FKR)<=FKR){
        temp<- which(Sub_FKR==min(Sub_FKR))
        if (length(temp)>1) temp<- sample(temp, 1)
        Tree<- Subbaum[[temp]]
    }

  # Ausgabe:
    return(Tree)
}
