Mutation_Cutoff100_4 <-
function(Tree, Elternbaeume, Mutationswsk, 
                                temp_mc_count, gen=100, alpha=0.85, sigma=10, Erfolg,
                                Daten, uS, oS, CV.Lauf){

  # Berechnung der maximalen FKR in der Elternpopulation (wichtig fuer Erfolgsmessung)
    FKR_Eltern<- c()
    for (l in 1:length(Elternbaeume)) FKR_Eltern[l]<- Elternbaeume[[l]]$FKR
    FKR_Eltern_max<- max(FKR_Eltern)

  # Soll mutiert werden?
    temp<- runif(n=1, min=-0.01, max=1.01)

    if (temp<=Mutationswsk){

      # Zaehlen der Cutoff-Mutationen (wichtig fuer Schrittweitenanpassung)
        temp_mc_count<- temp_mc_count+1 
  
      # Erstellung des mutierten Nachkommen
        Tree   <- Gauss_M100_4(Tree=Tree, M=temp_mc_count, gen=gen, alpha=alpha, 
                               sigma=sigma, FKR_max=FKR_Eltern_max, Erfolg=Erfolg,
                               Daten=Daten, uS=uS, oS=oS, CV.Lauf=CV.Lauf)
        sigma  <- Tree$sigma
        Erfolg <- Tree$Erfolg
        Tree   <- list(Zusammenhangstabelle=Tree$Zusammenhangstabelle, 
                       Knoten=Tree$Knoten, 
                       Varis=Tree$Varis, CO=Tree$CO)
    } # Ende Mutation

  # Ausgabe:
    Ausgabe              <- list()
    Ausgabe[[1]]         <- Tree
    Ausgabe$temp_mc_count<- temp_mc_count
    Ausgabe$Gauss_sigma  <- sigma
    Ausgabe$Gauss_Erfolg <- Erfolg
    return(Ausgabe)
}
