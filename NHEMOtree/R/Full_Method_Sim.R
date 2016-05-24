Full_Method_Sim <-
function(Knoten_max, N_Varis, uS=uS, oS=oS){

  # Erstellen des Wurzelknotens
    WK<- ZHT_Wurzelknoten(1)

  # Vervollstaendigung der Zusammenhangstabelle
    if (Knoten_max==1){
        WK$Zusammenhangstabelle[,3]<- -99
        WK$Knoten                  <- 1
        Baum                       <- WK
    }

    if (Knoten_max>1){
        Baum<- Vater_in_ZHT(ZHT=WK[[1]], zaehler=WK[[2]], 
                            Knoten_max=Knoten_max, Wsk=1, Wsk_Abnahme=0)
    }

  # Variablen- und Cutoff-Auswahl
    Baum$Varis<- Var_ZHT_Sim(Knoten=Baum$Knoten, N_Varis=N_Varis)
    Baum$CO   <- Cutoff_ZHT_Sim(Knoten=Baum$Knoten, uS=uS, oS=oS)

  # Ausgabe:
    return(Baum)
}
