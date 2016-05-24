Grow_Method_Sim <-
function(Knoten_max, N_Varis, Wsk, uS=uS, oS=oS){
  Wsk_Abnahme<- (Wsk-0.001)/Knoten_max

  # Zusammenhangstabelle und Knotenanzahl
    Baum<- Baum_zufaellig(Knoten_max=Knoten_max, Wsk=Wsk, Wsk_Abnahme=Wsk_Abnahme)

  # Variablen- und Cutoff-Auswahl
    Baum$Varis<- Var_ZHT_Sim(Knoten=Baum$Knoten, N_Varis=N_Varis)
    Baum$CO   <- Cutoff_ZHT_Sim(Knoten=Baum$Knoten, uS=uS, oS=oS)

  # Ausgabe:
    return(Baum)
}
