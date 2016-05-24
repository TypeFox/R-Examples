Baum_zufaellig_komplett_Sim <-
function(Knoten_max, N_Varis, Wsk=0.8, Wsk_Abnahme=0.1, uS=uS, oS=oS){
  # Tree
    Baum<- Baum_zufaellig(Knoten_max=Knoten_max, Wsk=Wsk, Wsk_Abnahme=Wsk_Abnahme)

  # Variables
    Baum$Varis<- Var_ZHT_Sim(Knoten=Baum$Knoten, N_Varis=N_Varis)
 
  # Cutoffs
    Baum$CO   <- Cutoff_ZHT_Sim(Knoten=Baum$Knoten, uS=uS, oS=oS)

  # Output
    return(Baum)
}
