Init_RHaH <-
function(N, Knoten_max, Wsk, N_Varis, uS=uS, oS=oS){
  # Aufteilen der Anzahl zu erstellender Baeume auf die beiden Methoden:
  N1<- round(N/2)
  
  # Baumerstellung
  Baum<- list()
  for (i in 1:N1){
    Knoten   <- sample(1:Knoten_max, 1)  # Verschiedene Baumgroessen erwuenscht!
    if (Knoten_max==1) Knoten<- 1
    Baum[[i]]<- Grow_Method_Sim(Knoten_max=Knoten, Wsk=Wsk, N_Varis=N_Varis, uS=uS, oS=oS)
  }
  
  for (i in (N1+1):N){
    Knoten   <- sample(1:Knoten_max, 1)  # Verschiedene Baumgroessen erwuenscht!        
    if (Knoten_max==1) Knoten<- 1
    Baum[[i]]<- Full_Method_Sim(Knoten_max=Knoten, N_Varis=N_Varis, uS=uS, oS=oS)
  }
  
  # Ausgabe:
  return(Baum)
}
