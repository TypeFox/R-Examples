Rekombination_VIM3 <-
function(Daten, Proteinkosten, 
                              CV_Laeufe, Max_Knoten, 
                              X_Wsk, X_Art, Brutgroesse,
                              Eltern_Baum1, Eltern_Baum2, VIM,
                              N_Varis, uS, oS){ 
  Nachkommen<- list()
 
  # Falls beide Eltern Wurzelbaeume --> Keine Rekombination
    if (Eltern_Baum1$Knoten==1 && Eltern_Baum2$Knoten==1){
        Nachkommen[[1]]<- Eltern_Baum1
        Nachkommen[[2]]<- Eltern_Baum2
    }else{ # Falls mind. ein Elternteil kein Wurzelbaum ist:
    
  # Soll Rekombination stattfinden?
    temp_crossover<- runif(n=1, min=-0.01, max=1.01)

  if (temp_crossover <= X_Wsk){

  # Standard-Rekombination:
    if (X_Art=="standard"){
    # 1 Wurzelbaum
      if (Eltern_Baum1$Knoten==1 || Eltern_Baum2$Knoten==1){
          Nachkommen<- X_SGP_WK1_VIM(E1=Eltern_Baum1, E2=Eltern_Baum2, VIM=VIM, 
                                     N_Varis=N_Varis)
      }
        
    # 0 Wurzelbaeume
      if (Eltern_Baum1$Knoten!=1 && Eltern_Baum2$Knoten!=1){
          Nachkommen<- X_SGP_VIM2(E1=Eltern_Baum1, E2=Eltern_Baum2, 
                                  N_Varis=N_Varis, uS=uS, oS=oS, 
                                  Max_Knoten=Max_Knoten, Daten=Daten, CV_Laeufe=CV_Laeufe,
                                  VIM=VIM)

      }
    }

  # Brood-Rekombination:
    if (X_Art=="brood"){
        Nachkommen<- X_Brood_VIM2(Daten=Daten, E1=Eltern_Baum1, E2=Eltern_Baum2, 
                                  Max_Knoten=Max_Knoten, Brutgroesse=Brutgroesse, 
                                  Kostenmatrix=Proteinkosten, 
                                  CV_Laeufe=CV_Laeufe, VIM=VIM, N_Varis=N_Varis, uS=uS, oS=oS)
    }

  # 1-Punkt-Rekombination nach Poli:
    if (X_Art=="poli"){ 
        Nachkommen<- X_1Point_VIM2(E1=Eltern_Baum1, E2=Eltern_Baum2, VIM=VIM,
                                   Max_Knoten=Max_Knoten, Daten=Daten, CV_Laeufe=CV_Laeufe,
                                   N_Varis=N_Varis, uS=uS, oS=oS)
    }
  } # temp_crossover <= X_Wsk
  
  if (temp_crossover > X_Wsk){
      Nachkommen[[1]]<- Eltern_Baum1
      Nachkommen[[2]]<- Eltern_Baum2 
  }
  } # Beide Eltern sind nicht Wurzelbaeume!    

  # Ausgabe:
    Ausgabe<- list()
    temp3  <- runif(n=1, min=0, max=1)
      if (temp3 <= 0.5) Ausgabe<- Nachkommen[[1]]
      if (temp3 >  0.5) Ausgabe<- Nachkommen[[2]]
    return(Ausgabe)
}
