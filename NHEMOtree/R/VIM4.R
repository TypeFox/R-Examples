VIM4 <-
function(ZHT, Varis, Gewichte, Knoten, Ebene, VIM){
  # Welche Variable gehoert zu dem Knoten
    Var<- Varis[Varis[,1]==Knoten, 2]
 
  # Berechnung des Ergebnises fuer diesen Knoten
    VIM[Var]<- VIM[Var] + (1/nrow(Varis))*Gewichte[Ebene]

  # Was ist der linke Knoten
    refKnotenl<- ZHT[ZHT[,1]==Knoten & ZHT[,2]==0, 3]
  if(refKnotenl != -99){
    # Wenn Knoten kein Blatt, dann rufe Berechnung auf fuer linken Knoten und speichere Ergebnis
    VIM<- VIM4(ZHT, Varis, Gewichte, refKnotenl, Ebene+1, VIM)
  }

  # Was ist der rechte Knoten
    refKnotenr<- ZHT[ZHT[,1]==Knoten & ZHT[,2]==1, 3]
    if(refKnotenr != -99){
       # Wenn Knoten kein Blatt, dann rufe Berechnung auf fuer rechten Knoten und speichere Ergebnis
    VIM<- VIM4(ZHT, Varis, Gewichte, refKnotenr, Ebene+1, VIM)
  }
  # Gebe Ergebnis aus
  return(VIM)
}
