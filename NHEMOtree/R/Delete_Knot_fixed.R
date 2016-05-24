Delete_Knot_fixed <-
function(Tree, Knoten, N_Varis, uS, oS){

# Falls Wurzelbaum, keine Reduktion sondern neuer zufaelliger Wurzelbaum!
  if (Tree$Knoten==1){
      WK<- Baum_zufaellig_komplett_Sim(Knoten_max=1, Wsk=0.8, N_Varis=N_Varis, uS=uS, oS=oS)
      return(WK)
  }

# Falls kein Wurzelbaum, Loeschen des Knotens!
  ZHT   <- Tree$Zusammenhangstabelle
  Zeilen<- which(ZHT[,1]==Knoten)

  ## Falls beide Tochterknoten des auszuschneidenen Knotens nicht leer sind, ABBRUCH!
  #  if (ZHT[Zeilen[1],3] != -99 && ZHT[Zeilen[2],3] != -99){
  #      stop(paste("'Delete_Knot_fixed': Deletion of node impossible!"))
  #  }

  # Falls beide Tochterknoten des auszuschneidenen Knotens Blaetter sind, loeschen des Knotens!
    if (ZHT[Zeilen[1],3] == -99 && ZHT[Zeilen[2],3] == -99){
        Tree_kurz<- Delete_Branch_fixed(Tree=Tree, Knoten=Knoten) 
        return(Tree_kurz)
    }      

  Tochterknoten<- max(ZHT[Zeilen[1],3], ZHT[Zeilen[2],3])
  ZHT_kurz     <- ZHT[-Zeilen,]
  ZHT_kurz[which(ZHT_kurz[,3]==Knoten), 3]<- Tochterknoten

  Varis_kurz<- Tree$Varis[-which(Tree$Varis[,1]==Knoten),]
  CO_kurz   <- Tree$CO[-which(Tree$CO[,1]==Knoten),]

  if (is.vector(Varis_kurz)==TRUE) Varis_kurz<- t(as.matrix(Varis_kurz))
  if (is.vector(CO_kurz)==TRUE)    CO_kurz   <- t(as.matrix(CO_kurz))
  
  # Ausgabe:
    Ausgabe<- list()
    Ausgabe$Zusammenhangstabelle<- ZHT_kurz          # Verkuerzte ZHT
    Ausgabe$Knoten              <- nrow(ZHT_kurz)/2  # Verringerte Knotenanzahl
    Ausgabe$Varis               <- Varis_kurz        # Verkuerzte Variablentabelle
    Ausgabe$CO                  <- CO_kurz           # Verkuerzte Cutofftabelle
    return(Ausgabe)
}
