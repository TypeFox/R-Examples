Baum_zufaellig <-
function(Knoten_max, Wsk, Wsk_Abnahme){

    WK<- ZHT_Wurzelknoten(Wsk)

    Erg<- Vater_in_ZHT(ZHT=WK[[1]], zaehler=WK[[2]], Wsk=Wsk, 
                           Wsk_Abnahme=Wsk_Abnahme, Knoten_max=Knoten_max)

  # If root tree
    if (Knoten_max==1){
        Erg$Zusammenhangstabelle[,3]<- -99
        Erg$Knoten<- 1
    }

    return(Erg)
}
