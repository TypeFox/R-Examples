Ctree_Sim <-
function(Data, Varis, CV.Laeufe){

  # Einstellungen
    apriori_v<- table(Data[,1])/nrow(Data)  # Geschaetzte prior

  # Bilden einer Input-Variable fuer die rpart-Funktion
    if (length(Varis)==1){
        if (Varis=="Alle") Varis<- names(Data)[-1]
    }
    Varis<- paste(Varis, sep="") 
    VVn  <- length(Varis)
    EV   <- Varis[1]
    if (VVn>1){
        for(j in 2:VVn) EV<- paste(EV,"+",Varis[j], sep="")
    }
    Modell<- paste(names(Data)[1],"~", EV, sep="")

  # Kreuzvalidierung - k Kreuzvalidierungsgruppen
    Total_GFR  <- 0
    Fehlerraten<- 0

    if (CV.Laeufe!="LOOCV"){
                                   fold<- sample(1:CV.Laeufe, nrow(Data), replace=TRUE)
        if (CV.Laeufe==nrow(Data)) fold<- 1:nrow(Data)                                  
    }
    if (CV.Laeufe=="LOOCV"){ fold     <- 1:nrow(Data)
                             CV.Laeufe<- nrow(Data)
    }
       
    for(k in 1:CV.Laeufe){
        Learning_set <- Data[which(fold!=k),]
        Test_set     <- Data[which(fold==k),]
        rpart.objects<- NULL
        Fehlerrate   <- c()
                
      # Pruning bzgl. aller beta-Werte als Komplexitaetsparameter:
      # Baum lernen
        rpart.objects<- rpart(Modell, 
                              data=Learning_set, method="class", 
                              parms=list(split="gini", prior=apriori_v))

      # Baum testen
        rpart.test <- predict(object=rpart.objects, newdata=Test_set, type="class")
        Tabelle    <- table(rpart.test, Test_set[,1])
        Fehlerrate <- 1-(sum(diag(Tabelle))/sum(Tabelle))
        Fehlerraten<- c(Fehlerraten, Fehlerrate)
    } # Ende Schleife k

    Total_GFR2<- mean(Fehlerraten[-1])  # Loeschen der fuehrenden Null

  # Original-Baum
    Var_Original<- as.vector(subset(rpart.objects$frame$var,rpart.objects$frame$var != "<leaf>"))

  # Geprunter Baum
    rpart_pruned<- prune(rpart.objects, cp=0.1)
    Var_pruned  <- as.vector(subset(rpart_pruned$frame$var,rpart_pruned$frame$var != "<leaf>"))

  # Ausgabe:
    Ausgabe             <- list()
    Ausgabe$Modell      <- Modell
    Ausgabe$Originalbaum<- Var_Original    # Variablen im ungepruntem Baum
    Ausgabe$PrunedBaum  <- Var_pruned      # Variablen im gepruntem Baum
    Ausgabe$FKR         <- Total_GFR2*100  # FKR fuer alle betas 
    Ausgabe$Input_plot_O<- rpart.objects   # Ausgabe des rpart-Objekts zum Plotten
    Ausgabe$Input_plot_P<- rpart_pruned    # Ausgabe des rpart-Objekts zum Plotten incl. Pruning
    return(Ausgabe)
}
