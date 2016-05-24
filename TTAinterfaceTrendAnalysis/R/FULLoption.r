FULLoption <- function (param, depth=NULL, sal = NULL, site=NULL, rawdata="NO", select="NO", resume.reg="NO",
test.normality="NO", plotB = "NO", selectBox="ByYears", log.trans="NO", plotZ="NO", datashow="NO",
help.timestep = "NO", auto.timestep = "NO", time.step = NULL, help.aggreg = "NO", auto.aggreg = "NO", aggreg = NULL ,  
mix = "YES", outliers.re = "NO", na.replace="NO", start = NULL, end = NULL, months = c(1:12), norm = "NO", npsu = 30,
autocorr = "NO", spectrum="NO", anomaly="NO", a.barplot="NO", zsmooth="NO", local.trend = "NO", test= "MK")

################################################################################################################################################################
################################################################# FULLoption arguments #########################################################################
################################################################################################################################################################
##
##  param = parameters on which analysis will be performed; site = sampling site to treat
##  depth = depths take into account; sal = salinities take into account (without missing values if sal is different from salinity max and salinity min)
##  rawdata = desciptive statistics on raw data ; select = desciptive statistics on selected site and parameter ; resume.reg = desciptive statistics on regularized time series
##  test.nomrality = perform a Shapiro-Wilk normality test on selected parameter ; plotB = display a boxplot (by year) of rawdata
##  plotZ = display a plot of the regularized time series ; datashow = show a table of the regularized data ; time.step = choice of the time step for data aggregation
##  help.timestep = time step advice ; auto.timestep = automaticaly apply the most relevant time step (mean time between two measurments)
##  aggreg = choice of the aggregation method of the data ; help.aggreg = method advice ; auto.aggreg = automaticaly apply the best method (p of Wilcoxon)
##  mix = aggregate the sampling site ; outliers.re = remove extreme values (boxplot method) ; na.replace = replace missing values
##  start and end = select the first and last year to take into account in the analysis ; months = same thing as for the years
##  norm and npsu = compute normalised values of nutrients at the salinity npsu
##  autocorr et spectrum = respectively display the autocorrelation diagramme and Fourrier spectrum of the regularized time series
##  anomaly = display a color boxplot of the anomalies calculated for each time step of each year (value at time step of the year n - the median of the time step of all year)
##  zsmooth = display a detrended plot of the time series ; local.trend = display the cusum plot of the time series ; test = perform the selected temporel trend test
##
################################################################################################################################################################
################################################################################################################################################################
##  How to prepare your data base before importation :
##  your datasheet must be in .txt file, if possible select only the interesting column before creation of the datasheet (parameters, Category names, dates...)
##  name exactly the dates column "Dates", the station one "Category", the depth ones "Depth" and the salinity one "Salinity"
##  date format must be as followed : YYYY-mm-dd  and decimal separator must be '.' (dote)
################################################################################################################################################################
################################################################################################################################################################

{
#____________________________________________________________________________________________________________Pour eviter des 'NOTES' supplementaire lors du check 
 {
  YEARS <- NULL
  DayYears <-NULL
  Category <- NULL
  Salinity <- NULL
  Depth <- NULL
 }
#____________________________________________________________________________________________________________Stat descriptive sur les donnees de base
 {
  if (rawdata == "YES") {                                                                                                                                         # appel la fonction si rawdata = YES
#_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ capture la sortie de summary(Envir$Data)
      sumdata <- capture.output(summary(Envir$Data))
#_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ insert le resultat de la capture dans la fenetre de texte 'Envir$txt2'
      tkinsert(Envir$txt2, "end", paste("-Descriptive statistics on raw data-", "\n"))                     # mise en page
      tktag.add(Envir$txt2, "titre", "end -2 lines linestart","end -2 lines lineend")
      tktag.configure(Envir$txt2, "titre", font=tkfont.create(family="courier",size=10,weight="bold"))
      sumdata <- paste(sumdata, collapse="\n")
      tkinsert(Envir$txt2,"end",paste(sumdata, "\n\n"))

      ll <- ncol(Envir$Data)                                                                                                                                     # comptage de la quantite de donnee par colonne
      res1 <- length(Envir$Data[1][!is.na(Envir$Data[1])])
      for (i in 2:ll) { res2 <- length(Envir$Data[i][!is.na(Envir$Data[i])])
      res1 <- data.frame(cbind(res1,res2)) }
      nn <- names(Envir$Data)
      names(res1) <- nn
      DataQ <- capture.output(res1)
      DataQ <- paste(DataQ, collapse="\n")

      A <- tapply(as.numeric(format(as.Date(Envir$Data$Dates, format="%Y-%m-%d"), format = "%Y")), list(Envir$Data$Category), max, na.rm=TRUE)  # nombre d'annee par Category
      B <- tapply(as.numeric(format(as.Date(Envir$Data$Dates, format="%Y-%m-%d"), format = "%Y")), list(Envir$Data$Category), min, na.rm=TRUE)
      C <- A-B
      Res.Per <- data.frame(cbind(B,A,C))
      names(Res.Per) <- c("Min.","Max.","Length")
      DataP <- capture.output(Res.Per)
      DataP <- paste(DataP, collapse="\n")

      tkinsert(Envir$txt2,"end", paste("Data quantity: " ,"\n"))
      tksee (Envir$txt2,"end")
      tkinsert(Envir$txt2,"end", DataQ)
      tkinsert(Envir$txt2,"end","\n\n")
      tkinsert(Envir$txt2,"end", paste("Years sampled by Categorical factor(s): " ,"\n"))
      tkinsert(Envir$txt2,"end", DataP)
      return(tkinsert(Envir$txt2,"end","\n\n"))
      }
  else{}
 }
#________________________________________________________________________________________preparation du tableau de donnees (TS) par station(s) et parametre
 {                                                          
  if (any(colnames(Envir$Data) == "Category")) {                                                                                          # effectue la commande si une colonne 'Category' existe
      Ts <- subset(Envir$Data, Category %in% site, drop =TRUE) }                                                                          # nouveau tableau avec les Category selectionnees
  else { tkmessageBox(message=paste("No Category column in your dataset")  ,icon = "warning", type = "ok", title="!Warning!") }                                                                           
  
  if (any(colnames(Envir$Data) == "Salinity") & any(!is.na(Envir$Data$Salinity)) == TRUE) {                                               # effectue la commande si une colonne 'Salinity' existe
      if (max(sal) >= round(max(Ts$Salinity, na.rm=TRUE)*2,0)/2 & min(sal) <= round(min(Ts$Salinity, na.rm=TRUE)*2,0)/2 ) { }             # si toutes les salinites sont selectionnees (de min a max) alors toutes les salinites, NA inclus, seront pris en compte.
      else { Sal <- na.omit(Ts$Salinity[(Ts$Salinity >= min(sal)) &  (Ts$Salinity <= max(sal))])                                          # selection des donnees comprises entre les salinites indiquees, NA exclues
             Ts <- subset(Ts, Salinity %in% Sal, drop = TRUE) } }                                                                                                                    # nouveau tableau avec les salinites selectionnees
  else{ }

  if (any(colnames(Envir$Data) == "Depth") & any(!is.na(Envir$Data$Depth)) == TRUE) {                                                     # effectue la commande si une colonne 'Depth' existe
      if (max(depth) >= round(max(Ts$Depth, na.rm=TRUE)*2,0)/2 & min(depth) <= round(min(Ts$Depth, na.rm=TRUE)*2,0)/2 ) { }               # si toutes les profondeurs sont selectionnees (de min a max) alors toutes les profondeurs, NA inclus, seront pris en compte.
      else { Depths <- na.omit(Ts$Depth[(Ts$Depth >= min(depth)) &  (Ts$Depth <= max(depth))])                                            # selection des donnees comprises entre les profondeurs indiquees, NA exclues
             Ts <- subset(Ts, Depth %in% Depths, drop = TRUE) }  }                                                                        # nouveau tableau avec les profondeurs selectionnees
  else{}
  
  if (any(colnames(Envir$Data) == "Depth") & length(Ts$Depth)==0 | any(colnames(Envir$Data) == "Salinity") & length(Ts$Salinity)==0) {
  return(tkmessageBox(message=paste("The selected combinaison of salinity and depth doesn't exist in your dataset")
                               , icon = "warning", type = "ok", title="!Warning!")) }
  else{}

  if (any(colnames(Envir$Data) == "Salinity")) {     
      TS <- data.frame(Ts$Category, Ts$Dates, Ts$Salinity)                                              # nouveau tableau avec les Category selectionnees, les dates et les salinites pour la standardisation des sels nut (les profondeurs sont prise en compte mais n'apparaissent pas)
      TS$param <- Ts[,names(Ts)==param]                                                                 # le parametre selectionne est inclu dans le nouveau tableau
      names(TS) <- c("Category", "Dates", "Salinity", "param") }                                        # les colonnes sont nommees  (le parametre restera "param" pour faciliter les analyses automatiques)
  else { TS <- data.frame(Ts$Category, Ts$Dates)                                                        # meme chose mais s'il n'y a pas de colonne 'Salinity'
      TS$param <- Ts[,names(Ts)==param]                                                       
      names(TS) <- c("Category", "Dates", "param") }                                           
 }
#________________________________________________________________________________________________________Transformer les donnees en log+1 
{ 
 if (log.trans == "YES") { Envir$param <- TS$param
 TS$param <- log10((TS$param)+1)                                                                        # transformation des donnees du parametre
 logtext <- "Log10(x+1)"                                                                                # on renomme le parametre pour les figures
 param <- paste("Log10(x+1)", param) } 
 else{ logtext <- "" }
}
#________________________________________________________________________________________________________Calcul du temps moyen entre 2 donnees en jour
 {
  if (is.numeric(TS$Dates) == TRUE) { }
  else{ TS$Dates <- as.Date(TS$Dates, format="%Y-%m-%d") }                                              # les dates sont formatees
  TS <- TS[!is.na(TS$Dates), ]                                                                          # selection des cas ou les dates sont bien presentes
  if (mix == "YES"){                                                                                    # trie dans l'ordre toutes les mesures si les Category sont mixees
      TS <- TS[order(TS$Dates), ] }
  else {
      TS <- TS[order(TS$Dates), ]
      TS <- TS[order(TS$Category), ]}                                                                   # trie dans l'ordre les mesures par Category si les Category ne sont pas mixees (plus d'actualite)
  TSS <- TS
  ecarts <- TS$Dates[2:nrow(TS)]-TS$Dates[1:(nrow(TS)-1)]                                               # calcul l'ecart en jour entre deux mesures successives
  xx <- as.numeric(ecarts)
  TS$time <- c(1, round(cumsum(xx)))                                                                    # somme cumulee des ecarts et integration au tableau (TS$time) : donne une echelle temporelle a la serie de donnees
  TS$time[TS$time==0] <- 1                                                                    
 }
#_____________________________________________________________________________________________Extraction des dates hors semaines (besoin du time step)
 {
  if (is.numeric(TS$Dates) == TRUE) { TS$YEARS <- TS$Dates }
  else {  
  TS$YEARS <- as.numeric(format(TS$Dates, format = "%Y"))                                               # extraction des annees
  TS$MONTHS <- as.numeric(format(TS$Dates, format = "%m"))                                              # extraction des mois
  TS$days <- as.numeric (format(TS$Dates, format = "%d"))                                               # extraction des jours 01 -> 31
  TS$DY <- as.numeric(format(TS$Dates, format="%j"))                                                    # extraction des jours/annee 1 -> 356
  TS$week.month <- TS$days }
 
  if (is.null(start)) {                                                                                 # selection de l'annee de debut de la serie a analyser
      start <- min(TS$YEARS) } else {}
  if (is.null(end)) {                                                                                   # selection de l'annee de fin de la serie a analyser
      end <- max(TS$YEARS) } else {} 
  TS <- TS[(TS$YEARS >= start & TS$YEARS <= end),]                                                      # tableau avec les donnees comprises entre les annees de debut et de fin
 }  
#______________________________________ Etabli la liste des Category pour les inclure dans les resultats sous forme de texte (3 category maxi ensuite '...') 
{
 if (length(site)==1) liste.stations <- site[1]
 if (length(site)==2) liste.stations <- paste(site[1], site[2], sep=", ")
 if (length(site)==3) liste.stations <- paste(site[1], site[2], site[3], sep=", ")
 if (length(site)>3) { for (i in 1 : 3) { liste.stations <- paste(site[1], site[2], site[3], sep=", ") }
                       liste.stations <- paste(liste.stations, "[...]", sep=", ") }
}
#_______________________________________Statistiques descriptives sur les donnees de base (option : select), prend en compte la selection des annees et mois
 ## meme principe que pour 'rawdata' ##
 {
  if (select == "YES"){
    if (is.numeric(TS$Dates) == TRUE) { } else {
    TS <- subset(TS, MONTHS %in% months, drop = TRUE) }
    Descriptive.Statistics <- summary(TS$param)
    sumdata2 <- capture.output(Descriptive.Statistics)  
    tkinsert(Envir$txt2, "end", paste("-Descriptive statistics on selected site and parameter-"
                                 , "\n", "Parameter: ", param , "\n", "Categorical factor(s): ", liste.stations
                                 , sep="" , "\n\n")) 
    tktag.add(Envir$txt2, "titre", "end -5 lines linestart","end -5 lines lineend")
    tktag.configure(Envir$txt2, "titre", font=tkfont.create(family="courier",size=10,weight="bold"))
    sumdata2 <- paste(sumdata2,collapse="\n")
    tkinsert(Envir$txt2,"end",paste(sumdata2, "\n\n"))
    tkinsert(Envir$txt2,"end", paste("Data quantity (",param,") : ", length(TS$param[!is.na(TS$param)]), "\n\n", sep=""))
    tksee (Envir$txt2,"end")
       
    Descriptive.Statistics <- as.matrix(Descriptive.Statistics)
    Descriptive.Statistics <- as.data.frame(Descriptive.Statistics)
    return(names(Descriptive.Statistics) <- param)
    }
  else{}
 }
#_________________________________________________________________________Aide pour trouver le bon time.step (option : help.timestep et auto.timestep)
 {
#_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _Calcul de l'ecart moyen, en jours, entre 2 mesures successives (les mesures prise le meme jour sont considerees comme 1 seule)
  
  TSS <- TSS[!is.na(TSS$param), ]                                                     # reprend le tableau TS sans les valeurs manquantes
  TSS$YEARS <- as.numeric(format(TSS$Dates, format = "%Y"))                           # extraction des annees
  TSS$MONTHS <- as.numeric(format(TSS$Dates, format = "%m"))                          # extraction des mois
  TSS <- subset(TSS, MONTHS %in% months, drop = TRUE)
  TSS <- TSS[(TSS$YEARS >= start & TSS$YEARS <= end),]                                # tableau avec les donnees comprises entre les annees de debut et de fin
  ecarts <- TSS$Dates[2:nrow(TSS)]-TSS$Dates[1:(nrow(TSS)-1)]                         # calcule les ecarts reels entre chaque mesure presente
  Ecarts <- ecarts[ecarts!=0]                             
  Mean.time <- data.frame(mean(Ecarts))
  mt <- as.numeric(Mean.time)

##_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _Proposition du meilleur time step (option : help.timestep) en fonction de mt
 {  
  if (help.timestep == "NO" & auto.timestep == "NO" & is.null(time.step)) {              
      return() } 
  
  if (time.step != "NULL") {                                                          # affiche la proposition mais effectue l'aggregation avec la methode choisi par l'utilisateur
    
    if (mt<=5) {                                                                      # si l'ecart moyen est inferieur a 5 jours
        T <- "Use daily frequency" }                                                                
     else { if (mt>5 & mt<=10) {                                                      # si l'ecart moyen est entre 5 et 10 jours
          T <- "Use semi-fortnightly frequency"}
      else { if (mt>10 & mt<=24) {                                                    # si l'ecart moyen est entre 10 jours et 24 jours
                 T <- "Use fortnightly frequency" }
             else { if (mt>24 & mt<=60) {                                             # si l'ecart moyen est entre 24 jours et 60 jours
                        T <- "Use monthly frequency" }
                    else { if (mt>60) {                                               # si l'ecart moyen est superieur a 60 jours
                               T <- "Use yearly frequency" }
                           else { return(print("No solution! Should be a problem somewhere...")) }}}}}           
      time.step <- time.step                                               
      tkinsert(Envir$txt, "end", paste("-Time step choice-" , "\n"))
      tktag.add(Envir$txt, "titre", "end -2 lines linestart","end -2 lines lineend")
      tktag.configure(Envir$txt, "titre", font=tkfont.create(family="courier",size=10,weight="bold"))
      tkinsert(Envir$txt, "end", paste("Parameter: ", param , "\n", "Categorical factor(s): ", liste.stations, sep="" , "\n"))  
      tkinsert(Envir$txt, "end", paste("-> Time step manually selected ", "(", time.step, ")", sep="" , "\n"))
      tkinsert(Envir$txt, "end", paste("Advice: ", T, "\n\n")) 
      tksee (Envir$txt,"end") }
                       
###_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _affiche uniquement la proposition
  else { if (help.timestep == "YES") {                                                                           
            Mean.time <- data.frame(rbind("Mean time between two measurements (days):", as.character(round(mean(Ecarts),2))))
            names(Mean.time) <- "T"
            Data.range <- data.frame(rbind("Time range (days):                        ", t(t(as.character(range(Ecarts))))))
            names(Data.range) <- "T"
            Help.timestep <- rbind(Mean.time, Data.range)
 
       if (mt<=5) {
           T <- as.character("-> Use daily frequency             ") }
           else { if (mt>5 & mt<=10) {
                T <- as.character("-> Use semi-fortnightly frequency             ") }
            else { if (mt>10 & mt<=24) {
                       T <- as.character("-> Use fortnightly frequency       ") }
                   else { if (mt>24 & mt<=60) {
                              T <- as.character("-> Use monthly frequency       ") }
                          else { if (mt>60) {
                                     T <- as.character("-> Use yearly frequency           ") }
                                 else {return(print("No solution! Should be a problem somewhere..."))}}}}}
            
            tt <- paste("Mean time between two measurements:", "\n", "  ", as.character(round(mt, 2))," days", "\n",
            "Time range (min - max):", "\n", "  ", as.character(range(Ecarts))[1]," - ", as.character(range(Ecarts))[2], " days", "\n\n",
            "Advice ", T , sep="")

            tkmessageBox(title = "Frequency Choice Guidance", message = tt, icon = "info", type = "ok")
            return()   }
 
###_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _Selection auto du time step (option : auto.timestep)
  else { if (auto.timestep == "YES")  {
          if (mt<=5) {
             Time.step <- "Daily frequency automatically use"
             time.step <- "Daily" }
          else { if (mt>5 & mt<=10) {
                  Time.step <- "Semi-fortnight frequency automatically use"
                  time.step <- "Semi-fortnight" }
             else { if (mt>10 & mt<=24) {
                        Time.step <- "Fortnight frequency automatically use" 
                        time.step <- "Fortnight" }
                    else { if (mt>24 & mt<=60) {
                               Time.step <- "Monthly frequency automatically use"
                               time.step <- "Monthly" }
                           else { if (mt>60) {
                                      Time.step <- "Annual frequency automatically use"
                                      time.step <- "Annual" }
                                  else {return(print("No solution! Should be a problem somewhere..."))}}}}}                

         tkinsert(Envir$txt, "end", paste("-Time step choice-", "\n"))
         tktag.add(Envir$txt, "titre", "end -2 lines linestart","end -2 lines lineend")                
         tktag.configure(Envir$txt, "titre", font=tkfont.create(family="courier",size=10,weight="bold"))
         tkinsert(Envir$txt, "end", paste("Parameter: ", param , "\n", "Categorical factor(s): ", liste.stations, sep="" , "\n"))
         tkinsert(Envir$txt, "end", paste(Time.step, "\n\n"))
         tksee (Envir$txt,"end")             
 }
  else {}  
 } } } 
 }
#_____________________________________________________________________________________________________Extraction des semaines en fonction du time.step
 {
  if (is.numeric(TS$Dates) == TRUE) { } else { 
  if (time.step == "Fortnight") {                                                                                                                                                     
      TS$week.month[TS$days <= 15] <- 1                                              # extraction du numero de quinzaine du mois (echelle hebdo = 2 quinzaines/mois)
      TS$week.month[TS$days > 15] <- 2 }                                                                                                                                           
  else {                                                                                                                                                                           
      TS$week.month[TS$days <= 8] <- 1                                               # extraction du numero de semi-quinzaine du mois (autres echelle = 4 semi-quinzaines/mois)
      TS$week.month[TS$days > 8 & TS$days <= 15] <- 2                                                                                                                              
      TS$week.month[TS$days > 15 & TS$days <= 23] <- 3                                                                                                                             
      TS$week.month[TS$days > 23] <- 4 }                                                                                                                                           
  TS$days <- NULL                                                                    # suppression de la colonne days pour plus de clarete (plus besoin)
  TSs <- TS
 } }
#_____________________________Valeurs normalisees en fonction de la salinite (npsu), pour les sels nutritifs. Effectue sur les donnees non regularisees
 {  
  if (norm == "YES") {
     if (any(colnames(Envir$Data) == "Salinity")) { 
      if (local.trend == "YES") { return(tkmessageBox(message="Work only for Seasonal and Global Trend", icon = "warning", type = "ok", title="Warning")) }
      TS <- subset(TS, MONTHS %in% months, drop = TRUE)                              # selection des mois a traiter
      dir.create(paste(Envir$save.WD,"/",Envir$File.Name,"/", liste.stations, "/", start,"-", end, "/", strsplit(param,"/")[[1]][1], "/", sep= ""), recursive = TRUE, showWarnings = FALSE)  # creation du repertoir de sauvegarde
      
      a <- min(TS$YEARS, na.rm=TRUE)
      b <- max(TS$YEARS, na.rm=TRUE)
      c <- sort(rep(a:b, length(months)))
      NormNutri <- matrix(c)                                                         # matrice vide que l'on va remplir avec les valeurs normalisees par annees (a:b)
      NormNutri <- data.frame(cbind(NormNutri, (rep(months, length(a:b)))))
      NormNutri <- data.frame(cbind(NormNutri, (rep(NA, length(c)))))
      names(NormNutri) <- c("YEARS", "MONTHS", "Normalized")                         # on nomme les colonnes par souci de clarete
      
      for (i in a:b) {                                                               # pour chaque annee i :
       for (j in months) {                                                           # et pour chaque mois j :
           nutrient <- TS$param[TS$YEARS==i & TS$MONTHS==j]                          #      on recupere les valeurs du parametre a traiter
           sal <- TS$Salinity[TS$YEARS==i & TS$MONTHS==j]                            #      on recupere les salinites a traiter
 
           if ( sum(nutrient*sal, na.rm=TRUE)== 0 | max(sal, na.rm=TRUE)<npsu ) { 
           NormNutri$Normalized[NormNutri$YEARS==i & NormNutri$MONTHS==j] <- NA }                                                 # si toute les conditions sont rempli (parametres present et salinites inferieures a la salinite de normalisation) :
           else { reg <- lm(nutrient~sal)                                                                                                                                #           on effectue la regression entre salinite et nutriment par annee
           NormNutri$Normalized[NormNutri$YEARS==i & NormNutri$MONTHS==j] <- (npsu * (reg$coefficients[2]))+(reg$coefficients[1])          #          et on calcule la valeur normalisee
            
            save.mixingdia.path <- paste(Envir$save.WD,"/",Envir$File.Name,"/", liste.stations, "/", start,"-",end,"/", strsplit(param,"/")[[1]][1], "/", Envir$File.Name,"_MixingDiagram_", strsplit(param,"/")[[1]][1],"_",j,",",i,".png", sep="")                        # sauve chaque mixing diagram (par annee)
            if (nchar(save.mixingdia.path)>259) { return(tkmessageBox(message= paste("The save path is too long (NTFS system limit)", "\n", "Your results cannot be saved properly", "\n", "Please consider shorter parameter and/or station name", sep=""), icon = "warning", type = "ok", title="!Warning!")) } else { }
            png(save.mixingdia.path)
            plot(nutrient~sal, xlab="Salinity", ylab=param, main=paste("Mixing diagram of ", param, " in ", j, "/", i, sep=""))
             if( !is.na(reg$coefficients[2]) ) { abline(lm(nutrient~sal)) }                                                        # ajoute la droite de regression si elle est realisable
             else{ }
            dev.off() }}}        
           
           save.mixingreg.path <- paste(Envir$save.WD,"/",Envir$File.Name,"/", liste.stations, "/", start,"-",end,"/", strsplit(param,"/")[[1]][1], "/", Envir$File.Name,"_Normalized_", strsplit(param,"/")[[1]][1],"_at_",npsu,".txt", sep="")
           if (nchar(save.mixingreg.path)>259) { return(tkmessageBox(message= paste("The save path is too long (NTFS system limit)", "\n", "Your results cannot be saved properly", "\n", "Please consider shorter parameter and/or station name", sep=""), icon = "warning", type = "ok", title="!Warning!")) } else { }
           write.table(NormNutri, sep="\t", row.names=FALSE, file=save.mixingreg.path)
 
           if (sum(NormNutri$Normalized, na.rm=TRUE)== 0 | max(TS$Salinity, na.rm=TRUE)<npsu ) {                                          # message d'avertissement si les conditions ne sont pas rempli
           return(tkmessageBox(message=paste("No parameters available or all the salinities in your database are under the 'selected psu'"
                               , "\n", "Try again by decreasing the 'selected psu' value", sep="")
                               , icon = "warning", type = "ok", title="!Warning!")) }
      
      z <- ts(NormNutri$Normalized, start = a, deltat=1/(length(months)))
      dev.new()
      plot(z, ylab="", xlab="", type = "o", pch = 20, col="aquamarine3", axes=FALSE, lwd=1.5, main=paste("Trend of", param,"at salinity", npsu))
      title(main=paste("\n\n", "Categorical factor(s): ", liste.stations), cex.main=0.8)
      box(which = "plot", col="steelblue2", lwd=2)
      axis(side=1, col="steelblue2", col.ticks="black", tck=-0.02, lwd=2)
      axis(side=1, at = seq(start(z)[1], end(z)[1]+1, by = 1), col="steelblue2", col.ticks="black", labels = FALSE, tck= -0.01, lwd=2)
      axis(side=2, col="steelblue2", col.ticks="black", tck=-0.02, lwd=2)
      axis(side=2, at = seq(par()$yaxp[1], par()$yaxp[2], by = ((par()$yaxp[2]-par()$yaxp[1])/(par()$yaxp[3]*2)))
           , labels=FALSE, col="steelblue2", col.ticks="black", tck=-0.01, lwd=2)
      mtext(paste0("                                                                           Years (at selected period)"), side=1, line=2.5, font=1, cex=1.2, col="black")
      mtext("Normalized concentration", side=2, line=2.5, font=1, cex=1.2, col="black")
           save.normnutri.path <- paste(Envir$save.WD,"/",Envir$File.Name,"/", liste.stations, "/", start,"-",end,"/", strsplit(param,"/")[[1]][1],"/", Envir$File.Name, "_NormalNutri_", strsplit(param,"/")[[1]][1],".png", sep = "")
            if (nchar(save.normnutri.path)>259) { return(tkmessageBox(message= paste("The save path is too long (NTFS system limit)", "\n", "Your results cannot be saved properly", "\n", "Please consider shorter parameter and/or station name", sep=""), icon = "warning", type = "ok", title="!Warning!")) } else { }
             png(save.normnutri.path, width =1300, height=700, res=100)
              plot(z, ylab="", xlab="", type = "o", pch = 20, col="aquamarine3", axes=FALSE, lwd=1.5, main=paste("Trend of", param,"at salinity", npsu))
              title(main=paste("\n\n", "Categorical factor(s): ", liste.stations), cex.main=0.8)
              box(which = "plot", col="steelblue2", lwd=2)
              axis(side=1, col="steelblue2", col.ticks="black", tck=-0.02, lwd=2)
              axis(side=1, at = seq(start(z)[1], end(z)[1]+1, by = 1), col="steelblue2", col.ticks="black", labels = FALSE, tck= -0.01, lwd=2)
              axis(side=2, col="steelblue2", col.ticks="black", tck=-0.02, lwd=2)
              axis(side=2, at = seq(par()$yaxp[1], par()$yaxp[2], by = ((par()$yaxp[2]-par()$yaxp[1])/(par()$yaxp[3]*2)))
                   , labels=FALSE, col="steelblue2", col.ticks="black", tck=-0.01, lwd=2)
              mtext(paste0("Years (at selected period)"), side=1, line=2.5, font=1, cex=1.2, col="black")
              mtext("Normalized concentration", side=2, line=2.5, font=1, cex=1.2, col="black")
             dev.off()
      tkinsert(Envir$txt2,"end", paste("-Global trend of",param,"at salinity", npsu,"-", "\n"))
      tktag.add(Envir$txt2, "titre", "end -2 lines linestart","end -2 lines lineend")
      tktag.configure(Envir$txt2, "titre", font=tkfont.create(family="courier",size=10,weight="bold"))
      tkinsert(Envir$txt2, "end", paste("Categorical factor(s): ", liste.stations, sep="" , "\n\n"))
      sK <- seaKen(z)
      save.normtrend.path <- paste(Envir$save.WD,"/",Envir$File.Name,"/", liste.stations, "/", start,"-",end,"/", strsplit(param,"/")[[1]][1],"/", Envir$File.Name, "_NormGlobalTrend_", strsplit(param,"/")[[1]][1],".txt", sep = "")
      if (nchar(save.normtrend.path)>259) { return(tkmessageBox(message= paste("The save path is too long (NTFS system limit)", "\n", "Your results cannot be saved properly", "\n", "Please consider shorter parameter and/or station name", sep=""), icon = "warning", type = "ok", title="!Warning!")) } else { }
      write.table(sK, sep="\t", row.names=FALSE, file=save.normtrend.path) 
      tkinsert(Envir$txt2,"end", paste("Trend (sen.slope): ", round(sK$sen.slope, 4), "  original units per year", "\n"
                                 , "%Trend (sen.slope.pct): ", round(sK$sen.slope.pct, 4), "  percent of mean quantity per year", "\n"  
                                 , "p.value: ", round(sK$p.value, 4), "\n\n\n", sep=""))
      if (sK$p.value <= 0.05) { tktag.add(Envir$txt2, "titre5", "end -4 lines linestart","end -4 lines lineend") 
      tktag.configure(Envir$txt2, "titre5", font=tkfont.create(family="courier",size=10,weight="bold"))} else {}
      tksee (Envir$txt2,"end")  
      return(showData(capture.output(sK), rowname.bar = NULL, colname.bar = NULL, title="Global Trend (Mann Kendall Test) results"))
 }
  else{ return(tkmessageBox(message="No salinity data available", icon = "warning", type = "ok", title="Warning")) } }
  else{ } 
 } 
#______________________________________________________________________________Suppression de la colonne salinite maintenant que l'on en a plus besoin
 { 
 if(any(colnames(TS) == "Salinity")) { TS$Salinity <- NULL }
 else { }                                   
 }
#_____________________________________________________________Boxplot et Extraction des outliers par annee: pas necessaire pour un test Mann-Kendall
 {   
  if (plotB == "YES") {
      if (is.numeric(TS$Dates) == TRUE) { } else {
      TS <- subset(TS, MONTHS %in% months, drop = TRUE) }
      dir.create(paste(Envir$save.WD,"/",Envir$File.Name,"/", liste.stations, "/", start,"-", end, "/", strsplit(param,"/")[[1]][1], "/", sep= ""), recursive = TRUE, showWarnings = FALSE)
     if (selectBox=="ByYears") {
      save.boxplot.path <- paste(Envir$save.WD,"/",Envir$File.Name,"/", liste.stations, "/", start,"-", end, "/", strsplit(param,"/")[[1]][1], "/", Envir$File.Name,"_Boxplot_byYears_", strsplit(param,"/")[[1]][1],".png", sep = "") 
      if (nchar(save.boxplot.path)>259) { return(tkmessageBox(message= paste("The save path is too long (NTFS system limit)", "\n", "Your results cannot be saved properly", "\n", "Please consider shorter parameter and/or station name", sep=""), icon = "warning", type = "ok", title="!Warning!")) } else { }
      png(save.boxplot.path, width=1000, height=1000, res=150)
      boxplot(TS$param~TS$YEARS, xlab="YEARS", ylab=paste(param) , main=paste("Boxplot of",param,"(o = outliers)"))
      grid(lwd=0.3, lty="dotted")
      title(main=paste("\n\n\n", "Categorical factor(s): ", liste.stations), cex.main=0.8) 
      dev.off()
      dev.new(width=8, height=5)
      boxplot(TS$param~TS$YEARS, xlab="YEARS", ylab=paste(param), col="grey90", main=paste("Boxplot of",param,"(o = outliers)"))                                     # affiche la boite a moustache par annees
      title(main=paste("\n\n\n", "Station(s): ", liste.stations), cex.main=0.8)
      return() }
     if (selectBox=="ByMonths") {
      save.boxplot.path <- paste(Envir$save.WD,"/",Envir$File.Name,"/", liste.stations, "/", start,"-", end, "/", strsplit(param,"/")[[1]][1], "/", Envir$File.Name,"_Boxplot_byMonths_", strsplit(param,"/")[[1]][1],".png", sep = "") 
      if (nchar(save.boxplot.path)>259) { return(tkmessageBox(message= paste("The save path is too long (NTFS system limit)", "\n", "Your results cannot be saved properly", "\n", "Please consider shorter parameter and/or station name", sep=""), icon = "warning", type = "ok", title="!Warning!")) } else { }
      png(save.boxplot.path, width=1000, height=1000, res=150)
      boxplot(TS$param~TS$MONTHS, xlab="MONTHS", ylab=paste(param), main=paste("Boxplot of",param,"(o = outliers)"))        
      title(main=paste("\n\n\n", "Categorical factor(s): ", liste.stations), cex.main=0.8)
      dev.off()
      dev.new(width=8, height=5)
      boxplot(TS$param~TS$MONTHS, xaxt="n", xlab="MONTHS", ylab=paste(param), axes=F, col="grey90", main=paste("Boxplot of",param,"(o = outliers)"))                 # affiche la boite a moustache par mois
      title(main=paste("\n\n\n", "Station(s): ", liste.stations, "   Years: ", start,"-", end), cex.main=0.8)
      box(which = "plot")
      axis(side=1, at =(1:12), labels = month.abb)
      axis(side=2)
      return() } }
  else{}

  if (outliers.re == "YES") {
      a <- min(TS$YEARS)                                                              # pose les arguments qui seront utilises
      b <- max(TS$YEARS)
      Tf <- NULL
      Text <- NULL

      for (i in a:b) {                                                                # extraction des donnees tous les ans de la premiere (a) a la derniere (b)
           T <- TS$param[TS$YEARS==i]                                                 # extrait les donnees a traite pour l'annee i
           EXT <- TS$param[TS$YEARS==i]                                          
           Q3 <- quantile (T, 0.75, na.rm=TRUE)                                       # calcul le 3eme quantile
           Q1 <- quantile (T, 0.25, na.rm=TRUE)                                       # calcul le 1er quantile
           Q <- Q3-Q1
           h <- Q3+(Q*1.5)                                                            # limite haute des outliers
           l <- Q1-(Q*1.5)                                                            # limite basse des outliers
           T[T>h]<-NA                                                                 # suppression des donnees au dessus de la limite haute
           T[T<l]<-NA                                                                 # suppression des donnees en dessous de la limite basse
           EXT[EXT<=h & EXT>=l]<-NA                                                   # conserve les donnees supprimees dans EXT

           t <- list(T)                                                                    
           Tf <- rbind.data.frame(Tf, t)                                              # construit le tableau (Tf) annee par annee (donnees avec outliers supprimes)
           ext <- list(EXT)                                                                 
           Text <- rbind.data.frame(Text, ext) }                                      # construit le tableau (Text) des donnees supprimees annee apres annee
                                                                                 
      tf <- Tf[ ,1]
      names(tf) <- c("param")
      TS$param <- tf                                                                  # remplace les donnees brute du tableau TS par les donnees (Tf) avec outliers retires
      Outliers <- data.frame(TS$Category, Text[ ,1], TS$Dates)
      names(Outliers) <- c("Category", "Outliers", "Dates")
      save.outl.path <- paste(Envir$save.WD,"/",Envir$File.Name,"/", liste.stations, "/", start,"-", end, "/", strsplit(param,"/")[[1]][1], "/", Envir$File.Name, "_Outliers_", strsplit(param,"/")[[1]][1],".txt", sep = "")
      if (nchar(save.outl.path)>259) { return(tkmessageBox(message= paste("The save path is too long (NTFS system limit)", "\n", "Your results cannot be saved properly", "\n", "Please consider shorter parameter and/or station name", sep=""), icon = "warning", type = "ok", title="!Warning!")) } else { }
      dir.create(paste(Envir$save.WD,"/",Envir$File.Name,"/", liste.stations, "/", start,"-", end, "/", strsplit(param,"/")[[1]][1], "/", sep = ""), recursive = TRUE, showWarnings = FALSE)
      write.table(Outliers, sep="\t", row.names=FALSE, file=save.outl.path) }         # sauve les donnees supprimees dans un fichier csv
         
  else{}
 }
#____________________________________________________________________________________________Supprime la colonne TS$DATES dans TS pour plus de clarete 
 {
  Dates <- TS$Dates
  TS$Dates <- NULL                                                                               
 }
#_________________________________________________________________________Ajout d'un tableau temporel fictif pour combler les mois/semaines manquants
 {
#_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ Tableau des annees

  yy <- (min(TS$YEARS):max(TS$YEARS))                                                 # suite de toutes les annees comprise entre l'annee de debut et de fin de la serie
  mm <- rep(1, length(yy))                                                            # remplis la nouvelle colonne mois avec des 1
  aa <- rep(NA, length(yy))
  if (is.numeric(Dates) == TRUE) {
  fic2 <- data.frame(aa,aa,aa,yy)                                                     # nouveau tableau fictif
  names(fic2) = c("Category", "param", "time", "YEARS") } else {                      # ...
  
  fic2 <- data.frame(aa,aa,aa,yy,mm,aa,aa)                                            # nouveau tableau fictif
  names(fic2) = c("Category", "param", "time", "YEARS", "MONTHS", "week.month", "DY")  }                

  TS <- rbind(fic2, TS)                                                               # fusion des tableaux
  
  if (is.numeric(Dates) ==TRUE) { } else {
  if (time.step == "Daily") {
    ddd <- c(1:366)
    www <- rep(NA, 366)
    
  fic3 <- data.frame(www,www,www,www,www,www,ddd)
  names(fic3) = c("Category", "param", "time", "YEARS", "MONTHS", "week.month", "DY")

  TS <- rbind(fic3, TS) } }
 
#_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ Tableau des mois et semaines (meme principe, mais 24 semaines pour semi-hebdomadaire et 48 pour les autres)

  if (is.numeric(Dates) ==TRUE) { } else {
  if (time.step == "Fortnight") {

      w <- rep(1:2, 12)
      m <- gl(12, 2)
      y <- rep(-1000, 24)
      a <- rep(NA, 24)
      b <- rep(NA, 24)
      c <- rep(NA, 24)}
 
  else {
      w <- rep(1:4, 12)
      m <- gl(12, 4)
      y <- rep(-1000, 48)
      a <- rep(NA, 48)
      b <- rep(NA, 48)
      c <- rep(NA, 48) }
 
  fic1 <- data.frame(a,b,c,y,m,w,a)
  names(fic1) = c("Category", "param", "time", "YEARS", "MONTHS", "week.month", "DY")

  TS <- rbind(fic1, TS)  }
 }
#________________________________________________________________ Tableaux croises avec toutes les methodes d'aggregation (temporelle et mathematique)
 {
 if(time.step == "Daily") {  name.freq <- c("Days")
        Mean <- tapply(TS$param, list(TS$YEARS, TS$DY), mean, na.rm=TRUE)
        Quantile <- tapply(TS$param, list(TS$YEARS, TS$DY), quantile, probs = (0.9), na.rm=TRUE)
        Median <- tapply(TS$param, list(TS$YEARS, TS$DY), median, na.rm=TRUE)
        Max <- tapply(TS$param, list(TS$YEARS, TS$DY), max, na.rm=TRUE) }
  else { if (time.step == "Fortnight") { name.freq <- c("Fortnights")
          Mean <- tapply(TS$param, list(TS$YEARS, TS$week.month, TS$MONTHS), mean, na.rm=TRUE)
          Quantile <- tapply(TS$param, list(TS$YEARS, TS$week.month, TS$MONTHS), quantile, probs = (0.9), na.rm=TRUE)
          Median <- tapply(TS$param, list(TS$YEARS, TS$week.month, TS$MONTHS), median, na.rm=TRUE)
          Max <- tapply(TS$param, list(TS$YEARS, TS$week.month, TS$MONTHS), max, na.rm=TRUE) }
     else { if (time.step=="Monthly") {  name.freq <- c("Months")
                 Mean <- tapply(TS$param, list(TS$YEARS, TS$MONTHS), mean, na.rm=TRUE)
                 Quantile <- tapply(TS$param, list(TS$YEARS, TS$MONTHS), quantile, probs = (0.9), na.rm=TRUE)
                 Median <- tapply(TS$param, list(TS$YEARS, TS$MONTHS), median, na.rm=TRUE)
                 Max <- tapply(TS$param, list(TS$YEARS, TS$MONTHS), max, na.rm=TRUE) }
         else { if (time.step=="Annual") {   name.freq <- c("Years")
                        Mean <- tapply(TS$param, list(TS$YEARS), mean, na.rm=TRUE)
                        Quantile <- tapply(TS$param, list(TS$YEARS), quantile, probs = (0.9), na.rm=TRUE)
                        Median <- tapply(TS$param, list(TS$YEARS), median, na.rm=TRUE)
                        Max <- tapply(TS$param, list(TS$YEARS), max, na.rm=TRUE)}
                else { if (time.step == "Mono-mensual") {  name.freq <- c("Months")
                               Mean <- tapply(TS$param, list(TS$MONTHS), mean, na.rm=TRUE)
                               Quantile <- tapply(TS$param, list(TS$MONTHS), quantile, probs = (0.9), na.rm=TRUE)
                               Median <- tapply(TS$param, list(TS$MONTHS), median, na.rm=TRUE)
                               Max <- tapply(TS$param, list(TS$MONTHS), max, na.rm=TRUE) }
                       else { if (time.step == "Semi-fortnight") {   name.freq <- c("Semi-fortnights")
                                      Mean <- tapply(TS$param, list(TS$YEARS, TS$week.month, TS$MONTHS), mean, na.rm=TRUE)
                                      Quantile <- tapply(TS$param, list(TS$YEARS, TS$week.month, TS$MONTHS), quantile, probs = (0.9), na.rm=TRUE)
                                      Median <- tapply(TS$param, list(TS$YEARS, TS$week.month, TS$MONTHS), median, na.rm=TRUE)
                                      Max <- tapply(TS$param, list(TS$YEARS, TS$week.month, TS$MONTHS), max, na.rm=TRUE) }     
                              else {return(print("You have to choose a time step to aggregate your data!"))}                                
  }}}}}
 }  
#__________________________________________Choix de la methode d'aggregation avec les p.values de l'ANOVA suivi d'un post-hoc de Dunett (options : help.aggreg et auto.aggreg)
 {
  if (help.aggreg == "YES" | auto.aggreg== "YES") {
 
      a <- data.frame(rep(NA, (max(TS$YEARS, na.rm=TRUE)+1-min(TS$YEARS[TS$YEARS > 0], na.rm=TRUE))*365.25))        # tableau vide avec nombre de jour total de la serie regularisee
      a <- data.frame(a, 1:nrow(a))                                                            # numerotation des lignes
      aa <- data.frame(TS$param[!is.na(TS$time)])                                              # on recupere les donnees parametres
      aa <- data.frame(aa, TS$time[!is.na(TS$time)])                                           #     et l'on y accole les donnees temporelles
      regraw <- merge(a, aa, by.x = 2, by.y = 2, all = TRUE)                                   # on fusionne les donnees brute (aa) et la serie temporelle
      zraw <- ts(regraw[, 3], start=min(TS$YEARS[TS$YEARS > 0], na.rm=TRUE), frequency=366)    #     reguliere pour obtenir une serie regularisee (echelle de la journee) avec des NA
   
      v1 <- as.data.frame(Mean)                                                                # on recupere le tableau croise
      v2 <- as.data.frame(Quantile)
      v3 <- as.data.frame(Median)
      v4 <- as.data.frame(Max)
      v4[v4=="-Inf"]<-NA

      MEA <- reshape(v1, direction = "long", varying = list(1: ncol(v1)), times = colnames(v1))    # on le decroise
      QUA <- reshape(v2, direction = "long", varying = list(1: ncol(v2)), times = colnames(v2))
      MED <- reshape(v3, direction = "long", varying = list(1: ncol(v3)), times = colnames(v3))
      MAX <- reshape(v4, direction = "long", varying = list(1: ncol(v4)), times = colnames(v4))
      
      score <- c(regraw[, 3], MEA[,2], QUA[,2],  MED[,2], MAX[,2])
      method <- c(rep(1, length(regraw[, 3])), rep(2, length(MEA[,2])), rep(3, length(QUA[,2])), rep(4, length(MED[,2])), rep(5, length(MAX[,2])))
      method <- factor(method)
      
      anova_results <- aov(score~method)
      test.dunnett <- glht(anova_results,linfct=mcp(method="Dunnett"))
      
      C1 <- "The means give better fit"
      C2 <- "The quantiles 0.9 give better fit"
      C3 <- "The medians give better fit"
      C4 <- "The maximums give better fit"
      choice <- rbind(C1,C2,C3,C4)                                                            # tableau des differentes propositions

      R <- cbind(data.frame(summary(test.dunnett)$test$pvalues),choice)
      r <- as.data.frame(R)                                                                   # fusion des tableaux
      names(r)[c(1)] <- c("p.values")                                                         # on fignole le tableau...
      names(r)[c(2)] <- c("choice") 
      r$choice <- as.character(r$choice)                                                                   
      r$p.values <- as.numeric(r$p.values) }
  else{}

###_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ Proposition de la methode d'aggregation (option : help.aggreg)
 
  if (help.aggreg=="YES") {
      tkmessageBox(title = "Method of Aggregation Guidance", message = as.character(subset(r$choice, r$p.values==max(r$p.values)), icon = "info", type = "ok"))
      return() }
  else{ }
 
###_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ Choix automatique de la methode d'aggregation (option :auto.aggreg)
 
  if (auto.aggreg == "YES") {
     if (which.max(r$p.values) == 1) { 
         aggreg <- "Mean" 
         Aggreg <- "Method of aggregation automatically use : mean" }
     else { if (which.max(r$p.values) == 2) { 
                aggreg <- "Quantile" 
                Aggreg <- "Method of aggregation automatically use : quantile(0.9)" }   
            else { if (which.max(r$p.values) == 3) { 
                       aggreg <- "Median" 
                       Aggreg <- "Method of aggregation automatically use : median" }
                   else { if (which.max(r$p.values) == 4) { 
                              aggreg <- "Max" 
                              Aggreg <- "Method of aggregation automatically use : maximum" }}}}
    
  tkinsert(Envir$txt,"end", paste("-Method of aggregation choice-","\n"))
  tktag.add(Envir$txt, "titre", "end -2 lines linestart","end -2 lines lineend")
  tkinsert(Envir$txt, "end", paste("Parameter: ", param, "\n", "Categorical factor(s): ", liste.stations, sep="" , "\n")) 
  tkinsert(Envir$txt,"end", paste(Aggreg, "\n\n"))
  tksee (Envir$txt,"end")
           
 }} 
#_____________________________________________________________________________________________Remplacement des valeurs manquantes, option (na.replace) 
 {   
  v <- as.data.frame(get(aggreg))
  v[v=="-Inf"]<- NA
  v[v=="NaN"] <- NA

  if (na.replace=="YES") {
     wx <- v[-1, ]
     test.wx <- as.data.frame(wx)                                                                                                                                # evite les problemes de format "array" pour les frequences annuelle et monomensuelle
     P.NA <- (length(test.wx[is.na(test.wx)])/(dim(test.wx)[1]*dim(test.wx)[2]))*100
     P.NA <- round(P.NA,2)
     
     if ( length(test.wx[is.na(test.wx)]) > ((dim(test.wx)[1]*dim(test.wx)[2])/20) ) {                          # affiche un message de warning si les NAs representent plus de 5% de la serie
         tkinsert(Envir$txt, "end", paste("-!Warning message!-", "\n"))
         tktag.add(Envir$txt, "titre2", "end -2 lines linestart","end -2 lines lineend") 
         tktag.configure(Envir$txt, "titre2", foreground="red", font=tkfont.create(family="courier",size=10,weight="bold"))
         tkinsert(Envir$txt, "end", paste("Missing values represent ", P.NA, "% of the regularised data,", "\n", sep=""))
         tkinsert(Envir$txt, "end", paste("replacing them is not a good idea.", "\n\n"))
         tksee (Envir$txt,"end") 
         tkmessageBox(message=paste("Missing values represent ", P.NA, "% of the regularised data,", "\n", "replacing them is not a good idea.", sep=""), icon = "warning", type = "ok", title="!Warning!") }
     else{} 
     
     cc <- matrix(NA, nrow(v), ncol(v))                                              # matrice vide a remplir
     vc <- v
     vc[is.na(vc)] <- -5000                                                          # remplacement des NA par la valeur -5000 (evite les problemes avec les NA et les 0)
     ccc <- impute(v, what=c("median"))                                              # matrice avec valeurs remplacees par la methode des medianes/cycle # impute = package 'e1071' ##
     
     if (time.step=="Annual" | time.step=="Mono-mensual") {
        for (i in 1:(nrow(v)-2)) {
           if ((vc[i,]+vc[i+1,]+vc[i+2,])== -15000) {                                # pour chaque serie de 3 valeurs successives dans vc, si leur somme = -15000 (valeurs manquantes) alors
                cc[i,] <- ccc[i,] }                                                  # la premiere donnee de cette serie est remplacee par la valeur correspondante dans la matrice des medianes
           else { cc[i,] <- v[i,] } } 
        cc[nrow(v)-1, ] <- v[nrow(v)-1, ]
        cc[nrow(v), ] <- v[nrow(v), ] } 
     
     else {                                                                        
        for (i in 1:nrow(v)) {
           for (j in 1:(ncol(v)-2)) {
              if ((vc[i,j]+vc[i,j+1]+vc[i,j+2])== -15000) {                          # pour chaque serie de 3 valeurs successives dans vc, si leur somme = -15000 (valeurs manquantes) alors
                  cc[i,j] <- ccc[i,j] }                                              # la premiere donnee de cette serie est remplacee par la valeur correspondante dans la matrice des medianes
              else { cc[i,j] <- v[i,j] } } }                                         # sinon la valeur est laissee telle quelle (meme si c'est une valeurs manquante)
        cc[ , ncol(v)-1] <- v[ , ncol(v)-1]
        cc[ , ncol(v)] <- v[ , ncol(v)]    }
  
     if (is.na(cc[1,1])) {                                                           # remplace la premiere valeur de la serie si celle-ci est toujours NA
         cc[2,1] <- ccc[2,1] }                            
     else{}                                             
     if (is.na(cc[nrow(cc),ncol(cc)])) {                
         cc[nrow(cc),ncol(cc)] <- ccc[nrow(cc),ncol(cc)] }
     else{} 
       
     row.names(cc) <- row.names(v)
     v <- as.data.frame(cc) 
     names(v) <- names(vc) }
##### autre methode possible: les splines ->  spline(v, n=length(v), method="natural")  , mais moyen si beaucoup de valeur manquante  #####
  
  else { v <- as.data.frame(v) }
  
 }
#_____________________________________________________________________________Re-formation du tableau de donnees regularisees en fonction du time.step
 {
  TimeSerie <- reshape(v, direction = "long", varying = list(1: ncol(v)), times = colnames(v))        # Transpose le tableau croise en tableau colonne
  
  if (time.step == "Daily"){
  
  TimeSerie$DayYears <- do.call("rbind", strsplit(TimeSerie$time,"\\.")) [,1]
  TimeSerie$DayYears <- as.numeric(TimeSerie$DayYears)
   
  TimeSerie$time <- NULL
  names(TimeSerie)[c(1)] <- c("param")
  names(TimeSerie)[c(2)] <- c("YEARS")
  Y <- as.numeric(row.names(v))
  TimeSerie$YEARS <- c(Y)
  TimeSerie <- TimeSerie[order(TimeSerie$YEARS, TimeSerie$DayYears), ]
  TimeSerie <- TimeSerie[TimeSerie$YEARS!=(-1000), ]
  row.names(TimeSerie) <- 1:(nrow(TimeSerie))
  TimeSerie$time <- 1:(nrow(TimeSerie))
  
  Ja <- as.matrix(rep(1, 31))
  Fe <- as.matrix(rep(2, 28))
  Fe.leap <- as.matrix(rep(2, 29))
  Ma <- as.matrix(rep(3, 31))
  Av <- as.matrix(rep(4, 30))
  My <- as.matrix(rep(5, 31))
  Ju <- as.matrix(rep(6, 30))
  Jl <- as.matrix(rep(7, 31))
  Ao <- as.matrix(rep(8, 31))
  Se <- as.matrix(rep(9, 30))
  Oc <- as.matrix(rep(10, 31))
  No <- as.matrix(rep(11, 30))
  De <- as.matrix(rep(12, 31))
  Mois <- rbind(Ja,Fe,Ma,Av,My,Ju,Jl,Ao,Se,Oc,No,De,NA)
  Mois.leap <- rbind(Ja,Fe.leap,Ma,Av,My,Ju,Jl,Ao,Se,Oc,No,De)

  for (i in yy)  {  if (leap_year(i) == TRUE ) {
   TimeSerie$MONTHS [TimeSerie$YEARS ==i] <- Mois.leap }
   else { TimeSerie$MONTHS [TimeSerie$YEARS==i]  <- Mois } }
  
  TimeSerie <- subset(TimeSerie, select=c(param, YEARS, MONTHS, DayYears, time)) } 
  
 else { if (time.step=="Annual"){

   TimeSerie$time <- NULL
   names(TimeSerie)[c(1)] <- c("param")
   names(TimeSerie)[c(2)] <- c("YEARS")
   Y <- as.numeric(row.names(v))
   TimeSerie$YEARS <- c(Y)
   TimeSerie <- TimeSerie[order(TimeSerie$YEARS), ]
   TimeSerie <- TimeSerie[TimeSerie$YEARS!=(-1000), ]
   row.names(TimeSerie) <- 1:(nrow(TimeSerie))
   TimeSerie$time <- 1:(nrow(TimeSerie)) }
 
  else { if (time.step == "Monthly"){
          TimeSerie$MONTHS <- do.call("rbind", strsplit(TimeSerie$time,"\\.")) [,1]
          TimeSerie$MONTHS <- as.numeric(TimeSerie$MONTHS)
          TimeSerie$time <- NULL
          names(TimeSerie)[c(1)] <- c("param")
          names(TimeSerie)[c(2)] <- c("YEARS")
          Y <- as.numeric(row.names(v))
          TimeSerie$YEARS <- c(Y)
          TimeSerie <- TimeSerie[order(TimeSerie$YEARS, TimeSerie$MONTHS), ]
          TimeSerie <- TimeSerie[TimeSerie$YEARS!=(-1000), ]
          row.names(TimeSerie) <- 1:(nrow(TimeSerie))
          TimeSerie$time <- 1:(nrow(TimeSerie))
          coefT <- 1 }
          
              else { if (time.step=="Fortnight"){
                     TimeSerie$MONTHS <- do.call("rbind", strsplit(TimeSerie$time,"\\.")) [,2]
                     TimeSerie$week.month <- do.call("rbind", strsplit(TimeSerie$time,"\\.")) [,1]
                     TimeSerie$MONTHS <- as.numeric(TimeSerie$MONTHS)
                     TimeSerie$week.month <- as.numeric(TimeSerie$week.month)
                     TimeSerie$time <- NULL
                     names(TimeSerie)[c(1)] <- c("param")
                     names(TimeSerie)[c(2)] <- c("YEARS")
                     Y <- as.numeric(row.names(v))
                     TimeSerie$YEARS <- c(Y)
                     TimeSerie <- TimeSerie[order(TimeSerie$YEARS, TimeSerie$MONTHS), ]
                     TimeSerie <- TimeSerie[TimeSerie$YEARS!=(-1000), ]
                     row.names(TimeSerie) <- 1:(nrow(TimeSerie))
                     TimeSerie$time <- 1:(nrow(TimeSerie))
                     TimeSerie$week.year <- rep(1:24, length(yy))
                     coefT <- 2 }
                     
                        else { if (time.step=="Semi-fortnight") {
                                TimeSerie$MONTHS <- do.call("rbind", strsplit(TimeSerie$time,"\\.")) [,2]
                                TimeSerie$week.month <- do.call("rbind", strsplit(TimeSerie$time,"\\.")) [,1]
                                TimeSerie$MONTHS <- as.numeric(TimeSerie$MONTHS)
                                TimeSerie$week.month <- as.numeric(TimeSerie$week.month)
                                TimeSerie$time <- NULL
                                names(TimeSerie)[c(1)] <- c("param")
                                names(TimeSerie)[c(2)] <- c("YEARS")
                                Y <- as.numeric(row.names(v))
                                TimeSerie$YEARS <- c(Y)
                                TimeSerie <- TimeSerie[order(TimeSerie$YEARS, TimeSerie$MONTHS), ]
                                TimeSerie <- TimeSerie[TimeSerie$YEARS!=(-1000), ]
                                row.names(TimeSerie) <- 1:(nrow(TimeSerie))
                                TimeSerie$time <- 1:(nrow(TimeSerie))
                                TimeSerie$week.year <- rep(1:48, length(yy)) 
                                coefT <- 4 }

                                   else { if (time.step == "Mono-mensual") {
                                               
                                               TimeSerie$time <- NULL
                                               names(TimeSerie)[c(1)] <- c("param")
                                               names(TimeSerie)[c(2)] <- c("MONTHS")
                                               row.names(TimeSerie) <- 1:(nrow(TimeSerie))
                                               TimeSerie$time <- 1:(nrow(TimeSerie))} 
                                               
                                          else {return(print("No options selected!"))}   }  }  }  }  }  }
#_________________________________________________________________________________________________________________________Selection des mois a traiter et suppression des valeurs extrapolees en debut et fin de serie 
 {
  if (time.step == "Mono-mensual") { F <- 1 }  else {
      F <- length(TimeSerie$param[TimeSerie$YEARS == min(TimeSerie$YEARS)])                                          # trouve la frequence de la serie de base
      
      if (time.step == "Monthly" | time.step == "Semi-fortnight" | time.step == "Fortnight") {
      Fm <- (TSs[order(TSs$YEARS,TSs$MONTHS,TSs$week.month), ]$MONTHS[1])-1                                          # trouve le premier mois de la serie temporelle originale
            if (time.step == "Monthly") { Fw <- 1 } 
                else{ Fw <- TSs[order(TSs$YEARS,TSs$MONTHS,TSs$week.month), ]$week.month[1] }                        # trouve la premiere semaine de la serie temporelle originale
            if (((Fm*coefT+Fw)-1) != 0 ) {                                                                           # si les dates sont differentes
            TimeSerie$param[1:((Fm*coefT+Fw)-1)] <- NA } else {}                                                     #    on remplace le debut de la ts par des NA
      Lm <- 12-((TSs[order(TSs$YEARS,TSs$MONTHS,TSs$week.month), ]$MONTHS)[length(na.omit(TSs$MONTHS))])             # idem avec la fin de la serie temporelle
            if (time.step == "Monthly") { Lw <- 1 } 
                else{ Lw <- (TSs[order(TSs$YEARS,TSs$MONTHS,TSs$week.month), ]$week.month)[length(na.omit(TSs$week.month))] }
            if (((Lm*coefT)+(coefT-Lw)) != 0 ) {
            TimeSerie$param[(((length(TimeSerie$param)+1)-((Lm*coefT)+(coefT-Lw)))):length(TimeSerie$param)] <- NA  } else {}  } else{}  }            
  
  if (any(colnames(TimeSerie) == "MONTHS")) {
      TimeSerie <- subset(TimeSerie, MONTHS %in% months, drop = TRUE) }                                                         # on extrait les mois voulu
  else{}
 }
#__________________________________________________________________________________________________________________Creation de la serie temporelle 'z'
 {          
#_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ choix de la frequence en fonction du time.step et des mois selectionnes
           
  if (time.step=="Mono-mensual") {
      start.year <- min(TimeSerie$MONTHS) }
  else { start.year <- min(TimeSerie$YEARS)  }
  
  start.time <- min(na.omit(TimeSerie)["time"])
   
  freq <- 1/F    
 
  if (!is.null(months) & any(colnames(TimeSerie) == "MONTHS")) {
      freq <- 1/((F/12)*(length(months))) }
  else{}                 
 }
#_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ Construction de la serie temporelle et remplacement des dernieres valeurs manquantes par regression lineaire
 {      
     z <- ts(as.numeric(TimeSerie$param), start = c(start.year), deltat = freq)
     z1 <- z
     z2 <- z
      if (na.replace=="YES") {
         #z <- interpNA(z, "linear")                                                              # interpolation des NAs restant en faisant une regression lineaire entre les valeur avant et apres
         z <- interpTs(z, "linear")
         z <- ts(z, start = c(start.year), deltat = freq)                                        # dans cette version il reste les NAs du debut et de la fin qui ne doivent pas etre extrapoles
         z1 <- z                                                                                 # garde la serie avec NAs au debut et a la fin (z1) pour afficher le tableau (evite les problemes de longueur de colonne)
         z2 <- na.omit(z, method = c("ir"), interp = c("linear"))                                # supprime les NAs en debut et fin de serie temporelle
         z2 <- ts(as.numeric(z2), start = c(start.year, start.time), deltat = freq)              # on reconstruit la serie pour eliminer les probemes de serie temporelle multivariee  (ss NA pour spectrum() et detrend() )
         }  else{} 
   
       if (time.step == "Mono-mensual") { time.xlab <- c("Months") } else { time.xlab <- c("Years")  }
        
       if (plotZ == "YES") {
         dir.create(paste(Envir$save.WD,"/",Envir$File.Name,"/", liste.stations, "/", start,"-", end, "/", strsplit(param,"/")[[1]][1], "/", "na.", na.replace, "-", "out.", outliers.re, "/", time.step, "-", aggreg, "/", sep= ""), recursive = TRUE, showWarnings = FALSE)
         save.figure.path <- paste(Envir$save.WD,"/",Envir$File.Name,"/", liste.stations, "/", start,"-", end, "/", strsplit(param,"/")[[1]][1], "/", "na.", na.replace, "-", "out.", outliers.re, "/", time.step, "-", aggreg, "/", Envir$File.Name, "_TimeSeries_", strsplit(param,"/")[[1]][1],".png", sep = "")
         if (nchar(save.figure.path)>259) { return(tkmessageBox(message= paste("The save path is too long (NTFS system limit)", "\n", "Your results cannot be saved properly", "\n", "Please consider shorter parameter and/or station name", sep=""), icon = "warning", type = "ok", title="!Warning!")) } else { }
        png(save.figure.path, width=1300,height=700, res=100)
         plot(z, ylab="", xlab="", type = "o", pch = 20, col="aquamarine3", panel.first =grid(lwd=0.3, lty="dotted"), axes=FALSE, lwd=1.5, main=paste("Regularised Time Series of",param,"\n"))
         #title(main=paste("\n\n", "Categorical factor(s): ", liste.stations, "\n", "Time step: ", time.step, "   Method of aggregation: ", aggreg), cex.main=0.8)
         box(which = "plot", col="steelblue2", lwd=2)
         if (time.step == "Mono-mensual") { axis(side=1, at = (1:12), labels=month.abb, col="steelblue2", col.ticks="black", tck=-0.01, lwd=1)
                                            title(main=paste("\n\n",  "Station(s):", liste.stations, "    Years:", start, "-", end), cex.main=0.8) }
         else { axis(side=1, at = seq(start(z)[1], end(z)[1]+1, by = 1), col="steelblue2", col.ticks="black", tck= -0.02, lwd=1)
                axis(side=1, at = seq(start(z)[1], end(z)[1]+1, by = 0.5), col="steelblue2", col.ticks="black", labels = FALSE, tck= -0.01, lwd=1)
                title(main=paste("\n\n", "Categorical factor(s):", liste.stations), cex.main=0.8)     }
         axis(side=2, col="steelblue2", col.ticks="black", tck=-0.02, lwd=2)
         axis(side=2, at = seq(par()$yaxp[1], par()$yaxp[2], by = ((par()$yaxp[2]-par()$yaxp[1])/(par()$yaxp[3]*2)))
              , labels=FALSE, col="steelblue2", col.ticks="black", tck=-0.01, lwd=1)
         mtext(paste0(time.xlab), side=1, line=2.5, font=1, cex=1.2, col="black")
         mtext(paste(param), side=2, line=2.5, font=1, cex=1.2, col="black")
        dev.off()

         dev.new(width=8, height=5)
         plot(z, ylab="", xlab="", xaxt="n", type = "o", pch = 20, col="black", panel.first =grid(lwd=0.3, lty="dotted"), axes=FALSE, lwd=1.5, main=paste("Regularised Time Series of",param,"\n"))
         #title(main=paste("\n\n", "Categorical factor(s): ", liste.stations, "\n", "Time step: ", time.step, "   Method of aggregation: ", aggreg), cex.main=0.8)
         box(which = "plot", col="black", lwd=1)
            if (time.step == "Mono-mensual") {  axis(side=1, at =(1:12), col="black", col.ticks="black", labels = month.abb, tck= -0.01, lwd=1)
                                                title(main=paste("\n\n", "Station(s):", liste.stations, "    Years:",  start, "-", end), cex.main=0.8)  }
            else { axis(side=1, at = seq(start(z)[1], end(z)[1]+1, by = 1), col="black", col.ticks="black", tck=-0.02, lwd=1)
                   axis(side=1, at = seq(start(z)[1], end(z)[1]+1, by = 0.5), col="black", col.ticks="black", labels = FALSE, tck= -0.01, lwd=1)
                   title(main=paste("\n\n", "Categorical factor(s):", liste.stations), cex.main=0.8)  }
         axis(side=2, col="steelblue2", col.ticks="black", tck=-0.02, lwd=1)
         axis(side=2, at = seq(par()$yaxp[1], par()$yaxp[2], by = ((par()$yaxp[2]-par()$yaxp[1])/(par()$yaxp[3]*2)))
              , labels=FALSE, col="black", col.ticks="black", tck=-0.01, lwd=1)
         mtext(paste0(time.xlab), side=1, line=2.5, font=1, cex=1.2, col="black")
         mtext(paste(param), side=2, line=2.5, font=1, cex=1.2, col="black")
         }
       else{ } 
 }            
#____________________________________________________________________________________Construction du tableau de donnees regularisees avec NA remplaces
 {     
     Regularised.Data <- TimeSerie                      
     Param <- as.numeric(z)                            
     Regularised.Data$param <- Param                    # tableau de donnees avec NA remplaces
     Regularised.data <- Regularised.Data               # on garde la serie avec le nom inchange pour les tests qui ont besoin de la colone standard 'param'
     Name <- names(Regularised.Data)                    # pour la sortie tableau
     Name[1] <- param 
     colnames(Regularised.Data) <- Name 
 
 dir.create(paste(Envir$save.WD,"/",Envir$File.Name,"/", liste.stations, "/", start,"-", end, "/", strsplit(param,"/")[[1]][1], "/", "na.", na.replace, "-", "out.", outliers.re, "/", time.step, "-", aggreg, "/", sep= ""), recursive = TRUE, showWarnings = FALSE)                              # sauve le tableau
 save.regdata.path <- paste(Envir$save.WD,"/",Envir$File.Name,"/", liste.stations, "/", start,"-", end, "/", strsplit(param,"/")[[1]][1], "/", "na.", na.replace, "-", "out.", outliers.re, "/", time.step, "-", aggreg, "/", Envir$File.Name, "_Regularised_data_", strsplit(param,"/")[[1]][1],".txt", sep = "") 
 if (nchar(save.regdata.path)>259) { return(tkmessageBox(message= paste("The save path is too long (NTFS system limit)", "\n", "Your results cannot be saved properly", "\n", "Please consider shorter parameter and/or station name", sep=""), icon = "warning", type = "ok", title="!Warning!")) } else { }
 write.table(Regularised.Data, sep="\t", row.names=FALSE, file=save.regdata.path)
     
 if (datashow=="YES") { showData(Regularised.Data) } else {}   # affiche le tableau      
 }
#______________________________________________________________________________________________________Statistique desciptive sur la serie regularisee
 {  
  if (resume.reg=="YES") {                                                                   
    summary <- capture.output(summary(Regularised.data$param))
    tkinsert(Envir$txt2, "end", paste("-Desciptive statistics on regularised data-"
                                 , "\n", "Parameter: ", param , "\n", "Categorical factor(s): ", liste.stations, "\n"
                                 , "Outliers.removed: ", outliers.re, "   NA.replaced: ", na.replace, "\n"
                                 , "Time.step: ", time.step, "   Method: ", aggreg, sep="" , "\n\n")) 
    tktag.add(Envir$txt2, "titre", "end -7 lines linestart","end -7 lines lineend")
    tktag.configure(Envir$txt2, "titre", font=tkfont.create(family="courier",size=10,weight="bold"))
    summary <- paste(summary,collapse="\n")
    tkinsert(Envir$txt2,"end",paste(summary, "\n\n"))
    tkinsert(Envir$txt2,"end", paste("Data quantity (",param,") : ", length(Regularised.data$param[!is.na(Regularised.data$param)]), "\n\n", sep=""))
    tksee (Envir$txt2,"end") 
    
    summary <- as.matrix(summary(Regularised.data$param))
    summary <- as.data.frame(summary)
    return(names(summary) <- param)
    }
    
  else{}
 }
#___________________________________________________________________________________Tests de normalite sur la serie temporelle (option: test.normality)
 { 
  if (test.normality == "YES") {
      dev.new()
      qqnorm(z)                                                                                     # Q-Q plot
      qqline(z)
      shap <- shapiro.test(z)                                                                  
      shap.p <- shap$p.value                                                             # extrait la valeur de p du test Shapiro et la comparer au niveau de significativite
      shap.stat <- round(shap$statistic, digits = 4)  
      if (shap.p <= 0.05) {  
          tkinsert(Envir$txt2, "end", paste("-Shapiro-Wilk normality test result-", "\n"))
          tktag.add(Envir$txt2, "titre", "end -2 lines linestart","end -2 lines lineend")
          tktag.configure(Envir$txt2, "titre", font=tkfont.create(family="courier",size=10,weight="bold"))
          tkinsert(Envir$txt2, "end", paste("Parameter: ", param, "\n", "Categorical factor(s): ", liste.stations, "\n"
                                      , "Outliers.removed: ", outliers.re, "   NA.replaced: ", na.replace, "\n"
                                      , "Time.step: ", time.step, "   Method: ", aggreg, sep="" , "\n\n"))                   
          tkinsert(Envir$txt2, "end", paste("Normality of the distribution is rejected","\n"))
          tkinsert(Envir$txt2, "end", paste("W =", shap.stat, "  p-value =",shap.p, "\n\n"))
          tksee (Envir$txt2,"end") 
          
          tkmessageBox (title="Shapiro-Wilk normailty test result",
                  icon = "info" ,
                  message = "Normality of the distribution is rejected",
                  detail = paste("W =", shap.stat, "p-value =",round(shap.p, 4)) ,
                  type = "ok")     }
       else{ 
          tkinsert(Envir$txt2, "end", paste("-Shapiro-Wilk normality test result-", "\n"))
          tktag.add(Envir$txt2, "titre", "end -2 lines linestart","end -2 lines lineend")
          tktag.configure(Envir$txt2, "titre", font=tkfont.create(family="courier",size=10,weight="bold"))
          tkinsert(Envir$txt2, "end", paste("Parameter: ", param, "\n", "Categorical factor(s): ", liste.stations, "\n"
                                      , "Outliers.removed: ", outliers.re, "   NA.replaced: ", na.replace, "\n"
                                      , "Time.step: ", time.step, "   Method: ", aggreg, sep="" , "\n\n"))                  
          tkinsert(Envir$txt2, "end", paste("Normality of the distribution cannot be rejected","\n\n"))
          tkinsert(Envir$txt2, "end", paste("W =", shap.stat, "  p-value =",shap.p, "\n\n"))
          tksee (Envir$txt2,"end") 
          
          tkmessageBox (title="Shapiro-Wilk normailty test result",
                  icon = "info" ,
                  message = "Normality of the distribution cannot be rejected",
                  detail = paste("W =", shap.stat, "p-value =",round(shap.p, 4)) ,
                  type = "ok")   }}
  else{ } 
 }
#_________________________________________________________________________________________________________________________________Affichage du diagramme d'autocorrelation
 { 
  if (autocorr == "YES") { 
          dev.new()
          acf(z, lag.max = ((nrow(TimeSerie))/2), na.action = na.pass, main=paste("Autocorrelogram of", param, "regularised Time Series", "\n"))
          title(main=paste("\n\n", "Categorical factor(s): ", liste.stations, "\n", "Time step: ", time.step, "   Method of aggregation: ", aggreg), cex.main=0.8)
             dir.create(paste(Envir$save.WD,"/",Envir$File.Name,"/", liste.stations, "/", start,"-", end, "/", strsplit(param,"/")[[1]][1], "/", "na.", na.replace, "-", "out.", outliers.re, "/", time.step, "-", aggreg, "/", sep= ""), recursive = TRUE, showWarnings = FALSE)
             save.autocor.path <- paste(Envir$save.WD,"/",Envir$File.Name,"/", liste.stations, "/", start,"-", end, "/", strsplit(param,"/")[[1]][1], "/", "na.", na.replace, "-", "out.", outliers.re, "/", time.step, "-", aggreg, "/", Envir$File.Name, "_AutoCorr_", strsplit(param,"/")[[1]][1],".png", sep = "")
             if (nchar(save.autocor.path)>259) { return(tkmessageBox(message= paste("The save path is too long (NTFS system limit)", "\n", "Your results cannot be saved properly", "\n", "Please consider shorter parameter and/or station name", sep=""), icon = "warning", type = "ok", title="!Warning!")) } else { }
             png(save.autocor.path, width=1000, height=1000, res=150)
             acf(z, lag.max = ((nrow(TimeSerie))/2), na.action = na.pass, main=paste("Autocorrelogram of", param, "regularised Time Series", "\n"))
             title(main=paste("\n\n", "Categorical factor(s): ", liste.stations, "\n", "Time step: ", time.step, "   Method of aggregation: ", aggreg), cex.main=0.8)  
             dev.off()     } 
  else{} 
 }
#________________________________________________________Affichage le spectre de Fourier avec selection manuelle des harmonics maximales (options : autocorr et spectrum)  
 ## !! Les spectres de Fourier ne supporte pas les valeurs manquantes !! ##
 { 
  if (spectrum == "YES") {
    if (length(z2[is.na(z2)]) == 0) {                                                                                  # ne fait pas de spectre si z contient des valeurs manquantes
              Spec <- spectrum(z2, method=c("pgram"), plot=FALSE)
              dev.new()
              plot(1/Spec$freq, Spec$spec, type="l", log="x", lwd=2, main=paste("Spectral density of", param, "regularised Time Series", "\n\n"), ylab="Spectral density", xlab="Frequency = Months / 12")              
              title(main=paste("\n\n", "Categorical factor(s): ", liste.stations, "\n", "Time step: ", time.step, "   Method of aggregation: ", aggreg), cex.main=0.8)
              Pt <-  identify(1/Spec$freq,Spec$spec, plot=TRUE)
              SelecFreq <- round(12/Spec$freq[Pt], digits=0)
              FreqList <- NULL
              for (i in 1:length(SelecFreq)) { FreqList <- paste(FreqList, SelecFreq[i], sep="  ") }
                        tkinsert(Envir$txt2,"end", paste("-Spectral density results-", "\n"))
                        tktag.add(Envir$txt2, "titre", "end -2 lines linestart","end -2 lines lineend")
                        tktag.configure(Envir$txt2, "titre", font=tkfont.create(family="courier",size=10,weight="bold"))
                        tkinsert(Envir$txt2, "end", paste( "Parameter: ", param , "\n", "   Categorical factor(s): ", liste.stations, "\n"
                                                    , "Outliers.removed: ", outliers.re, "   NA.replaced: ", na.replace, "\n"
                                                    , "Time.step: ", time.step, "   Method: ", aggreg, sep="" , "\n\n")) 
                        tkinsert(Envir$txt2,"end", paste("Selected frequencies [months]: ", FreqList, "    months", "\n\n", sep="")) 
                        tksee (Envir$txt2,"end")
                        
#_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ sauvegarde de la figure
                 dir.create(paste(Envir$save.WD,"/",Envir$File.Name,"/", liste.stations, "/", start,"-", end, "/", strsplit(param,"/")[[1]][1], "/", "na.", na.replace, "-", "out.", outliers.re, "/", time.step, "-", aggreg, "/", sep= ""), recursive = TRUE, showWarnings = FALSE)
                 save.spectrum.path <- paste(Envir$save.WD,"/",Envir$File.Name,"/", liste.stations, "/", start,"-", end, "/", strsplit(param,"/")[[1]][1], "/", "na.", na.replace, "-", "out.", outliers.re, "/", time.step, "-", aggreg, "/", Envir$File.Name, "_Spectrum_", strsplit(param,"/")[[1]][1],".png", sep = "")
                 if (nchar(save.spectrum.path)>259) { return(tkmessageBox(message= paste("The save path is too long (NTFS system limit)", "\n", "Your results cannot be saved properly", "\n", "Please consider shorter parameter and/or station name", sep=""), icon = "warning", type = "ok", title="!Warning!")) } else { }
                 png(save.spectrum.path, width=1000, height=700, res=120)
                 plot(1/Spec$freq, Spec$spec, type="l", log="x", lwd=2, main=paste("Spectral density of", param, "regularised Time Series", "\n\n"), ylab="Spectral density", xlab="Frequency = Months / 12")
                 title(main=paste("\n\n", "Categorical factor(s): ", liste.stations, "\n", "Time step: ", time.step, "   Method of aggregation: ", aggreg), cex.main=0.8)  
                 dev.off()    }
         else { return(tkmessageBox(message="Cannot detrend with missing values", icon = "warning", type = "ok", title="!Warning!")) } }         
  else{}
 }
#_________________________________________________________________________________________________________Montre un color plot des anomalies par annee  
 {
  if (anomaly=="YES") {
          if (time.step == "Monthly") {
              cc <- NULL                
              for (i in months) {          
                   l <- Regularised.data$param[Regularised.data$MONTHS==i]                                 # prend toutes les donnees au mois i
                   c <- l-(median(l, na.rm=TRUE))                                                          # soustrait la mediane de l'ensemble de ces donnees a chaque donnee
                   cc <- cbind(cc, c) }                                                                    # superpose les donnees obetnues pour chaque mois
              x <- (min(Regularised.data$YEARS):max(Regularised.data$YEARS))
              y <- (1:length(months))
              dev.new()
              filled.contour(x,y,cc, nlevels= 64, color.palette = colorRampPalette(c("blue","cyan", "green","yellow","orange","red")) , xlab="YEARS", ylab="MONTHS", main = paste("Time series anomaly of", param,"\n\n"), key.title=title(main="Param.\nunit"))
              mtext(paste("Categorical factor(s): ", liste.stations, "\n", "Time step: ", time.step, "   Method of aggregation: ", aggreg, "\n",""), line = 0, cex=0.8)
               dir.create(paste(Envir$save.WD,"/",Envir$File.Name,"/", liste.stations, "/", start,"-", end, "/", strsplit(param,"/")[[1]][1], "/", "na.", na.replace, "-", "out.", outliers.re, "/", time.step, "-", aggreg, "/", sep= ""), recursive = TRUE, showWarnings = FALSE)
               save.colorplot.path <- paste(Envir$save.WD,"/",Envir$File.Name,"/", liste.stations, "/", start,"-", end, "/", strsplit(param,"/")[[1]][1], "/", "na.", na.replace, "-", "out.", outliers.re, "/", time.step, "-", aggreg, "/", Envir$File.Name, "_ColorPlot_", strsplit(param,"/")[[1]][1],".png", sep = "")
               if (nchar(save.colorplot.path)>259) { return(tkmessageBox(message= paste("The save path is too long (NTFS system limit)", "\n", "Your results cannot be saved properly", "\n", "Please consider shorter parameter and/or station name", sep=""), icon = "warning", type = "ok", title="!Warning!")) } else { }
               png(save.colorplot.path, width=1300, height=1000, res=150)
               filled.contour(x,y,cc, nlevels= 64, color.palette = colorRampPalette(c("blue","cyan", "green","yellow","orange","red")) , xlab="YEARS", ylab="MONTHS", main=paste("Time series anomaly of", param, "\n\n"))
               title(main=paste("\n\n", "Categorical factor(s): ", liste.stations, "\n", "Time step: ", time.step, "   Method of aggregation: ", aggreg), cex.main=0.8)
               dev.off() }
              
          else { 
              if (time.step == "Fortnight" | time.step =="Semi-fortnight") {
              cc <- NULL                                      
              for (i in levels(as.factor(Regularised.data$week.year))) {                               
                   l <- Regularised.data$param[Regularised.data$week.year==i]       
                   c <- l-(median(l, na.rm=TRUE))                    
                   cc <- cbind(cc, c) }                            
              x <- (min(Regularised.data$YEARS):max(Regularised.data$YEARS))
              y <- (1:length(levels(as.factor(Regularised.data$week.year))))                                     
              dev.new()
              filled.contour(x,y,cc, nlevels= 64, color.palette = colorRampPalette(c("blue","cyan", "green","yellow","orange","red")) , xlab="YEARS", ylab="WEEKS", main=paste("\n","Time series anomaly of" , param, "\n\n"))
              title(main=paste("\n\n", "Categorical factor(s): ", liste.stations, "\n", "Time step: ", time.step, "   Method of aggregation: ", aggreg), cex.main=0.8) 
              dir.create(paste(Envir$save.WD,"/",Envir$File.Name,"/", liste.stations, "/", start,"-", end, "/", strsplit(param,"/")[[1]][1], "/", "na.", na.replace, "-", "out.", outliers.re, "/", time.step, "-", aggreg, "/", sep= ""), recursive = TRUE, showWarnings = FALSE)
              save.colorplot.path <- paste(Envir$save.WD,"/",Envir$File.Name,"/", liste.stations, "/", start,"-", end, "/", strsplit(param,"/")[[1]][1], "/", "na.", na.replace, "-", "out.", outliers.re, "/", time.step, "-", aggreg, "/", Envir$File.Name, "_Anomaly ColorPlot_", strsplit(param,"/")[[1]][1],".png", sep = "")
              if (nchar(save.colorplot.path)>259) { return(tkmessageBox(message= paste("The save path is too long (NTFS system limit)", "\n", "Your results cannot be saved properly", "\n", "Please consider shorter parameter and/or station name", sep=""), icon = "warning", type = "ok", title="!Warning!")) } else { }
              png(save.colorplot.path, width=1300, height=1000, res=150)
              filled.contour(x,y,cc, nlevels= 64, color.palette = colorRampPalette(c("blue","cyan", "green","yellow","orange","red")) , xlab="YEARS", ylab="WEEKS", main=paste("\n","Time series anomaly of" , param, "\n\n"))
              title(main=paste("\n\n", "Categorical factor(s): ", liste.stations, "\n", "Time step: ", time.step, "   Method of aggregation: ", aggreg), cex.main=0.8)
              dev.off() } 
              else { return(tkmessageBox(message="Cannot perform anomaly color.plot with Climato, yearly or daily time step", icon = "warning", type = "ok", title="!Warning!")) } }  }
  else{}
 } 
#____________________________________________________________________________________________________Histogramme d'anomalie
{
if (a.barplot=="YES") { 

 if (time.step=="Annual") { median.annual <- median(Regularised.data$param, na.rm=TRUE)
                            Median.annual <- rep(median.annual,length(Regularised.data$param)) 
                            hist.anomaly <- Regularised.data$param-Median.annual }
                          
 if (time.step=="Monthly") { m.median <- NULL                           
                             for (i in months) { m.Median <- median(Regularised.data$param[Regularised.data$MONTHS==i], na.rm=TRUE) 
                             m.median <- rbind(m.median,m.Median) }
                             Median.mensual <- rep(m.median, (max(Regularised.data$YEARS)-min(Regularised.data$YEARS))+1)
                             hist.anomaly <- Regularised.data$param-Median.mensual  }

 if (time.step=="Fortnight" | time.step =="Semi-fortnight" ) { m.median <- NULL                           
                              for (i in 1:F) { m.Median <- median(Regularised.data$param[Regularised.data$week.year==i], na.rm=TRUE) 
                              m.median <- rbind(m.median,m.Median)  }
                              Median.fort <- rep(m.median, (max(Regularised.data$YEARS)-min(Regularised.data$YEARS))+1) 
                              hist.anomaly <- Regularised.data$param-Median.fort }
                                                        
 dev.new()
 bp <- barplot(hist.anomaly, ylab= paste(param," anomalies"), xlab="Time", col=c("blue","red")[as.numeric(hist.anomaly>=0)+1], main=paste("Time series anomaly of" , param, "\n\n"))
 title(main=paste("\n\n", "Categorical factor(s): ", liste.stations, "\n", "Time step: ", time.step, "   Method of aggregation: ", aggreg), cex.main=0.8)
 axis(1, at=bp, labels=Regularised.data$YEARS)
 
 dir.create(paste(Envir$save.WD,"/",Envir$File.Name,"/", liste.stations, "/", start,"-", end, "/", strsplit(param,"/")[[1]][1], "/", "na.", na.replace, "-", "out.", outliers.re, "/", time.step, "-", aggreg, "/", sep= ""), recursive = TRUE, showWarnings = FALSE)
            save.barplot.path <- paste(Envir$save.WD,"/",Envir$File.Name,"/", liste.stations, "/", start,"-", end, "/", strsplit(param,"/")[[1]][1], "/", "na.", na.replace, "-", "out.", outliers.re, "/", time.step, "-", aggreg, "/", Envir$File.Name, "_Anomaly BarPlot_", strsplit(param,"/")[[1]][1],".png", sep = "")
            if (nchar(save.barplot.path)>259) { return(tkmessageBox(message= paste("The save path is too long (NTFS system limit)", "\n", "Your results cannot be saved properly", "\n", "Please consider shorter parameter and/or station name", sep=""), icon = "warning", type = "ok", title="!Warning!")) } else { }
            png(save.barplot.path, width=1300, height=700, res=135)
            bp <- barplot(hist.anomaly, ylab= paste(param," anomalies"), xlab="Time", col=c("blue","red")[as.numeric(hist.anomaly>=0)+1], main=paste("Time series anomaly of" , param, "\n\n"))
            title(main=paste("\n\n", "Categorical factor(s): ", liste.stations, "\n", "Time step: ", time.step, "   Method of aggregation: ", aggreg), cex.main=0.8)
            axis(1, at=bp, labels=Regularised.data$YEARS)
            dev.off() 
 }
else{}
} 
#____________________________________________________________________________________________________Montre la courbe detendence ainsi que les residus
 {
  if (zsmooth=="YES")  {
   if (time.step=="Daily" | time.step=="Monthly" | time.step=="Semi-fortnight" | time.step=="Fortnight") {
      if (length(z2[is.na(z2)])==0) {
          z.stl <- stl(z2, s.window="periodic", t.window=(F*10), na.action=na.fail)
      ###construction du tableau a sauvegarder     
          z.date <- as.character(time(z.stl$time.series[, 1])+0.00000001)
          z.Years <- as.numeric(do.call("rbind", strsplit(z.date,"[.]")) [,1])
          z.Time <- round(((as.numeric(paste("0.", do.call("rbind", strsplit(z.date,"[.]")) [,2], sep=""))*F)+1), 0)
          z.tbl <- cbind(z.Years, z.Time, z.stl$time.series[, 1], z.stl$time.series[, 2], z.stl$time.series[, 3])
          colnames(z.tbl) <- c("Years", "Time", "Seasonal", "Trend", "Remainder")
      ###figure principale
          dev.new()
          plot(z.stl)
          title(main=paste("Seasonal Decomposition of Time Series by Loess", "\n\n"), cex.main=1)
      ###figure de la saisonalite
          dev.new()
          z.tbl2 <- as.numeric(as.character(window(z.stl$time.series[, 1], c(min(Regularised.data$YEARS)+1, 1), c(min(Regularised.data$YEARS)+1, F))))
          plot(z.tbl2, ylab=paste(param), xlab=paste(name.freq), main=paste("Seasonality of the Time Series by Loess", "\n\n"), cex.main=1, type="l")
                     
            dir.create(paste(Envir$save.WD,"/",Envir$File.Name,"/", liste.stations, "/", start,"-", end, "/", strsplit(param,"/")[[1]][1], "/", "na.", na.replace, "-", "out.", outliers.re, "/", time.step, "-", aggreg, "/", sep= ""), recursive = TRUE, showWarnings = FALSE)
            save.stltbl.path <- paste(Envir$save.WD,"/",Envir$File.Name,"/", liste.stations, "/", start,"-", end, "/", strsplit(param,"/")[[1]][1], "/", "na.", na.replace, "-", "out.", outliers.re, "/", time.step, "-", aggreg, "/", Envir$File.Name, "_Deseasonal_", strsplit(param,"/")[[1]][1],".txt", sep = "")
            write.table(z.tbl, sep="\t", row.names=FALSE, file=save.stltbl.path)
            save.stlgr.path <- paste(Envir$save.WD,"/",Envir$File.Name,"/", liste.stations, "/", start,"-", end, "/", strsplit(param,"/")[[1]][1], "/", "na.", na.replace, "-", "out.", outliers.re, "/", time.step, "-", aggreg, "/", Envir$File.Name, "_Deseasonal_", strsplit(param,"/")[[1]][1],".png", sep = "")
            save.stlgr2.path <- paste(Envir$save.WD,"/",Envir$File.Name,"/", liste.stations, "/", start,"-", end, "/", strsplit(param,"/")[[1]][1], "/", "na.", na.replace, "-", "out.", outliers.re, "/", time.step, "-", aggreg, "/", Envir$File.Name, "_Seasonality_", strsplit(param,"/")[[1]][1],".png", sep = "")
            if (nchar(save.stlgr.path)>259) { return(tkmessageBox(message= paste("The save path is too long (NTFS system limit)", "\n", "Your results cannot be saved properly", "\n", "Please consider shorter parameter and/or station name", sep=""), icon = "warning", type = "ok", title="!Warning!")) } else { }
            png(save.stlgr.path, width=1200, height=1000, res=160)
            plot(z.stl)
            title(main=paste("Seasonal Decomposition of Time Series by Loess", "\n\n"), cex.main=1)
            dev.off() 
            png(save.stlgr2.path, width=1200, height=1000, res=160)
            plot(z.tbl2, ylab=paste(param), xlab=paste(name.freq), main=paste("Seasonality of the Time Series by Loess", "\n\n"), cex.main=1, type="l")
            dev.off() 
            }
      else { return(tkmessageBox(message="Cannot detrend with missing values", icon = "warning", type = "ok", title="!Warning!")) } }
     else { return(tkmessageBox(message="Cannot detrend with annual or monthly-climato time steps", icon = "warning", type = "ok", title="!Warning!")) } }      
  else { } 
 } 
#________________________________Affichage de la courbe des sommes cumulees pour reperer la tendance generale et les 'cassures' (option : local.trend) 
 {      
  if (local.trend == "YES") {
      if (test != "MK" & test!= "SMK") { return(tkmessageBox(message="Choose a Kendall family test", icon = "warning", type = "ok", title="!Warning!")) } 
       #z <- interpNA(z, "linear")                                                # interpolation des NA restant en faisant une regression lineaire entre les valeur avant et apres
       z <- interpTs(z, "linear")
       z <- ts(as.numeric(z), start = (start.year), deltat = freq)              
      Pna <- na.omit(z)
      dev.new()
      w <- local.trend(Pna)                                                      ## local.trend = package 'pastecs' ##
      pos <- identify(w)
      ww <- pos$pos
      N <- length(ww)-1
      periods <- list()
      lc.mk <- list()
      lc.mk2 <- list()
      dir.create(paste(Envir$save.WD,"/",Envir$File.Name,"/", liste.stations, "/", start,"-", end, "/", strsplit(param,"/")[[1]][1], "/", "na.", na.replace, "-", "out.", outliers.re, "/", time.step, "-", aggreg, "/", sep= ""), recursive = TRUE, showWarnings = FALSE)
    
#_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ Limite de detection des test Kendall (periodes de plus de 1 an minimum)
    LimitD <- ww[2:length(ww)]-ww[1:N]
    if (time.step == "Daily") { if (LimitD<=F) { return(tkmessageBox(message="Selected periods should be longer than 1 year", icon = "warning", type = "ok", title="!Warning!")) } else {} }
    else { if (time.step == "Annual") { if (LimitD<=F) { return(tkmessageBox(message="Selected periods should be longer than 1 year", icon = "warning", type = "ok", title="!Warning!")) } else {} }
           else{ if (time.step == "Monthly") { if (LimitD<=F) { return(tkmessageBox(message="Selected periods should be longer than 1 year", icon = "warning", type = "ok", title="!Warning!")) } else {} } 
                 else{ if (time.step == "Semi-fortnight") { if (LimitD<=F) { return(tkmessageBox(message="Selected periods should be longer than 1 year", icon = "warning", type = "ok", title="!Warning!")) } else {} } 
                       else{ if (time.step == "Fortnight") { if (LimitD<=F) { return(tkmessageBox(message="Selected periods should be longer than 1 year", icon = "warning", type = "ok", title="!Warning!")) } else {} } 
                             else{ if (time.step == "Mono-mensual") { if (LimitD<=F) { return(tkmessageBox(message="Selected periods should be longer than 1 year", icon = "warning", type = "ok", title="!Warning!")) } else {} } }}}}}
      
#_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ Suivi du test de detection de la tendance effectue sur chaque portion de la serie temporelle idenitfiee avec cusum (option : test)
  
    if (time.step == "Annual"){ for (i in 1:N) {
    periods[i] <- list(ts(Regularised.data$param[(ww[i]):(ww[i+1])], start = c(Regularised.data$YEARS[nrow=ww[i]]) , deltat = freq)) }}
    else { if (time.step == "Monthly"){ for (i in 1:N) {
      periods[i] <- list(ts(Regularised.data$param[(ww[i]):(ww[i+1])], start = c((Regularised.data$YEARS[nrow=ww[i]]), (Regularised.data$MONTHS[nrow=ww[i]])) , deltat = freq)) }}
        else { if (time.step == "Semi-fortnight") { for (i in 1:N) {
          periods[i] <- list(ts(Regularised.data$param[(ww[i]):(ww[i+1])], start = c((Regularised.data$YEARS[nrow=ww[i]]),(Regularised.data$week.year[nrow=ww[i]])) , deltat = freq)) }}
            else { if (time.step == "Fortnight") {  for (i in 1:N) {
            periods[i] <- list(ts(Regularised.data$param[(ww[i]):(ww[i+1])], start = c((Regularised.data$YEARS[nrow=ww[i]]),(Regularised.data$week.year[nrow=ww[i]])) , deltat = freq)) }}
                else { if (time.step == "Mono-mensual") { for (i in 1:N) {
                    periods[i] <- list(ts(Regularised.data$param[(ww[i]):(ww[i+1])], start = c(Regularised.data$MONTHS[nrow=ww[i]]) , deltat = freq)) }}
                       else { if (time.step == "Daily") { for (i in 1:N) {
                          periods[i] <- list(ts(Regularised.data$param[(ww[i]):(ww[i+1])], start = c((Regularised.data$YEARS[nrow=ww[i]]), (Regularised.data$DayYears[nrow=ww[i]])) , deltat = freq)) }}}}}}}
 
 
    if (test == "MK") { tkinsert(Envir$txt2,"end", paste("-Local trend results (global)-", "\n"))
                        tktag.add(Envir$txt2, "titre", "end -2 lines linestart","end -2 lines lineend")
                        tktag.configure(Envir$txt2, "titre", font=tkfont.create(family="courier",size=10,weight="bold"))
                        tkinsert(Envir$txt2, "end", paste( "Parameter: ", param, "\n", "Categorical factor(s): ", liste.stations, "\n"
                                                    , "Outliers.removed: ", outliers.re, "   NA.replaced: ", na.replace, "\n"
                                                    , "Time.step: ", time.step, "   Method: ", aggreg, sep="" , "\n\n"))
                        tksee (Envir$txt2,"end")                             
                        for (i in 1:N) { lc.mk[i] <- lapply(periods[[i]], seaKen)
                                         save.localMK.path <- paste(Envir$save.WD,"/",Envir$File.Name,"/", liste.stations, "/", start,"-", end, "/", strsplit(param,"/")[[1]][1], "/", "na.", na.replace, "-", "out.", outliers.re, "/", time.step, "-", aggreg, "/", Envir$File.Name, "_Local_Global Trend_", strsplit(param,"/")[[1]][1],"_"
                                                                    , start(periods[[i]])[1], ",", start(periods[[i]])[2], "-", end(periods[[i]])[1], ",", end(periods[[i]])[2], ".txt", sep = "")
                                         if (nchar(save.localMK.path)>259) { return(tkmessageBox(message= paste("The save path is too long (NTFS system limit)", "\n", "Your results cannot be saved properly", "\n", "Please consider shorter parameter and/or station name", sep=""), icon = "warning", type = "ok", title="!Warning!")) } else { }
                                         write.table(lc.mk[i], sep="\t", row.names=FALSE, file=save.localMK.path)
                                         tkinsert(Envir$txt2,"end", paste("Period ", i, " :  ", start(periods[[i]])[1], ",", start(periods[[i]])[2], " - ", end(periods[[i]])[1], ",", end(periods[[i]])[2]
                                                  , "\n", "Trend (sen.slope): ", round(lc.mk[[i]]$sen.slope, 4), "  original units per year"
                                                  , "\n", "%Trend (sen.slope.pct): ", round(lc.mk[[i]]$sen.slope.pct, 4), "  percent of mean quantity per year", "\n"
                                                  , "p.value: ", round(lc.mk[[i]]$p.value, 4)
                                                  , "\n", "Ref. value of the complete series: ",round(pos$k, 4), "   Mean trend compare to ref. value: ", round(pos$trends[i], 4),"\n\n", sep=""))
                                         tktag.add(Envir$txt2, "titre4", "end -7 lines linestart","end -7 lines lineend")
                                         tktag.configure(Envir$txt2, "titre4", font=tkfont.create(family="courier",size=10)) 
                                         if (lc.mk[[i]]$p.value <= 0.05) { tktag.add(Envir$txt2, "titre3", "end -4 lines linestart","end -4 lines lineend") 
                                                                           tktag.configure(Envir$txt2, "titre3", font=tkfont.create(family="courier",size=9,weight="bold"))} else {} }}
       else { if (test == "SMK") { tkinsert(Envir$txt2,"end", paste("-Local trend results (seasonal)-", "\n"))
                                   tktag.add(Envir$txt2, "titre", "end -2 lines linestart","end -2 lines lineend")
                                   tktag.configure(Envir$txt2, "titre", font=tkfont.create(family="courier",size=10,weight="bold"))
                                   tkinsert(Envir$txt2, "end", paste("Parameter: ", param, "\n", "Categorical factor(s): ", liste.stations, "\n"
                                                               , "Outliers.removed: ", outliers.re, "   NA.replaced: ", na.replace, "\n"
                                                               , "Time.step: ", time.step, "   Method: ", aggreg, "\n\n"
                                                               , "trend (Theil-Sen slope) = original units per year", "\n"
                                                               , "%trend = percent of mean quantity per year", "\n"
                                                               , "missing = proportion of missing slopes", "\n\n", sep=""))
                                   tksee (Envir$txt2,"end")                            
                                   for (i in 1:N) { lc.mk[i] <- lapply(periods[[i]], seasonTrend)
                                                    lc.mk2[i] <- lapply(periods[[i]], seasonTrend, type = c("slope.pct"))
                                                    smk <- cbind(lc.mk[[i]][-5], lc.mk2[[i]][1])
                                                    names(smk)[5]<- "%trend"
                                                    save.localSMK.path <- paste(Envir$save.WD,"/",Envir$File.Name,"/", liste.stations, "/", start,"-", end, "/", strsplit(param,"/")[[1]][1], "/", "na.", na.replace, "-", "out.", outliers.re, "/", time.step, "-", aggreg, "/", Envir$File.Name, "_Local_Seasonal Trend_", strsplit(param,"/")[[1]][1],"_"
                                                                               , start(periods[[i]])[1], ",", start(periods[[i]])[2],"-", end(periods[[i]])[1], ",", end(periods[[i]])[2], ".txt", sep = "")
                                                    if (nchar(save.localSMK.path)>259) { return(tkmessageBox(message= paste("The save path is too long (NTFS system limit)", "\n", "Your results cannot be saved properly", "\n", "Please consider shorter parameter and/or station name", sep=""), icon = "warning", type = "ok", title="!Warning!")) } else { }
                                                    write.table(smk, sep="\t", row.names=FALSE, file=save.localSMK.path) 
                                                    Lc.Mk <- paste(capture.output(smk))
                                                    tkinsert(Envir$txt2,"end", paste("Period ", i, " :  ", start(periods[[i]])[1], ",", start(periods[[i]])[2], " - "
                                                                              , end(periods[[i]])[1], "," , end(periods[[i]])[2], "\n"))
                                                    tkinsert(Envir$txt2,"end", paste(Lc.Mk[1], "\n"))
                                                    for (j in 1:F) { tkinsert(Envir$txt2,"end", paste(Lc.Mk[j+1], "\n")) 
                                                       if (is.na(lc.mk[[i]]$p[j])) {} 
                                                       else { if (lc.mk[[i]]$p[j] <= 0.05) { tktag.add(Envir$txt2, "titre3", "end -2 lines linestart","end -2 lines lineend") 
                                                       tktag.configure(Envir$txt2, "titre3", font=tkfont.create(family="courier", size=10, weight="bold")) } else {}   } }
                                                                                                   
                                                    tkinsert(Envir$txt2,"end", paste("\n", "Ref. value of the complete series: ", round(pos$k, 4), "\n", "Mean trend compare to ref. value: ", round(pos$trends[i], 4), "\n\n", sep="")) }}
           else { return(tkmessageBox(message="Choose a Kendall family test", icon = "warning", type = "ok", title="!Warning!")) }} }
 }                                    
#___________________________________________________________________________Test de detection de la tendance effectue sur la serie temporelle complete
 {
  if (local.trend == "NO") {
   if ( (start-end) == 0 ) { return(tkmessageBox(message="Selected periods should be longer than 1 year", icon = "warning", type = "ok", title="!Warning!")) } else {}
  
#_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ Seasonal Kendall
  if (test == "MK") {
      dir.create(paste(Envir$save.WD,"/",Envir$File.Name,"/", liste.stations, "/", start,"-", end, "/", strsplit(param,"/")[[1]][1], "/", "na.", na.replace, "-", "out.", outliers.re, "/", time.step, "-", aggreg, "/", sep= ""), recursive = TRUE, showWarnings = FALSE)
      save.MK.path <- paste(Envir$save.WD,"/",Envir$File.Name,"/", liste.stations, "/", start,"-", end, "/", strsplit(param,"/")[[1]][1], "/", "na.", na.replace, "-", "out.", outliers.re, "/", time.step, "-", aggreg, "/", Envir$File.Name, "_Global Trend_", strsplit(param,"/")[[1]][1],".txt", sep = "")
      if (nchar(save.MK.path)>259) { return(tkmessageBox(message= paste("The save path is too long (NTFS system limit)", "\n", "Your results cannot be saved properly", "\n", "Please consider shorter parameter and/or station name", sep=""), icon = "warning", type = "ok", title="!Warning!")) } else { }
      write.table(seaKen(z), sep="\t", row.names=FALSE, file=save.MK.path)
      tkinsert(Envir$txt2,"end", paste("-Global trend results-", "\n"))
      tktag.add(Envir$txt2, "titre", "end -2 lines linestart","end -2 lines lineend")
      tktag.configure(Envir$txt2, "titre", font=tkfont.create(family="courier",size=10,weight="bold"))
      tkinsert(Envir$txt2, "end", paste( "Parameter: ", param, "\n", "Categorical factor(s): ", liste.stations, "\n"
                                 , "Outliers.removed: ", outliers.re, "   NA.replaced: ", na.replace, "\n"
                                 , "Time.step: ", time.step, "   Method: ", aggreg, sep="" , "\n\n"))
      mk <- seaKen(z)                                                                                                            ## 'wq' package ##
      tkinsert(Envir$txt2,"end", paste("Trend (sen.slope): ", round(mk$sen.slope, 4), "  original units per year"
                                        , "\n", "%Trend (sen.slope.pct): ", round(mk$sen.slope.pct, 4), "  percent of mean quantity per year", "\n"
                                        , "p.value: ", round(mk$p.value, 4), "\n\n", sep=""))                                      
      if (mk$p.value <= 0.05) { tktag.add(Envir$txt2, "titre3", "end -3 lines linestart","end -3 lines lineend") 
      tktag.configure(Envir$txt2, "titre3", font=tkfont.create(family="courier",size=10,weight="bold"))} else {}
      tksee (Envir$txt2,"end") }
#_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ Seasonal Kendall with details
  else { if (test == "SMK") {
       dir.create(paste(Envir$save.WD,"/",Envir$File.Name,"/", liste.stations, "/", start,"-", end, "/", strsplit(param,"/")[[1]][1], "/", "na.", na.replace, "-", "out.", outliers.re, "/", time.step, "-", aggreg, "/", sep= ""), recursive = TRUE, showWarnings = FALSE)
       save.smk.path <- paste(Envir$save.WD,"/",Envir$File.Name,"/", liste.stations, "/", start,"-", end, "/", strsplit(param,"/")[[1]][1], "/", "na.", na.replace, "-", "out.", outliers.re, "/", time.step, "-", aggreg, "/", Envir$File.Name, "_Seasonal Trend_", strsplit(param,"/")[[1]][1],".txt", sep = "")
       if (nchar(save.smk.path)>259) { return(tkmessageBox(message= paste("The save path is too long (NTFS system limit)", "\n", "Your results cannot be saved properly", "\n", "Please consider shorter parameter and/or station name", sep=""), icon = "warning", type = "ok", title="!Warning!")) } else { }
       tkinsert(Envir$txt2,"end", paste("-Seasonal trend results-", "\n"))
       tktag.add(Envir$txt2, "titre", "end -2 lines linestart","end -2 lines lineend")
       tktag.configure(Envir$txt2, "titre", font=tkfont.create(family="courier",size=10,weight="bold"))
       tkinsert(Envir$txt2, "end", paste( "Parameter: ", param, "\n", "Categorical factor(s): ", liste.stations, "\n"
                                 , "Outliers.removed: ", outliers.re, "   NA.replaced: ", na.replace, "\n"
                                 , "Time.step: ", time.step, "   Method: ", aggreg, sep="" , "\n\n"))
       smk <- seasonTrend(z)
       smk2 <- seasonTrend(z, type=c("slope.pct"))
       SMK <- cbind(smk[-5], smk2[1])
       names(SMK)[5]<- "%trend"      
       write.table(SMK, sep="\t", row.names=FALSE, file=save.smk.path)
       SMK2 <- paste(capture.output(SMK)) 
       tkinsert(Envir$txt2,"end", paste( "trend (Theil-Sen slope) = original units per year", "\n"
                                         , "%trend = percent of mean quantity per year", "\n"
                                         , "missing = proportion of missing slopes", "\n", sep=""))
       tkinsert(Envir$txt2,"end", paste(SMK2[1], "\n"))
       tksee (Envir$txt2,"end")   
       for (i in 1:F) { tkinsert(Envir$txt2,"end", paste(SMK2[i+1], "\n")) 
         if (is.na(smk$p[i])) {}
         else{ if(smk$p[i]<= 0.05) { tktag.add(Envir$txt2, "titre3", "end -2 lines linestart","end -2 lines lineend") 
         tktag.configure(Envir$txt2, "titre3", font=tkfont.create(family="courier", size=10, weight="bold")) } else {}   } }
       tkinsert(Envir$txt2,"end", paste("\n\n", sep=""))  }
 
         else { if (test == "LOWESS") {
             if (time.step=="Daily" | time.step=="Monthly" | time.step=="Semi-fortnight" | time.step=="Fortnight") {
                   dir.create(paste(Envir$save.WD,"/",Envir$File.Name,"/", liste.stations, "/", start,"-", end, "/", strsplit(param,"/")[[1]][1], "/", "na.", na.replace, "-", "out.", outliers.re, "/", time.step, "-", aggreg, "/", sep= ""), recursive = TRUE, showWarnings = FALSE)
                   Loess <- loess(param ~ time, Regularised.data, family="gaussian", span=0.25, control = loess.control(surface = "direct"), na.action=na.exclude)
                   tsLoess <- ts(predict(Loess), start=(min(Regularised.data$YEARS)), deltat=freq)
                   save.loessplot.path <- paste(Envir$save.WD,"/",Envir$File.Name,"/", liste.stations, "/", start,"-", end, "/", strsplit(param,"/")[[1]][1], "/", "na.", na.replace, "-", "out.", outliers.re, "/", time.step, "-", aggreg, "/", Envir$File.Name, "_LOESSplot_", strsplit(param,"/")[[1]][1],".png", sep = "")
                 if (nchar(save.loessplot.path)>259) { return(tkmessageBox(message= paste("The save path is too long (NTFS system limit)", "\n", "Your results cannot be saved properly", "\n", "Please consider shorter parameter and/or station name", sep=""), icon = "warning", type = "ok", title="!Warning!")) } else { }
                  png(save.loessplot.path)
                   plot(z, xlab="Years", ylab=paste(param, "concentration"), main=paste("Trend of", param, "based on LOESS", "\n"))
                   title(main=paste("\n\n", "at station(s): ", liste.stations, "\n", "Time step: ", time.step, "   Method of aggregation: ", aggreg), cex.main=0.8)
                   lines(tsLoess, col="red")
                   legend("topright", inset = 0.01, c("Time Series", "LOESS")
                          , col = c("black", "red"), lwd = 1)
                  dev.off()
                  dev.new()
                       plot(z, xlab="Years", ylab=paste(param, "concentration"), main=paste("Trend of", param, "based on LOESS", "\n"))
                       title(main=paste("\n\n", "at station(s): ", liste.stations, "\n", "Time step: ", time.step, "   Method of aggregation: ", aggreg), cex.main=0.8)      
                       lines(tsLoess, col="red")
                       legend("topright", inset = 0.01, c("Time Series", "LOESS")
                             , col = c("black", "red"), lwd = 1)                                                                                                                                                          
                   mk2 <- seaKen(tsLoess)
                   save.loess.path <- paste(Envir$save.WD,"/",Envir$File.Name,"/", liste.stations, "/", start,"-", end, "/", strsplit(param,"/")[[1]][1], "/", "na.", na.replace, "-", "out.", outliers.re, "/", time.step, "-", aggreg, "/", Envir$File.Name, "_LoessGlobalTrend_", strsplit(param,"/")[[1]][1],".txt", sep = "")
                   if (nchar(save.loess.path)>259) { return(tkmessageBox(message= paste("The save path is too long (NTFS system limit)", "\n", "Your results cannot be saved properly", "\n", "Please consider shorter parameter and/or station name", sep=""), icon = "warning", type = "ok", title="!Warning!")) } else { }
                   write.table(mk2, sep="\t", row.names=FALSE, file=save.loess.path)                                                                                       ## Kendall package ##
                   tkinsert(Envir$txt2,"end", paste("-Global trend results on LOESS-", "\n"))
                   tktag.add(Envir$txt2, "titre", "end -2 lines linestart","end -2 lines lineend")
                   tktag.configure(Envir$txt2, "titre", font=tkfont.create(family="courier",size=10,weight="bold"))
                   tkinsert(Envir$txt2, "end", paste( "Parameter: ", param, "\n", "Categorical factor(s): ", liste.stations, "\n"
                                               , "Outliers.removed: ", outliers.re, "   NA.replaced: ", na.replace, "\n"
                                               , "Time.step: ", time.step, "   Method: ", aggreg, sep="" , "\n\n"))                                                                            
                   tkinsert(Envir$txt2,"end", paste("Trend (sen.slope): ", round(mk2$sen.slope, 4), "  original units per year", "\n"
                                                     , "%Trend (sen.slope.pct): ", round(mk2$sen.slope.pct, 4), "  percent of mean quantity per year", "\n"
                                                     , "p.value: ", round(mk2$p.value, 4), "\n\n", sep="")) 
                   if (mk2$p.value <= 0.05) { tktag.add(Envir$txt2, "titre3", "end -3 lines linestart","end -3 lines lineend") 
                   tktag.configure(Envir$txt2, "titre3", font=tkfont.create(family="courier",size=10,weight="bold"))} else {}
                   tksee (Envir$txt2,"end") } 
                  else { return(tkmessageBox(message="Cannot be perform with annual or monthly-climato time steps", icon = "warning", type = "ok", title="!Warning!")) } }              
                }}      
      }     
      else{ if (test == "NO") { return() } }
 }
#___________________________________________________________________________________________________________________________________fin de la fonction 
}