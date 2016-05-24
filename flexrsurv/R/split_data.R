# Suppression de l'option "relative" du programme de Remontet.
# On ne calcule que la survie relative.

#######################################################
#            UNITE TEMPORELLE : LE MOIS               #
#######################################################


split.data <- function(jeudata, bands, entry, exit, fail, include=names(jeudata), name.runningtime=".t")    
{
  
  
  # Préparation et split des données
  #----------------------------------

  # entry is 1e-12 to prevent further problems

   entry <- ifelse(entry <  1e-12, 1e-12, entry )
   Ljeudata <- Lexis( entry = list(.entry = entry), 
                        exit = list( .entry = exit ),
                  exit.status = fail,
                  data= jeudata[, include],
                  merge=TRUE )
#   without entry (v1.1.2)                  
#   Ljeudata <- Lexis(   exit = list( .entry=exit ),
#                  exit.status = fail,
#                  data=jeudata[, include],
#                  merge=TRUE )
   splitjeudata <- splitLexis( Ljeudata, breaks = bands, time.scale=".entry")

   names(splitjeudata) <- ifelse(names(splitjeudata)=="lex.Xst", ".fail", names(splitjeudata))
#   names(splitjeudata) <- ifelse(names(splitjeudata)=="lex.id",  ".expand", names(splitjeudata))
   names(splitjeudata) <- ifelse(names(splitjeudata)=="lex.dur",  "tik", names(splitjeudata))
   splitjeudata$.exit <- splitjeudata$.entry + splitjeudata$tik

   splitjeudata <- splitjeudata[, c(1:2, dim(splitjeudata)[2], 3:(dim(splitjeudata)[2]-1))]
   

   splitjeudata$tik <- splitjeudata$.exit-splitjeudata$.entry 
  
   splitjeudata$orate <- splitjeudata$rate
   splitjeudata$orate1 <- ifelse(splitjeudata$.fail==0, 0, splitjeudata$rate)
   splitjeudata$newrate<- rep(0, dim(splitjeudata)[1])
   splitjeudata$newrate[cumsum(table(splitjeudata$lex.id))] <- splitjeudata$orate[cumsum(table(splitjeudata$lex.id))]
   splitjeudata$rate <- splitjeudata$newrate
  
    
  #création de la variable interval en numérique !!!
  #-------------------------------------------------
  splitjeudata$intnum<- ifelse(splitjeudata$.fail==1,
                               splitjeudata$.exit,
                               (splitjeudata$.entry+splitjeudata$.exit)/2)

  # on renomme l'intervalle numérique time (utilisé dans la formule)
  names(splitjeudata)<- ifelse(names(splitjeudata)=="intnum", name.runningtime, names(splitjeudata))
  
  return(splitjeudata)  
}
 
