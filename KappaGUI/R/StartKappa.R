StartKappa <-
function(){

GoOn <- 1 # variable globale indiquant si le programme doit continuer ou se terminer

## FENETRE DE CHOIX DE LA PROCEDURE :
proc = Sel_procedure() # on r\'ecup\`ere ici le choix de l'utilisateur

if (is.na(proc$Nb_raters) | (proc$Nb_raters < 2)) {
 GoOn <- 0
 tkmessageBox(message="Invalid number of raters.")
}

## FENETRE DE CHARGEMENT DES DONNEES :
# Faire charger les donn\'ees
if (GoOn == 1) {
 fileName <- tclvalue(tkgetOpenFile(filetypes="{{CSV Files} {.csv}}"))
 if (!nchar(fileName)) {
  tkmessageBox(message = "No file was selected!", icon="error")
  GoOn <- 0  
 }
}

## FENETRE DE MISE EN FORME DES DONNEES :
# L'utilisateur donne des d\'etails sur la mise en forme de ses donn\'ees
if (GoOn == 1) { 
 details_fichier=datafile_details() # on lui demande des d\'etails sur les en-t\^etes de son jeu de donn\'ees
 bdd = read.csv(file=fileName, header=details_fichier[1], sep=";") # dans un premier temps on ne met pas row.names=1, pour v\'erifier que les noms d'individus ne comportent pas de doublons
 if (anyDuplicated(bdd[,1])>0 & (details_fichier[2]==TRUE)) {
  tkmessageBox(title="Duplicated names", message = "Warning : there are duplicates in the names of your individuals. The program will continue normally anyway, but please check your data.", icon = "warning", type = "ok") # on previent l'utilisateur si son fichier a des doublons
 }
 if (details_fichier[2]==TRUE) {bdd = bdd[,-1]} # on retire la premi\`ere colonne (i.e. les noms d'individus) le cas \'ech\'eant
}

## VERIFICATION QUE LE NB DE COLONNES DU FICHIER EST BIEN MULTIPLE DU NB D'OBSERVATEURS
if (GoOn == 1) {
 if (ncol(bdd) %% proc$Nb_raters != 0) {
  tkmessageBox(message="Invalid datafile: the number of columns is not a multiple of the number of raters.")
  GoOn <-0
 }
 #for (j in 1:ncol(bdd)) {
  #bdd[,j] = factor(bdd[,j])
 #}
}

## CALCUL DES KAPPA DE COHEN OU FLEISS :
if (GoOn == 1) {
 if (proc$Nb_raters == 2) { # Kappa Cohen
  MatRes = matrix(ncol =ncol(bdd)/2 , nrow = 2)
  colnames(MatRes) = substr(colnames(bdd)[seq(from=1, to=ncol(bdd), by=proc$Nb_raters)], 1, nchar(colnames(bdd)[seq(from=1, to=ncol(bdd), by=proc$Nb_raters)])-2)
  rownames(MatRes) = c("Values for Cohen's Kappa", "Nb_indiv")
  
  for (j in 1:ncol(MatRes)) {
   z=data.frame(bdd[,2*j], bdd[,2*j-1])
   zz = na.omit(z)
   if (nrow(z)>0) {
    MatRes[1,j] = kappa2(zz, weight=proc$Weighting)$value # ou weight = "unweighted", au choix
    MatRes[2,j] = nrow(zz)
   } else {
    MatRes[1,j] = NA
    MatRes[2,j] = 0
   }
  } 
   
 } else if (proc$Nb_raters > 2) { # Kappa Fleiss
  MatRes = matrix(ncol = ncol(bdd)/proc$Nb_raters, nrow = 1)
  colnames(MatRes) = substr(colnames(bdd)[seq(from=1, to=ncol(bdd), by=proc$Nb_raters)], 1, nchar(colnames(bdd)[seq(from=1, to=ncol(bdd), by=proc$Nb_raters)])-2)
  rownames(MatRes) = c("Value(s) for Fleiss' Kappa", "Nb_indiv")
  
  for (j in 1:ncol(MatRes)) {
   MatRes[1,j] = kappam.fleiss(data.frame(bdd[ , (1:proc$Nb_raters)+(j-1)*proc$Nb_raters]))$value
   MatRes[2,j] = nrow(na.omit(data.frame(bdd[ , (1:proc$Nb_raters)+(j-1)*proc$Nb_raters])))
  } 
 } 
}

# ET ON PROPOSE FINALEMENT A L'UTILISATEUR DE SAUVEGARDER SA MATRICE DE RESULTATS :
if (GoOn == 1) {
 print(MatRes)
 fileName<-tclvalue(tkgetSaveFile(filetypes="{{CSV Files} {.csv}}"))
 if (!nchar(fileName)) {
  tkmessageBox(message="No file was selected!", icon="error")
 } 
 write.csv2(MatRes, fileName)
}

}
