selvar <-
function(T) {
# T est un tableau binaire.
# La premiere colonne de T est un facteur determinant des groupes.
# Renvoie le descriptif des effectifs, variable par variable, groupe par groupe.

nvar = ncol(T)-1 # le nombre de variables
groupes = T[,1] # la premiere colonne est l'indicateur de groupes
nbgroupes = nlevels(groupes)
nomgroupe = levels(groupes)

sel = matrix(nrow=2*nbgroupes, ncol=nvar)
colnames(sel) = colnames(T)[-1]

for (j in 1:nvar) {
 x = split(T[,j+1], groupes)
 for (i in 1:nbgroupes) {
 sel[i,j] = length(x[[i]][(is.na(x[[i]])==FALSE)])
 }
 for (i in 1:nbgroupes) {
 sel[i+nbgroupes,j] = length(x[[i]][(is.na(x[[i]])==FALSE) & (x[[i]]==1)]) / sel[i,j]
 }
}

noms = rep(NA, 2*nbgroupes)
for (i in 1:nbgroupes) {
 noms[i] = nomgroupe[i]
}
for (i in 1:nbgroupes) {
 noms[i+nbgroupes] = paste("1_", nomgroupe[i], sep="")
}

rownames(sel) = noms
sel = round(sel,3)

# write.csv2(sel[1:nbgroupes, ], "effectifs_MMD.csv")
# write.csv2(sel[(nbgroupes+1):(2*nbgroupes), ], "prop_MMD.csv")
return(sel)
}

Prepa_MMD <-
function(tab, type="raw_data", k=10, all_vars=FALSE, idiosync=FALSE) {
 # il s'agit d'une fonction pr\'eparant un tableau de r\'esum\'e pour l'analyse par MMD d'un tableau binaire "tab" dont la premiere colonne est un indicateur de groupes.
 # type : sont-ce des donn\'ees brutes ou r\'esum\'ees ?
 # Signification de k : nombre minimal d'individus par groupes pour que le caract\`ere soit retenu dans le calcul du MMD.
 # all_vars : bool\'een indiquant si meme les variables avec que des 0 (ou que des 1) doivent etre retenues.
 # idiosync : bool\'een indiquant si meme les variables avec que des 0 (ou que des 1) *sauf sur un individu* doivent etre retenues.
 
if (type == "raw_data") { ## POUR UN TABLEAU DE DONNEES BRUTES :
 # On s'assure d'abord que le tableau est constitu\'e de facteurs :
 for (j in 1:ncol(tab)){
  tab[,j] = factor(tab[,j])
 }
 
 # On garde en m\'emoire le nombre de groupes presents :
 nb_grp = nlevels(tab[,1])
 # Et on transforme une premi\`ere fois les donn\'ees brutes en donn\'ees r\'esum\'ees :
 tac = selvar(tab)
} else if (type == "summarized_data") { # POUR UN TABLEAU DE DONN\'eES R\'eSUM\'eES : 
 nb_grp = nrow(tab)/2
 tac = tab
}

 # On ne retient que les colonnes suffisamment bien renseignees,
 # i.e. celles qui ont plus de k individus par groupe :
 sel = rep(NA,ncol(tac))
 for (j in 1:ncol(tac)) {
  if (all(tac[1:nb_grp,j]>=k)==TRUE) {
   sel[j] = TRUE
  } else {
   sel[j] = FALSE
  }
 }


if (sum(sel)>1) { # s'il y a au moins 2 colonnes selectionnees, on continue

 # On ne retient donc que les colonnes bien renseignees, et la premiere qui est l'indicateur de groupes dans le cas de donn\'ees brutes :
 if (type == "raw_data") { 
  tab = tab[ , c(TRUE, sel)]
  tac = selvar(tab)
 } else {
  tac = tab[ , sel]
 }
 
 # S\'election suppl\'ementaire, le cas \'ech\'eant, pour \'eliminer des caract\`eres trop similaires dans tous les groupes :
 if (all_vars==FALSE){
  niveaux = rep(NA, ncol(tac))
  for (j in 1:ncol(tac)) {
   niveaux[j] = ifelse(all(tac[(nb_grp+1):nrow(tac), j] == 0) | all(tac[(nb_grp+1):nrow(tac), j] == 1), FALSE, TRUE) 
  }
  tac = tac[ , niveaux]
 }
 
 if (idiosync==FALSE){
  avirer = rep(NA, ncol(tac))
  for (j in 1:ncol(tac)) {
   tabprov = round(tac[1:nb_grp, j] * tac[(nb_grp+1):nrow(tac), j],1)
   tabprov = abs(tac[1:nb_grp, j] - tabprov) 
   avirer[j] = ifelse(sum(tabprov)<=1 | sum(tabprov)>=(sum(tac[1:nb_grp, j])-1), FALSE, TRUE) 
  }
  tac = tac[ , avirer]
 }

 rownames(tac) = c(rownames(tac)[1:nb_grp], paste("1_",rownames(tac)[1:nb_grp],sep="")) 
 #write.csv2(tac, "sortie_prep.csv")
 return(tac)

} else { # aucune variable du tableau de donnees ne permet de remplir les criteres fix\'es par l'utilisateur
 return(0)
}
}
