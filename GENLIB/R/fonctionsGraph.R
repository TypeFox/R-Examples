#5- gen.graph									-> garde article

gen.graph <- function(gen, pro = gen.pro(gen), ancestors = gen.founder(gen), indVarAffected = gen.genout(gen)$ind, 
					varAffected = gen.genout(gen)$ind, cex = 1, col = 0, symbolsize = 1, width = 1, packed = F, align = T, ...)
{
	num.max = max(gen.genout(gen)$ind) + 1
	#On divise les ascendances pour qu'elles soient plus petites
	nouv.gen = gen.branching(gen, pro = pro, ancestors = ancestors)
	#On extrait le tableau d'ascendances
	nouv.asc = gen.genout(nouv.gen)
	nouv.asc.temp = nouv.asc
	#EXCEPTION: la fonction pedigree n'accepte pas qu'un des deux parents soit absent, donc, 
	#il faut gerer cette contrainte
	#Il faut extraire tous les individus mere ou pere dont le numero est a 0 et lui assigner un numero bidon
	nouv.asc$status = rep(0, length(nouv.asc$ind))
	father.0 = length(nouv.asc$father[nouv.asc$father == 0 & nouv.asc$mother != 0])
	mother.0 = length(nouv.asc$mother[nouv.asc$mother == 0 & nouv.asc$father != 0])
	nouv.asc$father[nouv.asc$father == 0 & nouv.asc$mother != 0] = c(num.max:(num.max + father.0 - 1))
	nouv.asc$mother[nouv.asc$mother == 0 & nouv.asc$father != 0] = c((num.max + father.0 + 1):((num.max + father.0 + 1) + mother.0 - 1))
	ligneFin.tab = (length(nouv.asc$ind) + 1)
	#MERE
	if(mother.0 != 0) {
		nouv.asc[ligneFin.tab:(ligneFin.tab + mother.0 - 1), 1] = c((num.max + father.0 + 1):((num.max + father.0 + 1) + mother.0 - 1))
		nouv.asc[ligneFin.tab:(ligneFin.tab + mother.0 - 1), 2] = 0
		nouv.asc[ligneFin.tab:(ligneFin.tab + mother.0 - 1), 3] = 0
		nouv.asc[ligneFin.tab:(ligneFin.tab + mother.0 - 1), 4] = "F"
		nouv.asc[ligneFin.tab:(ligneFin.tab + mother.0 - 1), 5] = 1
	}
	#PERE
	ligneFin.tab = (length(nouv.asc$ind) + 1)
	if(father.0 != 0) {
		nouv.asc[ligneFin.tab:(ligneFin.tab + father.0 - 1), 1] = c(num.max:(num.max + father.0 - 1))
		nouv.asc[ligneFin.tab:(ligneFin.tab + father.0 - 1), 2] = 0
		nouv.asc[ligneFin.tab:(ligneFin.tab + father.0 - 1), 3] = 0
		nouv.asc[ligneFin.tab:(ligneFin.tab + father.0 - 1), 4] = "H"
		nouv.asc[ligneFin.tab:(ligneFin.tab + father.0 - 1), 5] = 1
	}
	#On extrait la colonne sex pour la transformer en 0 ou 1 selon le sex
	nouv.asc.sex = as.character(nouv.asc$sex)
	nouv.asc.sex[nouv.asc.sex == "H" | nouv.asc.sex == "1"] = 0
	nouv.asc.sex[nouv.asc.sex == "F" | nouv.asc.sex == "2"] = 1
	#On convertit la colonne en numerique
	nouv.asc.sex = as.integer(nouv.asc.sex)
	#pedigree de S-PLUS 8
	nouv.asc.df <- data.frame(upn = nouv.asc$ind, dadid = nouv.asc$father, momid = nouv.asc$mother, sex = nouv.asc.sex, status = nouv.asc$status)
	#Appel de la fonction pedigree de S-PLUS 8
	nouv.asc.pedigree <- invisible(kinship2::pedigree(id = nouv.asc.df$upn, dadid = nouv.asc.df$dadid, momid = nouv.asc.df$momid, sex = nouv.asc.df$sex,
									status = nouv.asc.df$status))
	asc.init = gen.genout(gen)
	#Affichage des ascendances 
	if(sum(col) == 0) {
		col = rep(1, length(nouv.asc$ind))
		col[nouv.asc$ind %in% gen.pro(gen)] = 2
		col[nouv.asc$ind %in% gen.founder(gen)] = 5
	}
	texte = as.character(c(varAffected[match(nouv.asc.pedigree$id, indVarAffected)]))
	texte[texte == "NA"] = ""
	plot(nouv.asc.pedigree, symbolsize = symbolsize, width = width, col = col, cex = cex, id = texte, packed = packed, align = align, ...)
}

