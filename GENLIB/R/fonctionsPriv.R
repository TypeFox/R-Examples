# 3 - data.sum						-> garde article
# 5 - gen.detectionErreur			-> garde article
# 6 - gen.etiquetteGenMinuscule		-> garde article
# 10 - gen.implex3V					-> garde article
# 13 - gen.initImp3V				-> garde article
# 14 - gen.isGen3V					-> garde article
# 19 - gen.validationAsc				-> garde article
# 20 - gen.validationGen				-> garde article
# 21 - gen.validationGLgen			-> garde article
# 23 - GLapplyCG					-> garde article
# 24 - GLapplyF					-> garde article
# 25 - GLapplyGroup					-> garde article
# 26 - GLapplyPhi					-> garde article
# 27 - GLapplyPhi.mat				-> garde article
# 30 - GLCGGroup					-> garde article
# 32 - GLFGroup					-> garde article
# 38 - GLgen						-> garde article
# 39 - GLgroup					-> garde article
# 44 - GLmulti						-> garde article
# 45 - GLOverGroup2					-> garde article
# 46 - GLOverGroup3					-> garde article
# 47 - GLOverGroup4					-> garde article
# 48 - GLOverNumber1				-> garde article
# 49 - GLOverVector2				-> garde article
# 50 - GLPhiGroup					-> garde article
# 53 - GLPriv.completeness3V			-> garde article
# 55 - GLPriv.entropie3V				-> garde article
# 61 - GLPriv.implex3V				-> garde article
# 62 - GLPriv.initcomp3V				-> garde article
# 63 - GLPriv.initDemifon3V			-> garde article
# 64 - GLPriv.initfon3V				-> garde article
# 65 - GLPriv.initImp3V				-> garde article
# 69 - GLPriv.variance3V				-> garde article
# 72 - GLPrivCG					-> garde article
# 73 - GLPrivCGcumul				-> garde article
# 74 - GLPrivCGgroup				-> garde article
# 75 - GLPrivCGmoyen				-> garde article
# 76 - GLPrivCGPLUS					-> garde article
# 77 - GLPrivCGproduit				-> garde article
# 78 - GLPrivCGtotal				-> garde article
# 79 - GLPrivExtCG					-> garde article
# 80 - GLPrivExtF					-> garde article
# 81 - GLPrivExtFSINGLE				-> garde article
# 82 - GLPrivExtPHI					-> garde article
# 83 - GLPrivExtPHISINGLE			-> garde article
# 86 - GLPrivOcc					-> garde article
# 102 - gen.formatage...				-> garde article

data.sum = function(dfData, bLine = T)
{
	if(bLine == T)
		countOut = apply(dfData, 1, sum)
	else countOut = apply(dfData, 2, sum)
	data.frame(countOut)
}

gen.detectionErreur = function(gen, sorted, pro, named, ancestors, nouvtempsmax, individuals, halfSibling, output, genNo, typecomp,print.it, 
	nbgenerations, depthmin, depthmax, matricephi, prob, b, icmatricephi, inter, correct, correctinterval, label, grppro, vectF,
	typeCG, nogrp, etiquettep, symbole, info, maxindpage, fond, grapheg, cex, font, ..., check = 0)
{
	#Les test #2 l'emporte sur test # 1 si les deux sont present
	#Les test (#1 et #2) et (#3) sont mutuellement exclusive.
	if(( is.element(1, check) && is.element(3, check)) || (is.element(2, check) && is.element(3, check)))
	 	return(list(erreur = T, messageErreur = "Invalid 'check' parameter: test #3 can not be concurent with tests #1 and #2."))
	 	#"Parametre 'check' invalide : Le test #3 ne peut pas etre fait en meme temps que les tests #1 et #2"))
	#validation des parametres gerer dans ... presentement sont : father, mother et sex seulement
	par... = names(list(...))
	if(length(par...) > 0) {
		for(i in 1:length(par...))
			if(!is.element(par...[i], c("mother", "father", "sex")))
				return(list(erreur = T, messageErreur = "Invalid '...' parameter: unknown parameter name"))
				#"Parametre '...' invalide : nom d'un parametre inconnu"))
	}
	#Contient les objets en retour
	retour = list()
	#Parametre 'gen'
	if(is.element(1, check) || is.element(2, check)) {
		ret1 = gen.validationGen(gen = gen, ..., check = check)
	#return(1)
		if(ret1$erreur == T)
			return(ret1)
		else {
			retour$erreur = ret1$erreur
			retour$gen = ret1$gen
			gen = ret1$gen
		}
	}
	else if(is.element(3, check)) {
		ret3 = gen.validationGLgen(gen)
		if(ret3$erreur == T)
			return(ret3)
		else
		{
			retour$erreur = ret3$erreur
			retour$gen = ret3$gen
		}
	}
	#Parametre 'sorted'
	if(is.element(4, check)) {
		if(!is(sorted, "logical"))
			return(list(erreur = T, messageErreur = "Invalid 'sorted' parameter: must be a logical value"))
			#"Parametre 'sorted' invalide: doit etre une valeur logique"))
		retour$sorted = sorted
	}
	#Parametre 'pro'
	if(is.element(5, check)) {
		if(is(pro, "GLgroup"))
			retour$pro = pro
		else {
			if(sum(as.numeric(pro)) == 0)
				pro = gen.pro(gen)
			if(!is(pro, "numeric"))
				return(list(erreur = T, messageErreur = "Invalid 'pro' parameter: must be a numerical vector"))
				#"Parametre 'pro' invalide: doit etre un vecteur numerique"))
			if(is(gen, "GLgen")) {
				#print("check 5 -> gen.genout(gen, check = 0)")
				if(!is.na(match(NA, match(pro, gen.genout(gen)$ind)))) #, check = 0
					return(list(erreur = T, messageErreur = "Invalid 'pro' parameter: one of the proband is not part of the individuals list"))
						#"Parametre 'prop' invalide: L'un des proposants ne fait pas parti de la liste des individuals"))
			}
			else {
				if(!is.na(match(NA, match(pro, gen$ind))))
					return(list(erreur = T, messageErreur = "Invalid 'pro' parameter: one of the proband is not part of the individuals list"))
						#"Parametre 'prop' invalide: L'un des proposants ne fait pas parti de la liste des individuals"))
			}
			retour$pro = pro
		}
	}
	#Parametre 'pro' 
	if(is.element(6, check)) {
		#Plus facile a gerer en matrice
		vectF <- temp <- as.matrix(vectF)
		if(sum(as.numeric(pro)) == 0)
			pro <- as.integer(dimnames(vectF <- temp)[[1]])
		if(!is.numeric(pro) || length(pro) != dim(vectF <- temp)[1])
			return(list(erreur = T, messageErreur = "Invalid 'pro' parameter: must be a vector with the same probands as 'VectF'"))
				#"Parametre 'prop' invalide: doit etre un vecteur contenant les memes proposants que 'VectF'"))
		retour$pro = pro
	}
	#Parametre 'pro' 
	if(is.element(7, check)) {
		if(sum(as.numeric(pro)) == 0)
			pro = gen.pro(gen)
		if(!is(pro, "numeric"))
			return(list(erreur = T, messageErreur = "Invalid 'pro' parameter: must be a numerical vector"))
			#"Parametre 'prop' invalide: doit etre un vecteur numerique"))
		if(is(gen, "GLgen")) {
			if(!is.na(match(NA, match(pro, gen.genout(gen)$ind)))) #, check = 0
				return(list(erreur = T, messageErreur = "Invalid 'pro' parameter: one of the proband is not part of the individuals list"))
					#"Parametre 'prop' invalide: L'un des proposants ne fait pas parti de la liste des individuals"))
		}
		else {
			if(!is.na(match(NA, match(pro, gen$ind))))
				return(list(erreur = T, messageErreur = "Invalid 'pro' parameter: one of the proband is not part of the individuals list"))
					#"Parametre 'prop' invalide: L'un des proposants ne fait pas parti de la liste des individuals"))
		}
		retour$pro = pro
	}
	#Parametre 'pro' 
	if(is.element(8, check)) {
		if(sum(as.numeric(pro)) == 0) {
			if(is.null(dimnames(matricephi)))
				return(list(erreur = T, messageErreur = "Invalid 'pro' parameter: must be a vector with the same probands as 'matricephi'"))
					#"Parametre 'prop' invalide: doit etre un vecteur contenant les memes proposants que 'matricephi'"))
			pro <- as.integer(dimnames(matricephi)[[1]])
		}
		if(!is(pro, "numeric"))
			return(list(erreur = T, messageErreur = "Invalid 'pro' parameter: must be a numerical vector"))
			#"Parametre 'prop' invalide: doit etre un vecteur numerique"))
		retour$pro = pro
	}
	#Parametre 'grppro'
	if(is.element(9, check)) {
		if(is(grppro, "GLgroup"))
			retour$grppro = grppro
		else return(list(erreur = T, messageErreur = "Invalid 'grppro' parameter: must be a 'GLgroup' object"))
		#"Parametre 'grpprop' invalide: doit etre un objet 'GLgroup'"))
	}
	#Parametre 'named'
	if(is.element(10, check)) {
		if(!is(named, "logical"))
			return(list(erreur = T, messageErreur = "Invalid 'named' parameter: must be a logical value"))
			#"Parametre 'named' invalide: doit etre une valeur logique"))
		retour$named = named
	}
	#Parametre 'ancestors'
	if(is.element(11, check)) {
		if(sum(as.numeric(ancestors)) == 0) {
			if(is(gen, "GLgen")) {
				#print("check 11 genout(gen)")
				genTemp = gen.genout(gen)
				#ancestors = genTemp$ind
				ancestors = genTemp[genTemp[,"father"]==0 & genTemp[,"mother"]==0,"ind"]
			}
			else #ancestors = gen$ind
				ancestors = gen[gen[,"father"]==0 & gen[,"mother"]==0,"ind"]
		}
		if(!is(ancestors, "numeric"))
			return(list(erreur = T, messageErreur = "Invalid 'ancestors' parameter: must be a numerical vector"))
				#"Parametre 'ancestors' invalide: doit etre un vecteur numerique"))
		retour$ancestors = ancestors
	}
	#Parametre 'nouvtempsmax' 
	if(is.element(12, check)) {
		newTemps = 0
		if(is.na(nouvtempsmax))
			newTemps <- as.double(-1.)
		else {
			if(is(nouvtempsmax, "character") && nouvtempsmax == "infini")
				newTemps <- as.double(0.)
			else {
				if(!is(nouvtempsmax, "numeric"))
					return(list(erreur = T, messageErreur = "Invalid 'nouvtempsmax' parameter: must be a numerical vector"))
						#"Parametre 'nouvtempsmax' invalide: doit etre un vecteur numerique"))
				newTemps <- as.double(nouvtempsmax)[1]
			}
		}
		retour$newTemps = newTemps
	}
	#Parametre 'individuals'
	if(is.element(13, check)) {
		if(!is(individuals, "numeric"))
			return(list(erreur = T, messageErreur = "Invalid 'individuals' parameter: must be a numerical vector"))
				#"Parametre 'individuals' invalide: doit etre un vecteur numerique"))
		if(is(gen, "GLgen")) {
			genTmp = gen.genout(gen)#, check = 0)
			posTrouve = match(individuals, unique(c(genTmp$ind, genTmp$mother, genTmp$father)))
		}
		else posTrouve = match(individuals, unique(c(gen$ind, gen$mother, gen$father)))
		if(length(posTrouve[is.na(posTrouve)]) > 0)
			return(list(erreur = T, 
					messageErreur = "Invalid 'individuals' parameter: all individuals must be present in the ascendance table"))
				#"Parametre 'individuals' invalide: tous les individuals doivent etre presents dans la table d'ascendance"))
		retour$individuals = individuals
	}
	#Parametre 'halfSibling'
	if(is.element(14, check)) {
		if(!is(halfSibling, "logical"))
			return(list(erreur = T, messageErreur = "Invalid 'halfSibling' parameter: must be a logical value"))
				#"Parametre 'halfSibling' invalide: doit etre une valeur logique"))
		retour$halfSibling = halfSibling
	}
	#Parametre 'output' 
	if(is.element(15, check)) {
		if(!(output == "Fa" || output == "Mo" || output == "FaMo"))
			return(list(erreur = T, messageErreur = "Invalid 'output' parameter: choices are 'Fa', 'Mo', or 'FaMo' (see documentation)"))
				#"Parametre 'sortie' invalide: les choix disponibles sont 'P' , 'M' et 'PM' (voir la documentation)"))
		retour$output = output
	}
	#Parametre 'genNo'
	if(is.element(16, check)) {
		if(!is(genNo, "numeric"))
			return(list(erreur = T, messageErreur = "Invalid 'genNo' parameter: must be a numerical vector"))
			#"Parametre 'genNo' invalide: doit etre un vecteur numerique"))
		depthMax = gen.depth(gen.genealogy(gen))
		if(length(genNo) == 1 && genNo == -1)
			genNo = 0:(depthMax - 1)
		if(max(genNo) > (depthMax - 1))
			return(list(erreur = T, messageErreur = "Invalid 'genNo' parameter: can not be deeper than maximal depth"))
				#"Parametre 'genNo' invalide: ne doit pas depasser la depth maximale"))
		retour$genNo = genNo
	}
	#Parametre 'typeComp'
	if(is.element(171, check)) {
		#print(typecomp)
		if(!(typecomp == "MEAN" || typecomp == "IND")) # || typecomp == "MOYSUJETS")) #typecomp == "CUM" || typecomp == "REL" || typecomp == "EGO" || 
			return(list(erreur = T, 
				messageErreur = "Invalid 'typecomp' parameter: choices are 'IND' and 'MEAN' (see documentation)"))#, 'CUM', 'REL' et 'EGO' "))
				#"Parametre 'typecomp' invalide: les choix disponibles sont 'MOYSUJETS','SUJETS','BRUT', 'CUM', 'REL' et 'EGO' (voir la documentation)"))
		retour$typecomp = typecomp
	}
	if(is.element(17, check)) {
		#print(typecomp)
		if(!(typecomp == "ALL" || typecomp == "IND" || typecomp == "MEAN")) #typecomp == "CUM" || typecomp == "REL" || typecomp == "EGO" || 
			return(list(erreur = T, 
				messageErreur = "Invalid 'typecomp' parameter: choices are 'MEAN','IND','ALL' (see documentation)"))#, 'CUM', 'REL' et 'EGO'"))
				#"Parametre 'typecomp' invalide: les choix disponibles sont 'MOYSUJETS','SUJETS','BRUT', 'CUM', 'REL' et 'EGO' (voir la documentation)"))
		retour$typecomp = typecomp
	}
	#Parametre 'print.it'
	if(is.element(18, check)) {
		if(!is(print.it, "logical"))
			return(list(erreur = T, messageErreur = "Invalid 'print.it' parameter: must be a logical value"))
			#"Parametre 'print.it' invalide: doit etre une valeur logique"))
		retour$print.it = print.it
	}
	#Parametre 'nbgenerations'
	if(is.element(19, check)) {
		if(!is(nbgenerations, "numeric"))
			return(list(erreur = T, messageErreur = "Invalid 'nbgenerations' parameter: must be a numerical vector"))
				#"Parametre 'nbgenerations' invalide: doit etre un vecteur numerique"))
		retour$nbgenerations = nbgenerations
	}
	#Parametre depthmin et depthmax (un ne peux pas aller sans l'autre)
	if(is.element(20, check)) {
		if(!is(depthmin, "numeric"))
			return(list(erreur = T, messageErreur = "Invalid 'depthmin' parameter: must be a numerical vector"))
				#"Parametre 'depthmin' invalide: doit etre un vecteur numerique"))
		if(!is(depthmax, "numeric"))
			return(list(erreur = T, messageErreur = "Invalid 'depthmax' parameter: must be a numerical vector"))
				#"Parametre 'depthmax' invalide: doit etre un vecteur numerique"))
		if(as.integer(depthmax) < as.integer(depthmin))
			return(list(erreur = T, messageErreur = "'depthmax' must be bigger than 'depthmin'"))
			#"'depthmax' doit etre plus grand que 'depthmin'"))
		retour$depthmin = depthmin
		retour$depthmax = depthmax
	}
	#Parametre prob
	if(is.element(21, check)) {
		if(!is(prob, "numeric"))
			return(list(erreur = T, messageErreur = "Invalid 'prob' parameter: must be a numerical vector"))
				#"Parametre 'prob' invalide:  doit etre un vecteur numerique"))
		retour$prob = prob
	}
	#Parametre b
	if(is.element(22, check)) {
		if(!is(b, "numeric") || any(b <= 0))
			return(list(erreur = T, messageErreur = "Invalid 'b' parameter: must be an integer"))
			#"Parametre 'b' invalide:  doit etre un nombre entier positif"))
		retour$b = b
	}
	#Parametre 'icmatricephi'
	if(is.element(23, check)) {
		if((!is(icmatricephi, "GLmultiVector")) & (!is.null(icmatricephi)) & (!is(icmatricephi, "named")))
			return(list(erreur = T, messageErreur = "Invalid 'icmatricephi' parameter: must be a 'GLmultiVector' object"))
				#"Parametre 'icmatricephi' invalide:  doit etre un objet 'GLmultiVector' "))
		retour$icmatricephi = icmatricephi
	}
	#Parametre 'inter'
	if(is.element(24, check)) {
		if(!is(inter, "logical"))
			return(list(erreur = T, messageErreur = "Invalid 'inter' parameter: must be a logical value"))
			#"Parametre 'inter' invalide: doit etre une valeur logique"))
		retour$inter = inter
	}
	#Parametre 'correct'
	if(is.element(25, check)) {
		if(!is(correct, "logical"))
			return(list(erreur = T, messageErreur = "Invalid 'correct' parameter: must be a logical value"))
			#"Parametre 'correct' invalide: doit etre une valeur logique"))
		retour$correct = correct
	}
	#Parametre 'correctinterval'
	if(is.element(26, check)) {
		if(!is(correctinterval, "numeric"))
			return(list(erreur = T, messageErreur = "Invalid 'correctinterval' parameter: must be a numerical vector"))
				#"Parametre 'correctinterval' invalide:  doit etre un vecteur numerique"))
		retour$correctinterval = correctinterval
	}
	#Parametre 'label'
	if(is.element(27, check)) {
		if(is.null(group(matricephi))) {
			if(is.null(label))
				label <- dimnames(matricephi)[[1]]
			if(length(label) != dim(matricephi)[1])
				return(list(erreur = T, messageErreur = "Invalid 'label' parameter: label vector length must equal proband number in 'matricephi'"))
					#"Parametre 'label' invalide: La longueur du vecteur 'label' doit etre egale au nombre de proposants dans 'matricephi'"))
		}
		else {
			if(is.null(label))
				label <- names(group(matricephi))
			if(length(label) != length(matricephi@group))
				return(list(erreur = T, messageErreur = "Invalid 'label' parameter: label number must equal group number in 'matricephi'"))
					#"Parametre 'label' invalide: Le nombre d'etiquette doit etre egale au nombre de group dans 'matricephi'"))
		}
		retour$label = label
	}
	#Parametre 'matricephi' (4 cas sur 4)(
	if(is.element(28, check)) {
		if(is(matricephi, "GLmultiPhiGroupSingle")) {
			retour$matricephi = matricephi
		}
		else if(is(matricephi, "GLmultiPhiGroup")) {
			retour$matricephi = matricephi
		}
		else if(is(matricephi, "GLmultiMatrix")) {
			retour$matricephi = matricephi
		}
		else if(is(matricephi, "matrix")) {
			if(!is.array(matricephi))
				return(list(erreur = T, messageErreur = "Invalid 'matricephi' parameter: must be a kinship matrix for only one depth"))
					#"Parametre 'matricephi' invalide: doit etre une matrice d'kinship pour une seule depth"))
			if(!is.numeric(matricephi) || is.null(dim(matricephi)) || dim(matricephi)[1] != dim(matricephi)[2])
				return(list(erreur = T, messageErreur = "Invalid 'matricephi' parameter: line number must equal column number"))
					#"Parametre 'matricephi' invalide: doit contenir autant de lignes que de colonnes (matrice carree)"))
			retour$matricephi = matricephi
		}
		else {
			return(list(erreur = T, messageErreur = "Invalid 'matricephi' parameter: must be one of 'matrix','GLmultimatrix', 'GLmultiPhiGroup' or 'GLmultiPhiGroupSingle' valid object"))
				#"Parametre 'matricephi' invalide: doit etre un objet 'matrix','GLmultimatrix', 'GLmultiPhiGroup' ou 'GLmultiPhiGroupSingle' valide"))
		}
	}
	#Parametre 'matricephi' (2 cas sur 4)(
	if(is.element(29, check)) {
		if(is(matricephi, "GLmultiMatrix")) {
			retour$matricephi = matricephi
		}
		else if(is(matricephi, "matrix")) {
			if(!is.array(matricephi))
				return(list(erreur = T, messageErreur = "Invalid 'matricephi' parameter: must be a kinship matrix for only one depth"))
					#"Parametre 'matricephi' invalide: doit etre une matrice d'kinship pour une seule depth"))
			if(!is.numeric(matricephi) || is.null(dim(matricephi)) || dim(matricephi)[1] != dim(matricephi)[2])
				return(list(erreur = T, messageErreur = "Invalid 'matricephi' parameter: line number must equal column number"))
					#"Parametre 'matricephi' invalide: doit contenir autant de lignes que de colonnes (matrice carree)"))
			retour$matricephi = matricephi
		}
		else {
			return(list(erreur = T, messageErreur = "Invalid 'matricephi' parameter: must be one of 'matrix' or 'GLmultimatrix' valid object"))
				#"Parametre 'matricephi' invalide: doit etre un objet 'matrix' ou 'GLmultimatrix' valide"))
		}
	}
	#Parametre 'matricephi' (2 cas sur 2)
	if(is.element(30, check)) {
		if(is(matricephi, "GLmultiPhiGroupSingle")) {
			retour$matricephi = matricephi
		}
		else if(is(matricephi, "matrix")) {
			if(!is.array(matricephi))
				return(list(erreur = T, messageErreur = "Invalid 'matricephi' parameter: must be a kinship matrix for only one depth"))
					#"Parametre 'matricephi' invalide: doit etre une matrice d'kinship pour une seule depth"))
			if(!is.numeric(matricephi) || is.null(dim(matricephi)) || dim(matricephi)[1] != dim(matricephi)[2])
				return(list(erreur = T, messageErreur = "Invalid 'matricephi' parameter: line number must equal column number"))
					#"Parametre 'matricephi' invalide: doit contenir autant de lignes que de colonnes (matrice carree)"))
			retour$matricephi = matricephi
		}
		else return(list(erreur = T, messageErreur = "Invalid 'matricephi' parameter: must be one of 'matrix' or 'GLmultiPhiGroupSingle' valid object"))
				#"Parametre 'matricephi' invalide: doit etre un objet 'matrix' ou 'GLmultiPhiGroupSingle' valide"))
	}
	#Parametre 'matricephi' (2 cas sur 2)
	if(is.element(31, check)) {
		if(is(matricephi, "GLmultiMatrix"))
			retour$matricephi = matricephi
		else if(is(matricephi, "GLmultiPhiGroup"))
			retour$matricephi = matricephi
		else return(list(erreur = T, messageErreur = "Invalid 'matricephi' parameter: must be one of 'GLmultiMatrix' or 'GLmultiPhiGroup' valid object"))
				#"Parametre 'matricephi' invalide: doit etre un objet 'GLmultiMatrix' ou 'GLmultiPhiGroup' valide"))
	}
	#Parametre 'vectF'  
	if(is.element(32, check)) {
		if(!is.numeric(vectF))
			return(list(erreur = T, messageErreur = "Invalid 'vectF' parameter: must be a numeric vector or a 'GLmultiVector'"))
				#"Parametre 'vectF' invalide: doit etre un vecteur numerique ou un 'GLmultiVector'"))
		if(sum(as.numeric(pro)) == 0)
			if(!(is(vectF, "named") || is(vectF, "GLmultiVector")))
				return(list(erreur = T,
					messageErreur = "Invalid 'vectF' parameter: must be a labeled numeric vector or the 'pro' parameter becomes obligatory"))
					#"Parametre 'vectF' invalide: doit etre un vecteur numerique etiquette ou le parametre 'prop' devient obligatoire"))
		retour$vectF = vectF
	}
	#Parametre 'vectF'  
	if(is.element(33, check)) {
		if(is.numeric(vectF) || is(vectF, "GLmultiVector") || is(vectF, "GLmultiFGroupSingle") || is(vectF, "GLmultiFGroup")
			)
			retour$vectF = vectF
		else return(list(erreur = T, 
			messageErreur = "Invalid 'vectF' parameter: must be a numeric vector or one of 'GLmultiVector', 'GLmultiFGroupSingle' or 'GLmultiFGroup'"))
			#"Parametre 'vectF' invalide: doit etre un vecteur numerique, un 'GLmultiVector', un 'GLmultiFGroupSingle' ou un 'GLmultiFGroup'"))
	}
	#Parametre 'typeCG'
	if(is.element(34, check)) {
		if(!(typeCG == "IND" || typeCG == "MEAN" || typeCG == "CUMUL" || typeCG == "TOTAL" || typeCG == "PRODUCT"))
			return(list(erreur = T, 
				messageErreur = "Invalid 'typeCG' parameter: choices are 'IND', 'MEAN', 'CUMUL' , 'TOTAL' or 'PRODUCT' (see documentation)"))
				#"Parametre 'typeCG' invalide: les choix disponibles sont 'BRUT', 'MOYEN', 'CUMUL' , 'TOTAL' et 'PRODUIT' (voir la documentation)"))
		retour$typeCG = typeCG
	}
	#Parametre 'nogrp'
	if(is.element(35, check)) {
		if(!is(nogrp, "numeric"))
			return(list(erreur = T, messageErreur = "Invalid 'nogrp' parameter: must be a numerical vector"))
				#"Parametre 'nogrp' invalide:  doit etre un vecteur numerique"))
		retour$nogrp = nogrp
	}
	#Parametre 'pro'
	if(is.element(36, check)) {
		if(!is(pro, "numeric"))
			return(list(erreur = T, messageErreur = "Invalid 'pro' parameter: must be a numerical vector"))
				#"Parametre 'prop' invalide: doit etre un vecteur numerique"))
		retour$pro = pro
	}
	#Parametre 'ancestors'
	if(is.element(37, check)) {
		if(!is(ancestors, "numeric"))
			return(list(erreur = T, messageErreur = "Invalid 'ancestors' parameter: must be a numerical vector"))
				#"Parametre 'ancestors' invalide: doit etre un vecteur numerique"))
		retour$ancestors = ancestors
	}
	#Parametre 'etiquettep' a finir
	if(is.element(38, check)) {
		retour$etiquettep = etiquettep
	}
	#Parametre 'symbole'
	if(is.element(39, check)) {
		if(!is(symbole, "numeric"))
			return(list(erreur = T, messageErreur = "Invalid 'symbole' parameter: must be a numerical vector"))
				#"Parametre 'symbole' invalide: doit etre un vecteur numerique"))
		retour$symbole = symbole
	}
	#Parametre 'info'
	if(is.element(40, check)) {
		retour$info = info
	}
	#Parametre 'maxindpage'
	if(is.element(41, check)) {
		if(!(is(maxindpage, "numeric") && length(maxindpage) == 1))
			return(list(erreur = T, messageErreur = "Invalid 'maxindpage' parameter: must be a numerical vector of size 1"))
				#"Parametre 'maxindpage' invalide: doit etre un vecteur numerique de longueur 1"))
		retour$maxindpage = maxindpage
	}
	#Parametre 'fond'
	if(is.element(42, check)) {
		if(!is(fond, "logical"))
			return(list(erreur = T, messageErreur = "Invalid 'fond' parameter: must be a logical value"))
				#"Parametre 'fond' invalide: doit etre une valeur logique"))
		retour$fond = fond
	}
	#Parametre 'grapheg'
	if(is.element(43, check)) {
		if(!is(grapheg, "logical"))
			return(list(erreur = T, messageErreur = "Invalid 'grapheg' parameter: must be a logical value"))
			 #"Parametre 'grapheg' invalide: doit etre une valeur logique"))
		retour$grapheg = grapheg
	}
	#Parametre 'cex'
	if(is.element(44, check)) {
		if(!(is(cex, "numeric") && length(cex) == 1))
			return(list(erreur = T, messageErreur = "Invalid 'cex' parameter: must be a numerical vector of size 1"))
				#"Parametre 'cex' invalide: doit etre un vecteur numerique de longueur 1"))
		retour$cex = cex
	}
	#Parametre 'font'
	if(is.element(45, check)) {
		if(!(is(font, "numeric") && length(font) == 1))
			return(list(erreur = T, messageErreur = "Invalid 'font' parameter: must be a numerical vector of size 1"))
				#"Parametre 'font' invalide: doit etre un vecteur numerique de longueur 1"))
		retour$font = font
	}
	#valeur de retour
	# S'il y a une erreur
	#  1- erreur
	#  2- message d'erreur
	# Si tout est valide
	#  1- erreur
	#  2- dataframe valide
	# ...
	#Si cette fonction est rendu la, il n'y  a pas d'erreur
	retour$erreur = F
	return(retour)
}

gen.etiquetteGenMinuscule = function(gen)
{
	retour = list()
	if(length(gen$IND) > 0 || length(gen$ind) > 0) {
		if(length(gen$IND) > 0)
			names(gen)[names(gen) == "IND"] = "ind"
	}
	else {
		retour$erreur = T
		retour$messageErreur = "Invalid 'gen' parameter: ascendance table must contain one column named 'ind' or 'IND'"
			#"Parametre 'gen' invalide : la table d'ascendance doit contenir une colonne se nommant 'ind' ou 'IND'"
		return(retour)
	}
	if(length(gen$FATHER) > 0 || length(gen$father) > 0) {
		if(length(gen$FATHER))
			names(gen)[names(gen) == "FATHER"] = "father"
	}
	else {
		retour$erreur = T
		retour$messageErreur = "Invalid 'gen' parameter: ascendance table must contain one column named 'father' or 'FATHER'"
			#"Parametre 'gen' invalide : la table d'ascendance doit contenir une colonne se nommant 'pere' ou 'PERE'"
		return(retour)
	}
	if(length(gen$MOTHER) > 0 || length(gen$mother) > 0) {
		if(length(gen$MOTHER) > 0)
			names(gen)[names(gen) == "MOTHER"] = "mother"
	}
	else {
		retour$erreur = T
		retour$messageErreur = "Invalid 'gen' parameter: ascendance table must contain one column named 'mother' or 'MOTHER'"
			#"Parametre 'gen' invalide : la table d'ascendance doit contenir une colonne se nommant 'mere' ou 'MERE'"
		return(retour)
	}
	if(length(gen$SEX) > 0 || length(gen$sex) > 0) {
		if(length(gen$SEX) > 0)
			names(gen)[names(gen) == "SEX"] = "sex"
	}
	retour$erreur = F
	retour$gen = gen
	return(retour)
}

gen.implex3V = function(ind, father, mother, pro = gen.pro(ind, father, mother), genNo = NULL, named = T)
{
	if(!gen.isGen3V(ind, father, mother))
		stop("at least one of ind, father or mother parameter is invalid")
		#stop("L'un des parametres ind, pere ou mere (au moins) est invalide")
	if(!is.null(genNo) && !is(genNo, "numeric"))
		stop("Invalid parameter: second parameter (genNo) must be an integer.")
		#stop("Parametre invalide : le deuxieme parametre (genNo) doit etre un entier.")
	if(!is.null(genNo))
		num = gen.initImp3V(ind, father, mother, pro, named = F)[genNo + 1]
	else {
		num = gen.initImp3V(ind, father, mother, pro, named = F)
		genNo <- 0:(length(num) - 1)
	}
	den = length(pro) * (2^genNo)
	complet = (100 * num)/den
	complet[is.na(complet)] = 0
	if(named)
		names(complet) <- genNo
	complet
}

gen.initImp3V = function(ind, father, mother, pro = gen.pro(ind, father, mother), named = T)
{
	if(!gen.isGen3V(ind, father, mother))
		stop("at least one of ind, father or mother parameter is invalid")
	if(!is.na(match(NA, match(pro, ind))))
		stop("One of the given proband is not part of the individuals list")
	taille = 0
	liste = pro
	nvx = ind
	# liste des individuals pas encore rencontrer
	retour = NULL
	while(length(liste) != 0) {
		taille = taille + 1
		retour[taille] = length(liste)
		tmp = NULL
		tmp = c(father[match(liste, ind)], mother[match(liste, ind)])
		tmp = intersect(tmp, nvx)
		nvx = setdiff(nvx, tmp)
		liste = tmp
	}
	if(named)
		names(retour) <- 0:(taille - 1)
	retour
}

gen.isGen3V = function(ind, father, mother)
{
	return(is.numeric(ind) && is.numeric(father) && is.numeric(mother) && length(ind) == length(father) && length(ind) == length(mother) &&
		all(is.element(father[father != 0], ind)) && all(is.element(mother[mother != 0], ind)) && 0 == intersect(mother, father))
}

gen.validationAsc = function(gen, pro = 0)
{
	asc = gen.genout(gen)
	if(sum(as.numeric(pro)) == 0)
		pro = gen.pro(gen)
	#VeRIFICATION DES ASCENDANCES
	#1. Verification s'il y existe des individuals egal a 0
	if(sum(asc$ind == 0) != 0) print("Error: some individuals are equal to 0 in the table d'ascendance")
	# print("Erreur: Il y a des individuals egal a 0 dans la table d ascendance")
	#2. Le nombre des fathers et des mothers == 0 doit etre identique
	#if(sum(asc$father == 0 & asc$mother == 0) != sum(asc$father == 0)) print(
	#"Erreur: Le nombre de couples 0 n egale pas le nombre de fathers 0"
	#)
	#if(sum(asc$father == 0 & asc$mother == 0) != sum(asc$mother == 0)) {
	# print("Erreur: Le nombre de couples = 0 n egale pas le nombre de mothers = 0"
	# )
	#}
	#if(sum(asc$father == 0) != sum(asc$mother == 0)) {
	# print("Erreur: Le nombre de fathers = 0 n egale pas le nombre de mothers = 0"
	#  )
	#}
	#3. Les parents doivent se retrouver dans la liste des individuals de la table d'ascendance
	if(sum(is.na(match(c(asc$father[asc$father != 0], asc$mother[asc$mother != 0]), asc$ind))) != 0) {
		print("Error: some parents are not found in the table d'ascendance")
		#print("Erreur: Il y a des parents qui ne retrouvent pas dans la table d ascendance")
	}
	#4. Les fathers et mothers ne peuvent pas etre les memes individuals ...bon c'est vrai aujourd'hui c'est possible!!! ;)
	if(sum(!is.na(match(asc$father[asc$father != 0], asc$mother[asc$mother != 0]))) != 0) {
		print("Error: some fathers and some mothers have the same individual number")
		#print("Erreur: Certains peres et certaines meres portent les memes no d indiviuds")
	}
	#5. Les ascendances doivent avoir une fin, on veut voir apparaitre le message "fin des ascendances"
	probands <- gen.pro(asc)
	pp <- unique(c(asc$father[match(probands, asc$ind)], asc$mother[match(probands, asc$ind)]))
	pp <- pp[pp != 0]
	for(i in 1:20) {
		pp <- unique(c(asc$father[match(pp, asc$ind)], asc$mother[match(pp, asc$ind)]))
		pp <- pp[pp != 0]
		if(length(pp == 0))
			break
	}
	#VeRIFICATIONS DES PROPOSANTS
	#6. Il ne doit pas y avoir de doublons de probands
	#if(length(pro$ind) != length(unique(pro$ind))) {
	if(length(pro) != length(unique(pro))) {
		print("Error: there are probands duplicates")
		#print("Erreur: Il existe des doublons de proposants")
	}
	#7. Les probands doivent faire partie des individuals de la table d'ascendance
	#if(sum(is.na(match(pro$ind, asc$ind))) != 0) {
	if(sum(is.na(match(pro, asc$ind))) != 0) {
		print("Error: some probands are not in the table d'ascendance")
		#print("Erreur: Il y a des proposants qui ne retrouvent pas dans la table d'ascendance")
	}
	#8. Les probands doivent etre les seuls individuals a ne pas avoir d'enfant dans la table d\'ascendance
	if(length(gen.pro(gen)) != length(pro)) {
		print("Error: probands are not the only ones without children in the table d'ascendance")
		#print("Erreur: Les proposants ne sont pas les seuls individuals a ne pas avoir d'enfant dans la table d'ascendance")
	}
}

gen.validationGen = function(gen, ..., check)
{
	#flag pour indiquer qu'un dataframe provient de l'utilisateur alors, il faudra mettre les etiquettes en minuscule
	gen.dataframe.utilisateur = 0
	#si gen est un objet GLgen, le dataframe correspondant a deja ete valide
	objet.glgen = 0
	#Verification de la nature de gen --
	#Peut etre un GLgen, dataframe (ind,father,mother) ou un vecteur (ind)
	#On s'assure que gen doit etre converti en dataframe
	if(is(gen, "GLgen")) {
		ret = gen.validationGLgen(gen)
		if(!ret$erreur) gen = gen.genout(gen)
		else			 return(ret)
		objet.glgen = 1
		#return(1)
	}
	else if(is(gen, "vector") && is.numeric(gen)) {
		#print("gen est vector et numeric")
		#on recherche les vecteurs father,mother,sex(optionnel)
		par... = gen.formatage...(...)
		posMere = match("mother", par...[[1]])
		posPere = match("father", par...[[1]])
		posSex = match("sex", par...[[1]])
		if(is.na(posMere) || is.na(posPere))
			return(list(erreur = T, messageErreur="Invalid '...' parameter: 'father' and 'mother' parameter names are obligatory"))
			#"Parametre '...' invalide : indication du nom des parametres 'pere' et 'mere' est obligatoire"))
		ind = gen
		father = par...[[2]][[posPere]]
		mother = par...[[2]][[posMere]]
		if(!(length(ind) == length(father) && length(ind) == length(mother)))
			return(list(erreur = T, messageErreur = "Invalid 'gen' parameter: ind, father and mother columns must have the same size"))
				#"Parametre 'gen' invalide : les colonnes ind, pere, mere doivent etre de meme taille"))
		if(is.na(posSex))
			gen = data.frame(ind = ind, father = par...[[2]][[posPere]], mother = par...[[2]][[posMere]])
		else {
			if(!(length(ind) == length(par...[[2]][[posSex]])))
				return(list(erreur = T, messageErreur = "Invalid 'gen' parameter: ind, father and mother columns must have the same size"))
					#"Parametre 'gen' invalide : les colonnes ind, pere, mere doivent etre de meme taille"))
			gen = data.frame(ind = ind, father = par...[[2]][[posPere]], mother = par...[[2]][[posMere]], sex = par...[[2]][[posSex]])
		}
	}
	else if(is(gen, "data.frame")) {
		#print("gen est data.frame")
		gen = gen.etiquetteGenMinuscule(gen)
		if(gen$erreur == T)
			return(gen)
		else if(gen$erreur == F)
			gen = gen$gen
		else return(list(erreur = T, messageErreur = "Invalid 'gen.etiquetteGenMinuscule' function return value"))
				#"La valeur de retour de la fonction 'gen.etiquetteGenMinuscule' n'est pas valide"))
		if(!(length(gen$ind) == length(gen$father) && length(gen$ind) == length(gen$mother)))
			return(list(erreur = T, messageErreur = "Invalid 'gen' parameter: ind, father and mother columns must have the same size"))
				#"Parametre 'gen' invalide : les colonnes ind, pere, mere doivent etre de meme taille"))
	}
	else
	 return(list(erreur = T, messageErreur = "Invalid 'gen' parameter: must be one of GLgen object, dataframe (ind,father,mother) or numeric vector"))
	 	#"Parametre 'gen' invalide : doit etre un objet GLgen, un dataframe (ind, pere, mere), ou un vecteur numerique (numeros d'individu)"))
	if(objet.glgen == 0) {
		ind = gen$ind
		father = gen$father
		mother = gen$mother
		if(!(is.numeric(ind) && is.numeric(father) && is.numeric(mother)))
			return(list(erreur = T, messageErreur = "Invalid 'gen' parameter: ind, father and mother columns must be numeric vectors"))
				#"Parametre 'gen' invalide : les colonnes ind, pere, mere doivent etre des vecteurs numeriques"))
		if(!(1 == length(intersect(mother, father))))
			return(list(erreur = T, messageErreur = "Invalid 'gen' parameter: identical individual number for both 'father' and 'mother'"))
				#"Parametre 'gen' invalide : un meme numero d'individu se retrouve dans les colonnes 'father' et 'mere'"))
		if(!is.element(2, check))
			if(!(all(is.element(father[father != 0], ind)) && all(is.element(mother[mother != 0], ind))))
				return(list(erreur = T, messageErreur = "Invalid 'gen' parameter: some father or mother number are not in the 'ind' column"))
					#"Parametre 'gen' invalide : des numeros de peres ou de meres ne se retrouvent pas dans la colonne 'ind'"))
		if(length(gen$sex) != 0) {
			if(!(length(gen$ind) == length(gen$sex)))
				return(list(erreur = T, messageErreur = "Invalid 'gen' parameter: sex column must have the same size as ind, father and mother"))
					#"Parametre 'gen' invalide : la colonne sex doit etre de meme taille que celle des colonnes ind, pere et mere"))
			sex = gen$sex
			tmp <- factor(sex, levels = c("H", "h", 1, "F", "f", 2))
			tmp2 <- as(tmp, "integer")
			tmp2[tmp2 == 2 | tmp2 == 3] <- 1
			tmp2[tmp2 == 4 | tmp2 == 5 | tmp2 == 6] <- 2
			if(any(is.na(tmp2)))
				return(list(erreur = T, messageErreur = "Invalid 'sex' parameter: bad data type"))
					#"Parametre 'sex' invalide : mauvais type de donnees"))
			if(is.element(2, check)) {
				if(all(tmp2[match(father[father != 0], ind, nomatch = T)] == 1) && all(tmp2[match(mother[mother != 0], ind,
					nomatch = T)] == 2))
					return(list(erreur = F, gen = gen))
				else return(list(erreur = T,messageErreur= "Error: some father or mother sex information do not concur with the 'sex' vector"))
						#"Erreur : le sex de certains individuals (meres ou peres) ne concorde pas avec les informations du vecteur 'sex'"))
			}
			else {
				if(all(tmp2[match(father[father != 0], ind)] == 1) && all(tmp2[match(mother[mother != 0], ind)] == 2))
					return(list(erreur = F, gen = gen))
				else return(list(erreur = T,messageErreur= "Error: some father or mother sex information do not concur with 'sex' vector"))
						#"Erreur : le sex de certains individuals (meres ou peres) ne concorde pas avec les informations du vecteur 'sex'"))
			}
		}
	}
	return(list(erreur = F, gen = gen))
}

gen.validationGLgen = function(object)
{
	#Verifie que 'object' est de la classe 'GLgen'
	if(!is(object, "GLgen")) return(list(erreur = T, messageErreur = "Invalid 'gen' parameter: must be a GLgen object"))
			#"Parametre 'gen' invalide : doit etre un objet GLgen"))
	
	#La signature (MD5) des genealogies est calculee et comparee avec celles contenues dans l'objet GLgen a valider.
	#prototype de la fonction c++ : ValidateGenealogie(long* Genealogie,long* isValid);
	isValid <- integer(1)
	ret <- try(.Call("SPLUSValidateGenealogie", object@.Data, isValid))

	isValid = ret$isValid
	if(class(ret) == "Error")
		return(list(erreur = T, messageErreur = "Error from dll return value")) #"Erreur de retour provenant de la dll"))
	if(as.logical(isValid))
		return(list(erreur = F, gen = object))
	else return(list(erreur = T, messageErreur = "Invalid 'gen' parameter: GLgen object altered after creation"))
			#"Parametre 'gen' invalide : objet GLgen modifie apres la creation"))
}

GLapplyCG = function(x, FUN, ..., FunReturnLength = 1, mirror = T, named = T, namesVector = NULL)
{
	#Si c'est une matrice applique la fonction a la matrice entiere
	#Si c'est un vecteur applique la fonction au vecteur en entier
	#Si c'est un matrice avec group ... 
	#applique la fonction a chaque group de proband et pour chaque ancestor
	#Bon determine si c'est pour plusieurs depth
	#Group
	if(is.null(group(x))) {
		#Non pas de groups
		#applique la formule pour chaque ancestor (chaque colonne)
		x <- as(x, "matrix")
		if(!named)
			dimnames(x) <- NULL
		ret <- apply(x, 2, function(x, FUN2, ...)
		{
			FUN2(x, ...)
		}
		, FUN2 = FUN, ...)
		dim(ret) <- c(FunReturnLength, dim(x)[2])
		if(named)
			dimnames(ret) <- list(namesVector, dimnames(x)[[2]])
		drop(ret)
	}
	else {
		#Oui des groups 
		#  si FUN retourne un vecteur alors retourne tableau 3 dim
		#si FUN retourne valeur alors retourne une matrice
		xgrindex <- x@grindex
		xgroupe <- x@group
		x <- as(x, "matrix")
		ng <- length(xgroupe)
		#nombre de group
		ret <- array(0., c(ng, dim(x)[2], FunReturnLength))
		#Calcul de la valeur (probablement moyen de faire plus efficace)    
		for(a in 1:ng) {
			ret[a,  ,  ] <- sapply(1:dim(x)[2], function(x, mat, FUN2, ...)
			{
				FUN2(mat[, x], ...)
			}
			, mat = x[xgrindex[[a]],  , drop = F], FUN2 = FUN, ...)
		}
		if(named)
			dimnames(ret) <- list(names(xgroupe), dimnames(x)[[2]], namesVector)
		drop(ret)
	}
}

GLapplyF = function(x, FUN, ..., FunReturnLength = 1, mirror = T, named = T, namesVector = NULL)
{
	#Si c'est un matrice, applique la fonction a la matrice entiere
	#Si c'est un vecteur, applique la fonction au vecteur entier
	#Si c'est un vecteur avec group ... applique la fonction a chaque sous-vecteur de group....
	#Et ce pour chaque depth si necessaire
	if(!is.null(depth(x))) {
		if(is.null(group(x))) {
			#Pas de group 
			#Un fois la fonction par depth (resultat est un gugus....) 
			depth <- x@depth
			x <- as(x, "array")
			#Boucle pour chaque depth
			ret <- apply(x, MARGIN = length(dim(x)), function(x, FUN2, ...)
			{
				FUN2(x, ...)
			}
			, FUN2 = FUN, ...)
			#Formatage
			if(is.matrix(ret) && named) dimnames(ret) <- list(namesVector, NULL)
			#elimine dimension inutile & construit l'objet approprier
			GLmulti(drop(ret), depth)
		}
		else {
			#Oui des groups
			#Retourne un array de dimension 4 (gr1,gr2,vecteurresult,depth) 
			xgrindex <- x@grindex
			xgroupe <- x@group
			ng <- length(xgroupe)
			#nombre de group
			depth <- x@depth
			x <- as(x, "array")
			#Creation du tableau original
			ng <- length(xgroupe)
			#nombre de group
			ret <- array(0., c(ng, FunReturnLength, length(depth)))
			#Calcul de la valeur (probablement moyen de faire plus efficace)    
			for(p in 1:length(depth))
				for(a in 1:ng) {
					ret[a,  , p] <- FUN(x[xgrindex[[a]], p, drop = F], ...)
				}
			if(named)
				dimnames(ret) <- list(names(xgroupe), namesVector, NULL)
			#elimine dimension inutile & construit l'objet approprier
			GLmulti(drop(ret), depth)
		}
	}
	else {
		#Une depth
		#Group
		if(is.null(group(x))) {
			#Pas de group (simple un resultat)
			x <- as(x, "array")
			ret = FUN(x, ...)
			if(named)
				names(ret) <- namesVector
			return(ret)
		}
		else {
			#Oui des groups 
			#  si FUN retourne un vecteur alors retourne tableau 3 dim
			#si FUN retourne valeur alors retourne une matrice
			xgrindex <- x@grindex
			xgroupe <- x@group
			x <- as(x, "array")
			ng <- length(xgroupe)
			#nombre de group
			ret <- array(0., c(ng, FunReturnLength))
			#Calcul de la valeur (probablement moyen de faire plus efficace)    
			for(a in 1:ng)
				ret[a,  ] <- FUN(x[xgrindex[[a]], drop = F], ...)
			if(named)
				dimnames(ret) <- list(names(xgroupe), namesVector)
			return(drop(ret))
		}
	}
}

GLapplyGroup = function(Matrice, GroupIndex, FUN, ..., named = T)
{
	#Applique une fonction a une grande matrice en fct du group
	Matrice = unclass(Matrice)
	#Creation de la matrice resultat
	n <- length(GroupIndex)
	ret <- matrix(0, ncol = n, nrow = n)
	#Calcul de la valeur (probablement moyen de faire plus efficace)
	for(a in 1:n)
		for(b in a:n) {
			ret[a, b] <- FUN(Matrice[GroupIndex[[a]], GroupIndex[[b]], drop = F], ...)
			#Mirroir
			ret[b, a] <- ret[a, b]
		}
	if(named)
		dimnames(ret) <- list(names(GroupIndex), names(GroupIndex))
	return(ret)
}

GLapplyPhi = function(x, FUN, ..., FunReturnLength = 1, mirror = T, named = T, namesVector = NULL)
{
	#Si c'est une matrice appliquer la fonction a la matrice entiere
	#Si c'est un vecteur appliquer la fonction au vecteur en entier
	#Si c'est un matrice avec group ... applique la fonction a chaque sous-matrice entre group....
	#Et ce pour chaque depth, si necessaire
	if(!is.null(depth(x))) {
		if(is.null(group(x))) {
			depth <- x@depth
			x <- as(x, "array")
			#Boucle pour chaque depth
			ret <- apply(x, MARGIN = length(dim(x)), function(x, FUN2, ...)
			{
				FUN2(x, ...)
			}
			, FUN2 = FUN, ...)
			#Formatage
			if(is.matrix(ret) && named) dimnames(ret) <- list(namesVector, NULL)
			#elimine dimension inutile & construit l'objet approprie
			GLmulti(drop(ret), depth)
		}
		else {
			#Retourne un array de dimension 4 (gr1,gr2,vecteurresult,depth) 
			xgrindex <- x@grindex
			xgroupe <- x@group
			ng <- length(xgroupe)
			#nombre de group
			depth <- x@depth
			x <- as(x, "array")
			#Creation du tableau original
			ng <- length(xgroupe)
			#nombre de group
			ret <- array(0., c(ng, ng, FunReturnLength, length(depth)))
			#Calcul de la valeur (probablement moyen de faire plus efficace)    
			for(p in 1:length(depth))
				for(a in 1:ng)
					for(b in a:ng) {
						ret[a, b,  , p] <- FUN(as(x[xgrindex[[a]], xgrindex[[b]], p], "matrix"), ...)
						if(mirror)
							ret[b, a,  , p] <- ret[a, b,  , p]
						else ret[b, a,  , p] <- FUN(as(x[xgrindex[[b]], xgrindex[[a]], p], "matrix"), ...)
					}
			if(named)
				dimnames(ret) <- list(names(xgroupe), names(xgroupe), namesVector, NULL)
			#elimine dimension inutile & construit l'objet approprie
			GLmulti(drop(ret), depth)
		}
	}
	else {
		if(is.null(group(x))) {
			#Pas de group (simple un resultat)
			x <- as(x, "array")
			ret = FUN(x, ...)
			if(named)
				names(ret) <- namesVector
			ret
		}
		else {
			#Oui des groups 
			#si FUN retourne un vecteur, alors retourne tableau 3 dim
			#si FUN retourne valeur, alors retourne une matrice
			xgrindex <- x@grindex
			xgroupe <- x@group
			x <- as(x, "array")
			ng <- length(xgroupe)
			#nombre de group
			ret <- array(0., c(ng, ng, FunReturnLength))
			#Calcul de la valeur (probablement moyen de faire plus efficace)    
			for(a in 1:ng)
				for(b in a:ng) {
					ret[a, b,  ] <- FUN(x[xgrindex[[a]], xgrindex[[b]], drop = F], ...)
					if(mirror)
						ret[b, a,  ] <- ret[a, b,  ]
					else ret[b, a,  ] <- FUN(x[xgrindex[[b]], xgrindex[[a]], drop = F], ...)
				}
			if(named)
				dimnames(ret) <- list(names(xgroupe), names(xgroupe), namesVector)
			#elimine dimension inutile
			drop(ret)
		}
	}
}

GLapplyPhi.mat = function(x, FUN, ..., FunReturnLength = 1, mirror = T, named = T, namesVector = NULL)
{
	#Applique a une matrice de phi ayant une seule depth
	#Pas de group
	x <- as(x, "array")
	ret = FUN(x, ...)
	if(named)
		names(ret) <- namesVector
	return(ret)
}

GLCGGroup = function(MatriceCG, Group, proband)
{
	#un cg group c'est une matrice de group x ancestor...
	if(!is(Group, "GLgroup")) stop("Invalid parameter: Group must a group of valid probands")
	#TROUVE LES PROPOSANTS CORRESPONDANTS A LA MATRICE A L'AIDE DES ETIQUETTES
	if(missing(proband)) proband <- as.integer(dimnames(MatriceCG)[[1]])
	#La dimension 1 c'est les ancestors
	if(!is.numeric(proband) || length(proband) != dim(MatriceCG)[[1]]) 
		stop("Invalid parameter: proband must be the probands list used to generate the MatricePhi")
			#"Parametre invalide: proposant doit etre la liste de proposants utilises pour generer la MatricePhi")
	if(length(dim(MatriceCG)) == 2) {
		#COMPACTAGE
		ind <- unique(unlist(Group, use.names = T))
		ind2 <- match(ind, proband)
		if(any(is.na(ind2)))
			stop("Invalid parameter: all probands of Group must be part of the probands")
			#stop("Parametre invalide: tous les probands du Group doit faire partie des probands")
		MatriceCG <- MatriceCG[ind2,  , drop = F]
		#pro,anc (tous les ancestors sont conserve)
		#GENERATION DES INDICES DES PROPOSANTS
		indice <- lapply(Group, function(x, pro)
		{
			match(x, pro)
		}
		, pro = ind)
		#CREATION D'OBJET
		new("GLCGMatrixGroupSingle", MatriceCG, group = Group, grindex = indice)
	}
	else stop("Invalid parameter: MatriceCG must one or more phi matrix or a GLCGMatrixGroupSingle object")
		#stop("Parametre invalide: MatriceCG doit-etre une ou plusieurs matrices phi ou un objet GLCGMatrixGroupSingle")
}

GLFGroup = function(VecteurF, Group, depth = NULL, proband)
{
	#Prend.. vecteur ou matrix x par 1   -> GLmultiFgroupSingle
	#Prend.. matrice x par y et depth!= y     -> GLmultiFgroup
	#Prend.. GLmultiVector          -> GLmultiFgroup  
	#Plus facile a gerer en matrice
	VecteurF <- as.matrix(VecteurF)
	if(is.numeric(VecteurF) && is.matrix(VecteurF) && dim(VecteurF)[2] > 1) {
		#Plusieurs depths
		if(is.null(depth) && !is.null(depth(VecteurF))) depth = depth(VecteurF)
		VecteurF <- as(VecteurF, "matrix")
		#Compactage
		ind <- unique(unlist(Group, use.names = T))
		ind2 <- match(ind, proband)
		if(any(is.na(ind2)))
			stop("Invalid grpPro parameter (function GLFGroup): all probands of 'grpPro' must be part of probands")
			#stop("Parametre 'grpPro' invalide (Fct : GLFGroup): tous les probands de 'grpPro' doivent faire partie des probands")
		VecteurF <- VecteurF[ind2,  , drop = F]
		#Generation des indices des probands
		indice <- lapply(Group, function(x, pro)
		{
			match(x, pro)
		}
		, pro = ind)
		#Creation de l'objet
		new("GLmultiFGroup", GLmulti(VecteurF, depth), group = Group, grindex = indice)
	}
	else if(dim(VecteurF)[2] == 1) {
		#C'est un vecteur
		#une depth, compactage
		ind <- unique(unlist(Group, use.names = T))
		ind2 <- match(ind, proband)
		if(any(is.na(ind2)))
			stop("Invalid grpPro parameter (function GLFGroup): all probands of 'grpPro' must be part of probands")
			#stop("Parametre 'grpPro' invalide (Fct : GLFGroup): tous les probands de 'grpPro' doivent faire partie des probands")
		VecteurF <- VecteurF[ind2,  , drop = F]
		#Generation des indices des probands
		indice <- lapply(Group, function(x, pro)
		{
			match(x, pro)
		}
		, pro = ind)
		#Creation de l'objet
		new("GLmultiFGroupSingle", VecteurF, group = Group, grindex = indice)
	}
	else	stop("Invalid 'VecteurF' parameter (function GLFGroup): must be one or more F vectors or a 'GLmultiFGroup' object")
		#stop("Parametre 'VecteurF' invalide (Fct : GLFGroup): doit etre un ou plusieurs vecteurs F ou un objet 'GLmultiFGroup'")
}

GLgen = function(...)
{
	gen.genealogy(...)
}

is.all.white <- function(listeNoms)
{
  unlist(lapply(listeNoms,function(n){gsub(" ", "", n, fixed=T)==""}))
}

GLgroup = function(liste)
{
	#creation du group et complementation
	if(!is.list(liste)) stop("Invalid parameter: liste must be a valid list") #stop("parametre invalide : liste doit etre une liste valide")
	defaultname <- sapply(1:length(liste), function(i)
	paste("Group", i))
	if(is.null(names(liste)))
		toreplace <- rep(T, length(liste))
	else toreplace <- is.all.white(names(liste)) #is.all.white(names(liste), empty = T)
	names(liste)[toreplace] <- defaultname[toreplace]
	return(new("GLgroup", liste))
}

GLmulti = function(Array, depth, drop = T, addDim = F)
{
	#construit un objet GLmultiMatrix ou GLmultiVector a l'aide de l'array fournie
	#trouve la dimension correspondante
	#Si Array n'est pas une liste 
	if(!is.list(Array)) {
		dimen <- length(dim(Array))
		if(length(depth) > 1 || drop == F) {
			#s'il y a plus d'une depth alors c'est un objet GLmulti 
			if(drop == T && is.array(Array) && dim(Array)[dimen] != length(depth)) 
				stop("Error: depth size must be the Array's last dimension")
				#stop("Erreur: la taille de depth doit correspondre \340 la derniere dimension d'Array")
			if(drop == T && !is.array(Array) && length(Array) != length(depth))
				stop("Error: depth size must be the Array's last dimension")
				#stop("Erreur: la taille de depth doit correspondre \340 la derniere dimension d'Array")
			if(drop == F && addDim == T && length(depth) == 1) {
				if(is.null(dim(Array))) {
					tmpnom <- list(names(Array))
					names(Array) <- NULL
					dim(Array) <- c(length(Array), 1)
				}
				else {
					tmpnom <- dimnames(Array)
					dimnames(Array) <- NULL
					dim(Array) <- c(dim(Array), 1)
				}
				dimnames(Array) <- c(tmpnom, list(NULL))
				dimen <- length(dim(Array))
			}
			if(dimen == 4)
				return(new("GLmultiArray4", Array, depth = depth))
			else if(dimen == 3)
				return(new("GLmultiMatrix", Array, depth = depth))
			else if(dimen == 2)
				return(new("GLmultiVector", Array, depth = depth))
			else if(dimen == 0 || dimen == 1 || is.null(dimen)) {
				#dans ce cas, il n'y a rien a retourne sauf des depth
				tmp <- new("GLmultiNumber", as.numeric(Array), depth = depth)
				#Verifie s'il y a des nom 
				if(is(Array, "named")) new("GLmultiNumber", as.numeric(Array), depth = depth, .Names = names(
						Array)) else new("GLmultiNumber", as.numeric(Array), depth = depth)
			}
			else stop("Impossible to create the object with the given parameters")
				#stop("Impossible de creer l'objet avec les parametres specifies")
		}
		else return(Array)
	}
	else {
		#Si Array est une liste
		#Ici qu'il faudrait creer un objet GLmultiList
		#S'il y a plus de 1 depth
		if(length(depth) > 1) {
			names(Array) <- paste(rep("Gener ", length(Array)), as.character(c(depth)), sep = "")
			return(new("GLmultiList", Array))
		}
		else return(Array)
	}
}

GLnoone = -999

GLOverGroup2 = function(dim1, dim2, ..., named)
{
	if(missing(dim1))
		dim1 <- GLnoone
	if(missing(dim2))
		dim2 <- GLnoone
	if(missing(named)) {
		n <- T
		p <- nargs()
	}
	else {
		n <- named
		p <- nargs() - 1
	}
	return(list(dim1 = dim1, dim2 = dim2, param = p, named = n))
}

GLOverGroup3 = function(pro, dim1, dim2, ..., named, abs)
{
     #print(paste( missing(pro),missing(dim1),missing(dim2),missing(named),missing(abs) ))
     #print(nargs())
	if(missing(pro))
		pro <- GLnoone
	if(missing(dim1))
		dim1 <- GLnoone
	if(missing(dim2))
		dim2 <- GLnoone
	if(missing(named)) {
		n <- T
		p <- nargs()
	}
	else {
		n <- named
		p <- nargs() - 1
	}
	if(missing(abs)) {
		ab <- F
	}
	else {
		ab <- abs
		p <- p - 1
	}
	return(list(pro = pro, dim1 = dim1, dim2 = dim2, param = p, named = n, abs = ab))
}

GLOverGroup4 = function(pro, dim1, dim2, dim3, ..., named, abs)
{
	if(missing(pro))
		pro <- GLnoone
	if(missing(dim1))
		dim1 <- GLnoone
	if(missing(dim2))
		dim2 <- GLnoone
	if(missing(dim3))
		dim3 <- GLnoone
	if(missing(named)) {
		n <- T
		p <- nargs()
	}
	else {
		n <- named
		p <- nargs() - 1
	}
	if(missing(abs)) {
		ab <- F
	}
	else {
		ab <- abs
		p <- p - 1
	}
	list(pro = pro, dim1 = dim1, dim2 = dim2, dim3 = dim3, param = p, named = n, abs = ab)
}

GLOverNumber1 = function(pro, ..., named, abs)
{
	if(missing(pro))
		pro <- GLnoone
	if(missing(named)) {
		n <- T
		p <- nargs()
	}
	else {
		n <- named
		p <- nargs() - 1
	}
	if(missing(abs)) {
		ab <- F
	}
	else {
		ab <- abs
		p <- p - 1
	}
	list(pro = pro, param = p, named = n, abs = ab)
}

GLOverVector2 = function(pro, dim1, ..., named, abs)
{
	if(missing(pro))
		pro <- GLnoone
	if(missing(dim1))
		dim1 <- GLnoone
	if(missing(named)) {
		n <- T
		p <- nargs()
	}
	else {
		n <- named
		p <- nargs() - 1
	}
	if(missing(abs)) {
		ab <- F
	}
	else {
		ab <- abs
		p <- p - 1
	}
	return(list(pro = pro, dim1 = dim1, param = p, named = n, abs = ab))
}

GLPhiGroup = function(MatricePhi, Group, depth = NULL, proband)
{
	#Nombre de depth
	if(is.numeric(MatricePhi) && is.array(MatricePhi) && length(dim(MatricePhi)) == 3) {
		#plusieurs depths
		if(is.null(depth) && !is.null(depth(MatricePhi))) depth = depth(MatricePhi)
		MatricePhi <- as(MatricePhi, "array")
		#compactage
		ind <- unique(unlist(Group, use.names = T))
		ind2 <- match(ind, proband)
		if(any(is.na(ind2)))
			stop("Invalid parameter: all probands of Group must be part of proband")
			#stop("Parametre invalide: tous les probands du Group doit faire partie des probands")
		MatricePhi <- MatricePhi[ind2, ind2,  ]
		#Generation des indices des probands
		indice <- lapply(Group, function(x, pro)
		{
			match(x, pro)
		}
		, pro = ind)
		#Creation d'objet
		new("GLmultiPhiGroup", GLmulti(MatricePhi, depth), group = Group, grindex = indice)
	}
	else if(length(dim(MatricePhi)) == 2) {
		#compactage
		ind <- unique(unlist(Group, use.names = T))
		ind2 <- match(ind, proband)
		if(any(is.na(ind2)))
			stop("Invalid parameter: all probands of Group must be part of proband")
			#stop("Parametre invalide: tous les probands du Group doit faire partie des probands")
		MatricePhi <- MatricePhi[ind2, ind2]
		#Generation des indices des probands
		indice <- lapply(Group, function(x, pro)
		{
			match(x, pro)
		}
		, pro = ind)
		#Creation d'objet
		new("GLmultiPhiGroupSingle", MatricePhi, group = Group, grindex = indice)
	}
	else stop("Invalid 'MatricePhi' parameter: must be a phi matrix with one or more depths")
		#stop("Parametre 'MatricePhi' invalide: doit etre un matrice phi avec une ou plusieurs depths")
}

GLPriv.completeness3V = function(ind, father, mother, pro, genNo, named)
{
	num = GLPriv.initcomp3V(ind, father, mother, pro, named = F)
	num = num[genNo + 1]
	den = length(pro) * (2^genNo)
	complet = (100 * num)/den
	complet[is.na(complet)] = 0
	if(named)
		names(complet) <- genNo
	return(complet)
}

GLPriv.entropie3V = function(ind, father, mother, pro)
{
	#Nombre de fondateurs
	vctF <- GLPriv.initfon3V(ind, father, mother, pro, named = F)
	#Nombre de demi-fondateurs
	vctDF <- GLPriv.initDemifon3V(ind, father, mother, pro, named = F)
	genNoF <- 0:(length(vctF) - 1)
	#Nombre de generations des fondateurs
	genNoDF <- 0:(length(vctDF) - 1)
	#Nombre de generations des demi-fondateurs
	#Calcul de l'entropie
	vctComp <- sum((genNoF * vctF)/(length(pro) * (2^genNoF))) + sum(((genNoDF * vctDF)/(length(pro) * (2^genNoDF))) * 0.5)
	#Si le calcul de l'entropie est vide, 0 lui est assigne
	if(is.na(vctComp)) vctComp <- 0
	#Retourne le resultat
	return(vctComp)
}

error.bar <- function(x, y, type, lty, lwd, col, upper, lower=upper, length=0.1, ...){
    if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper)) stop("vectors must be same length")
    arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, type=type, lty = lty, lwd=lwd, col=col, ...)
}

GLPriv.implex3V = function(ind, father, mother, pro, genNo, named)
{
	num = GLPriv.initImp3V(ind, father, mother, pro, named = F)
	num = num[genNo + 1]
	den = length(pro) * (2^genNo)
	complet = (100 * num)/den
	complet[is.na(complet)] = 0
	if(named)
		names(complet) <- genNo
	return(complet)
}

GLPriv.initcomp3V = function(ind, father, mother, pro, named)
{
	taille = 0
	liste = pro
	retour = NULL
	while(length(liste) != 0) {
		taille = taille + 1
		retour[taille] = length(liste)
		tmp = NULL
		tmp = c(father[match(liste, ind)], mother[match(liste, ind)])
		tmp = tmp[(tmp != 0)]
		liste = tmp
	}
	if(named)
		names(retour) <- 0:(taille - 1)
	return(retour)
}

GLPriv.initDemifon3V = function(ind, father, mother, pro = gen.pro(ind, father, mother), named = T)
{
	#Validation des parametres (a mettre dans gen.detectionErreur)
	if(!gen.isGen3V(ind, father, mother)) stop("At least one of the ind, father or mother parameter is invalid")
		#stop("L'un des parametres ind, father ou mother (au moins) est invalide.")
	if(!all(is.element(pro, ind))) stop("One of the given probands is not in the individuals list")
		#stop("L'un des probands fourni ne fait pas parti de la liste d'individuals")
	#Trouve les demi-fondateurs a chaque generation
	taille <- 0
	liste <- pro
	retour <- rep(0, length(ind))
	while(length(liste) != 0) {
		taille <- taille + 1
		a <- length(liste[is.na(liste)])
		liste <- liste[!is.na(liste)]
		b <- b <- length(c((liste[father[match(liste, ind)] != 0 & mother[match(liste, ind)] == 0]), (liste[father[match(liste,
			ind)] == 0 & mother[match(liste, ind)] != 0])))
		retour[taille] = a + b
		tmp <- NULL
		tmp <- c(father[match(liste, ind)], mother[match(liste, ind)])
		tmp <- tmp[(tmp != 0)]
		liste <- tmp
	}
	#Renvois le resultat dans un vecteur
	retour <- retour[1:taille]
	if(named)
		names(retour) <- 0:(taille - 1)
	retour
}

GLPriv.initfon3V = function(ind, father, mother, pro, named)
{
	taille = 0
	liste = pro
	retour = rep(0, length(ind))
	while(length(liste) != 0) {
		taille = taille + 1
		a = length(liste[is.na(liste)])
		liste = liste[!is.na(liste)]
		b = length(liste[father[match(liste, ind)] == 0 & mother[match(liste, ind)] == 0])
		retour[taille] = a + b
		tmp = NULL
		tmp = c(father[match(liste, ind)], mother[match(liste, ind)])
		tmp = tmp[(tmp != 0)]
		liste = tmp
	}
	retour = retour[1:taille]
	if(named)
		names(retour) <- 0:(taille - 1)
	return(retour)
}

GLPriv.initImp3V = function(ind, father, mother, pro, named)
{
	taille = 0
	liste = pro
	retour = NULL
	while(length(liste) != 0) {
		taille = taille + 1
		retour[taille] = length(unique(liste))
		tmp = NULL
		tmp = c(father[match(liste, ind)], mother[match(liste, ind)])
		tmp = tmp[(tmp != 0)]
		liste = tmp
	}
	if(named)
		names(retour) <- 0:(taille - 1)
	return(retour)
}

GLPriv.variance3V = function(ind, father, mother, pro)
{
	#Nombre de fondateurs
	vctF <- GLPriv.initfon3V(ind, father, mother, pro, named = F)
	#Nombre de demi-fondateurs
	vctDF <- GLPriv.initDemifon3V(ind, father, mother, pro, named = F)
	genNoF <- 0:(length(vctF) - 1)
	#Nombre de generations des fondateurs
	genNoDF <- 0:(length(vctDF) - 1)
	#Nombre de generations des demi-fondateurs
	#Calcul l'entropie
	P <- sum((genNoF * vctF)/(length(pro) * (2^genNoF))) + sum(((genNoDF * vctDF)/(length(pro) * (2^genNoDF))) * 0.5)
	#Calcule de la variance de l'entropie
	varP <- (sum((genNoF^2 * vctF)/(length(pro) * (2^genNoF))) + sum((genNoDF^2 * vctDF)/(length(pro) * (2^genNoDF)) * 0.5)) -
		(P^2)
	#Returne le resultat
	return(varP)
}

GLPrivCG = function(gen, pro, ancestors, print.it = F, named = T)
{
	#Structure necessaire pour emmagasiner le resultat la fonction de la dll
	tmp <- double(length(ancestors) * length(pro))
	#Call de la fonction en C
	.Call("SPLUSConGen", gen@.Data, pro, length(pro), ancestors, length(ancestors), tmp, print.it, specialsok = T)	
	#Creation de la matrice de resultat
	dim(tmp) <- c(length(pro), length(ancestors))
	if(named)
		dimnames(tmp) <- list(pro, ancestors)
	if(print.it) {
		argument <- c(deparse(substitute(gen)), deparse(substitute(pro)), deparse(substitute(ancestors)))
		header.txt <- paste("\n   ***   Calls : gen.GC(", argument[1], ",", argument[2], ",", argument[3], ")  ***\n\n")
		cat(header.txt)
	}
	return(invisible(tmp))
}

GLPrivCGcumul = function(CG, named = T)
{
	#Calcule la somme par group
	somme <- GLapplyCG(CG, sum, named = named)
	if(!is.matrix(somme))
		return(cumsum(rev(sort(somme))))
	else return(t(apply(somme, 1, function(x)
		cumsum(rev(sort(x))))))
}

GLPrivCGgroup = function(CG, grppro, pro)
{
	#Applique un group a une matrice CG
	#Accepte: Une Matrice et un group
	#Dans tous les cas, pro peut-etre omis si CG est etiquette dans le cas contraire
	#il faut fournir une liste de proband de meme taille que la premiere dimension de CG 
	if(missing(pro)) GLCGGroup(CG, grppro) else GLCGGroup(CG, grppro, proband = pro)
}

GLPrivCGmoyen = function(CG, named = T)
{
	GLapplyCG(CG, mean, named = named)
}

GLPrivCGPLUS = function(gen, pro, ancestors, vctProb, print.it = F, named = T)
{
	#Structure necessaire pour emmagasiner le resultat la fonction de la dll
	tmp <- double(length(ancestors) * length(pro))
	#Call de la fonction en C
	.Call("SPLUSConGenPLUS", gen@.Data, pro, length(pro), ancestors, length(ancestors), vctProb, tmp, print.it, specialsok = T)
	#Creation de la matrice de resultat
	dim(tmp) <- c(length(pro), length(ancestors))
	if(named)
		dimnames(tmp) <- list(pro, ancestors)
	if(print.it) {
		argument <- c(deparse(substitute(gen)), deparse(substitute(pro)), deparse(substitute(ancestors)))
		header.txt <- paste("\n   ***   Calls : gen.GC(", argument[1], ",", argument[2], ",", argument[3], ")  ***\n\n")
		cat(header.txt)
	}
	return(invisible(tmp))
}

GLPrivCGproduit = function(CG, named = T)
{
	if(is.null(group(CG))) {
		#La formule   
		resultat <- drop(exp(t(rep(1, dim(CG)[1])) %*% log(CG)))
		if(named)
			names(resultat) <- dimnames(CG)[[2]]
		#else names(resultat) <- NULL
		return(resultat)
	}
	else {
		#Oui des groups 
		xgrindex <- CG@grindex
		xgroupe <- CG@group
		x <- as(CG, "matrix")
		ng <- length(xgroupe)
		#nombre de group
		#Matrice de resultat
		ret <- array(0., dim = c(ng, dim(x)[2]))
		#Calcul de la valeur pour chaque group
		for(a in 1:ng) {
			#Trouve la matrice sur laquel applique la fonction
			ind = xgrindex[[a]]
			mat = x[ind,  , drop = F]
			#La formule
			ret[a,  ] <- drop(exp(t(rep(1, length(ind))) %*% log(mat)))
		}
		if(named)
			dimnames(ret) <- list(names(xgroupe), dimnames(x)[[2]])
		return(drop(ret))
	}
}

GLPrivCGtotal = function(CG, named = T)
{
	GLapplyCG(CG, sum, named = named)
}

GLPrivExtCG = function(x, ..., drop)
{
	#extrait pour certaine depth
	l2 <- GLOverGroup2(...)
	dim1 <- l2$dim1
	#Group
	dim2 <- l2$dim2
	#Ancestor
	if(missing(drop)) drop <- T
	#fin if un parametre
	#Valeur par defaults
	if(l2$param < 3) {
		#Un seul param        
		if(is(dim1, "GLnothing")) class(dim1) <- "missing"
		if(is(dim2, "GLnothing"))
			class(dim2) <- "missing"
		if(drop) {
			#Si drop=T alors retourne la matrice de resultat
			#Genere le tableau de CGmoyen
			m <- GLapplyCG(x, mean, named = l2$named)
			#Peut-importe dim1 et dim2 ca passe a l'operateur
			if(l2$param == 1) getMethod("[", "array")(m, dim1, drop = drop) else getMethod("[", "array")(m, dim1, dim2,
					drop = drop)
		}
		else {
			#Si drop = F return GLCGMatrixGroupSingle
			xgrindex <- x@grindex
			xgroupe <- x@group
			x <- unclass(x)
			if(is(dim1, "missing"))
				dim1 <- 1:length(xgroupe)
			#validation
			tmp <- 1:length(xgroupe)
			names(tmp) <- names(xgroupe)
			lg <- tmp[dim1]
			#Nouveau Group a construire
			if(any(is.na(lg))) stop("Invalid extraction: one of the element was not found")
				#stop("Extraction invalide, un des elements n'a pas ete trouve")
			#extraction
			lg <- xgroupe[dim1]
			#Regeneration de la liste de proband
			gind <- unlist(xgrindex, use.names = F)
			ggrou <- unlist(xgroupe, use.names = F)
			pro <- ggrou[match(1:(dim(x)[1]), gind)]
			if(any(is.na(pro))) stop("Can not use Drop=T for this particular object")
				#stop("Vous ne pouvez pas utilise Drop=T pour cette objet en particulier")
			#Nouvelle matrice reduite
			class(lg) <- "GLgroup"
			#Creation de l'objet 
			#dim2 est le subscript pour les ancestors
			GLCGGroup(x[, dim2, drop = F], lg, proband = pro)
		}
	}
	else stop("You can only use one or two subscript (just like a matrix)")
		#stop("Vous ne pouvez utiliser qu'un ou deux subscript (exactement comme une matrice)")
}

GLPrivExtF = function(x, ..., drop)
{
	#Ressemble beaucoup a une GLmultiVector.... (C'en est un en fait)
	#extrait pour certaine depth
	#list(pro=pro,dim1=dim1,param=p,named=n,abs=ab)
	l2 <- GLOverVector2(...)
	pro <- l2$pro
	#depth
	dim1 <- l2$dim1
	#Group
	if(missing(drop)) drop <- T
	#Valeur par defaults
	if(is(dim1, "GLnothing")) class(dim1) <- "missing"
	if(drop) {
		if(l2$param == 0 || l2$param == 1 || l2$param == 2) {
			#Si drop=T alors retourne un GLMultiMatrix de resultat
			#declassement
			xgrindex <- x@grindex
			xgroupe <- x@group
			xdepth <- x@depth
			x <- as(x, "array")
			ng <- length(xgroupe)
			#nombre de group   
			#Genere le tableau de phi moyen pour chaque depth 
			m <- array(0., c(ng, length(xdepth)))
			#Calcul de la valeur (probablement moyen de faire plus efficace)    
			for(p in 1:length(xdepth))
				for(a in 1:ng) {
					m[a, p] <- mean(x[xgrindex[[a]], p])
				}
			if(l2$named)
				dimnames(m) <- list(names(xgroupe), NULL)
			m <- drop(m)
			#Fin du calcul de m         
			#Cas particulier depth
			if(is(pro, "GLnothing")) class(pro) <- "missing" else if(is.numeric(pro) && !l2$abs) {
				#Si c'est pas valeur absolue    
				#Recherche comme prevu dans la liste de depth
				pos <- match(pro, xdepth)
				if(any(is.na(pos)))
					stop(cat("Some depths were not found: ", pro[is.na(pos)], "\n"))
					#stop(cat("Certaine(s) depth(s) demandees n'ont pas put etre trouvees :", pro[is.na(pos)], "\n"))
				pro <- pos
			}
			#Else dans ce cas garde pro comme il etait      
			if(is(dim1, "GLnothing")) class(dim1) <- "missing"
			#Peut-importe dim1 et dim2 ca passe a l'operateur 
			GLmulti(getMethod("[", "array")(m, dim1, pro, drop = T), xdepth[pro])
		}
		else stop("You can only use one or two subscript (depth alone or depth, group)")
			#stop("Vous ne pouvez utiliser qu'un ou deux subscripts (depth seul ou depth,group)")
	}
	else {
		#Si drop = F alors retourne un objet GLmultiPhiGroupSingle modifier
		xgrindex <- x@grindex
		xgroupe <- x@group
		xdepth <- x@depth
		if(l2$param == 0 || l2$param == 1 || l2$param == 2) {
			#Cas particulier depth
			if(is(pro, "GLnothing")) class(pro) <- "missing" else if(is.numeric(pro) && !l2$abs) {
				#Si c'est pas valeur absolue
				#Recherche comme prevu dans la liste de depth
				pos <- match(pro, xdepth)
				if(any(is.na(pos)))
					stop(cat("Some depths were not found: ", pro[is.na(pos)], "\n"))
					#stop(cat("Certaine(s) depth(s) demandees n'ont pas put etre trouvees :", pro[is.na(pos)], "\n"))
				pro <- pos
			}
			#Else dans ce cas garde pro comme il etait      
			#validation du group
			if(is(dim1, "missing")) dim1 <- 1:length(xgroupe)
			tmp <- 1:length(xgroupe)
			names(tmp) <- names(xgroupe)
			lg <- tmp[dim1]
			#Nouveau Group a construire
			if(any(is.na(lg))) stop("Invalid extraction: one of the lement was not found")
				#stop("Extraction invalide, un des elements n'a pas ete trouve")
			#extraction
			lg <- xgroupe[dim1]
			#Nouveau Group a construire
			#Regeneration de la liste de proband
			gind <- unlist(xgrindex, use.names = F)
			ggrou <- unlist(xgroupe, use.names = F)
			pro <- ggrou[match(1:(dim(x)[1]), gind)]
			if(any(is.na(pro)))
				stop("Can not use Drop=T for this particular object")
				#stop("Vous ne pouvez pas utiliser Drop=T pour cet objet en particulier")
			#Nouvelle matrice reduite
			class(lg) <- "GLgroup"
			#Creation de l'objet
			GLFGroup(getMethod("[", "matrix")(x,  , pro), lg, xdepth[pro], pro)
		}
		else	stop("If Drop=T, you can only use one or two subscript (depth and optionally groups")
			#stop("Si Drop=T, vous ne pouvez utiliser qu'un ou deux subscript (depth et optionnellement les groups)")
	}
}

GLPrivExtFSINGLE = function(x, ..., drop)
{
	#Ressemble beaucoup a un vecteur #Depend des parametres a drop
	#Extrait le seul parametre 
	#GLOverNumber1 <- function(pro,...,named,abs)
	#list(pro=pro,param=p,named=n,abs=ab) 
	l2 <- GLOverNumber1(...)
	pro <- l2$pro
	if(missing(drop))
		drop <- T
	#Valeur par defaults
	if(l2$param <= 1) {
		#Pour un seul parametre = Sous-ensemble d'element           
		if(drop) {
			#Si 'drop=T', alors retourne la matrice de resultat
			if(is(pro, "GLnothing")) class(pro) <- "missing"
			#Genere le tableau de F moyen par group
			#m <- GLapplyF(x,mean,named=l2$named) #Modifier si sa tartine
			#La suite est equivalent a la ligne ci-haut
			xgrindex <- x@grindex
			xgroupe <- x@group
			x <- as(x, "array")
			ng <- length(xgroupe)
			#nombre de group
			m <- double(ng)
			#Calcul de la valeur (probablement moyen de faire plus efficace)    
			for(a in 1:ng)
				m[a] <- mean(x[xgrindex[[a]]])
			if(l2$named)
				names(m) <- names(xgroupe)
			#Fin construction de m
			getMethod("[", "matrix")(m, pro, drop = T)
		}
		else {
			#Si 'drop = F', alors retourne un objet 'GLmultiPhiGroupSingle' modifie
			#declassement
			xgrindex <- x@grindex
			xgroupe <- x@group
			x <- unclass(x)
			if(l2$param <= 1) {
				if(is(pro, "GLnothing"))
					pro <- 1:length(xgroupe)
				#validation
				tmp <- 1:length(xgroupe)
				names(tmp) <- names(xgroupe)
				lg <- tmp[pro]
				#Nouveau Group a construire
				if(any(is.na(lg))) stop("Invalid extraction: one of the element was not found")
					#stop("Extraction invalide: un des elements n'a pas ete trouve")
				#extraction
				lg <- xgroupe[pro]
				#Nouveau Groupe a construire
				#Regeneration de la liste de proband
				gind <- unlist(xgrindex, use.names = F)
				ggrou <- unlist(xgroupe, use.names = F)
				pro <- ggrou[match(1:(dim(x)[1]), gind)]
				if(any(is.na(pro)))
					stop("you can not use Drop=T for this particular object")
					#stop("Vous ne pouvez pas utiliser 'Drop=T' pour cet objet en particulier")
				#Nouvelle matrice reduite
				class(lg) <- "GLgroup"
				#Creation de l'objet
				GLFGroup(x, lg, proband = pro)
			}
			else	stop("if Drop=T, you can only use one subscript (for the wanted groups)")
				#stop("Si 'Drop=T', vous ne pouvez utiliser qu'un seul subscript (pour les groups voulus)")
		}
	}
	else	stop("You can only use one subscript (as for a vector)")
		#stop("Vous ne pouvez utiliser qu'un subscript (exactement comme un vecteur)")
}

GLPrivExtPHI = function(x, ..., drop)
{
	#Ressemble beaucoup a une GLmultiMatrix....
	#Sauf que 'drop' n'est pas la par defaut... inversement a une matrice....
	#extrait pour certaine depth
	l2 <- GLOverGroup3(...)
	pro <- l2$pro
	dim1 <- l2$dim1
	dim2 <- l2$dim2
	if(missing(drop))
		drop <- T
	#Valeur par defaults
	if(is(dim1, "GLnothing")) class(dim1) <- "missing"
	if(is(dim2, "GLnothing"))
		class(dim2) <- "missing"
	#declassement
	xgrindex <- x@grindex
	xgroupe <- x@group
	xdepth <- x@depth
	x <- unclass(x)
	if(drop) {
		if(l2$param == 0 || l2$param == 1 || l2$param == 3) {
			#Si 'drop=T', alors retourne un GLMultiMatrix 
			#Genere le tableau de phi moyen pour chaque depth            
			m <- array(0, c(length(xgroupe), length(xgroupe), length(xdepth)))
			for(p in 1:length(xdepth)) {
				#Matrice pour une depth donne
				m[,  , p] <- GLapplyGroup(x[,  , p], xgrindex, gen.phiMean, check = 0, named = F)
			}
			if(l2$named)
				dimnames(m) <- list(names(xgroupe), names(xgroupe), NULL)
			else dimnames(m) <- NULL
			#Cas particulier depth
			if(is(pro, "GLnothing")) class(pro) <- "missing" else if(is.numeric(pro) && !l2$abs) {
				#Si c'est pas valeur absolue
				#Recherche comme prevu dans la liste de depth
				pos <- match(pro, xdepth)
				if(any(is.na(pos)))
					stop(cat("Some depths were not found: ", pro[is.na(pos)], "\n"))
					#stop(cat("Certaine(s) depth(s) demandees n'ont pas pu etre trouvees :", pro[is.na(pos)], "\n"))
				pro <- pos
			}
			if(is(dim1, "GLnothing"))
				class(dim1) <- "missing"
			if(is(dim2, "GLnothing"))
				class(dim2) <- "missing"
			#Peut-importe 'dim1' et 'dim2' ca passe a l'operateur 
			GLmulti(getMethod("[", "array")(m, dim1, dim2, pro, drop = drop), xdepth[pro])
		}
		else	stop("You can only use one or three subscript (depth alone or depth, dim1, dim2)")
			#stop("Vous ne pouvez utiliser qu'un ou trois subscript (depth seul ou depth,dim1,dim2)")
	}
	else {
		#Si 'drop = F', alors retourne un objet GLmultiPhiGroup(Single) modifie
		if(l2$param >= 0 && l2$param <= 2) {
			#Cas particulier depth
			if(is(pro, "GLnothing")) class(pro) <- "missing" else if(is.numeric(pro) && !l2$abs) {
				#Recherche comme prevu dans la liste de depth
				pos <- match(pro, xdepth)
				if(any(is.na(pos)))
					stop(cat("Some depth were not found: ", pro[is.na(pos)], "\n"))
					#stop(cat("Certaine(s) depth(s) demandees n'ont pas put etre trouvees :", pro[is.na(pos)], "\n"))
				pro <- pos
			}
			#Else dans ce cas garde pro comme il etait      
			#validation du group
			if(is(dim1, "missing")) dim1 <- 1:length(xgroupe)
			tmp <- 1:length(xgroupe)
			names(tmp) <- names(xgroupe)
			lg <- tmp[dim1]
			#Nouveau Group a construire
			if(any(is.na(lg))) stop("Invalid extraction: one of the element was not found")
				#stop("Extraction invalide, un des elements n'a pas ete trouve")
			#extraction
			lg <- xgroupe[dim1]
			#Nouveau Group a construire
			#Regeneration de la liste de proband
			gind <- unlist(xgrindex, use.names = F)
			ggrou <- unlist(xgroupe, use.names = F)
			pro <- ggrou[match(1:(dim(x)[1]), gind)]
			if(any(is.na(pro)))
				stop("You can not use Drop=T for this particular object")
				#stop("Vous ne pouvez pas utiliser Drop=T pour cet objet en particulier")
			#Nouvelle matrice reduite
			class(lg) <- "GLgroup"
			#Creation de l'objet
			GLPhiGroup(getMethod("[", "array")(x,  ,  , pro), lg, xdepth[pro], pro)
		}
		else	stop("if Drop=T, you can only use one or two subscript (depth and optionally groups)")
		 	#stop("Si Drop=T, vous ne pouvez utiliser qu'un ou deux subscript (depth et optionnellement les groups)")
	}
}

GLPrivExtPHISINGLE = function(x, ..., drop)
{
	#Ressemble beaucoup a une matrice...
	#Depend des parametres a drop
	#extrait pour certaine depth
	l2 <- GLOverGroup2(...)
	dim1 <- l2$dim1
	dim2 <- l2$dim2
	if(missing(drop))
		drop <- T
	#fin if un parametre
	#Valeur par defaut
	if(l2$param < 3) {
		#Pour un seul parametre = Sous-ensemble d'element
		if(is(dim1, "GLnothing")) class(dim1) <- "missing"
		if(is(dim2, "GLnothing"))
			class(dim2) <- "missing"
		#declassement
		xgrindex <- x@grindex
		xgroupe <- x@group
		x <- unclass(x)
		if(drop) {
			#Si drop=T, alors retourne la matrice de resultat
			#Genere le tableau de phi moyen
			m <- GLapplyGroup(x, xgrindex, gen.phiMean, check = 0, named = F)
			if(l2$named)
				dimnames(m) <- list(names(xgroupe), names(xgroupe))
			else dimnames(m) <- NULL
			#Peut-importe dim1 et dim2 ca se passe a l'operateur
			#S'il n'y a qu'un 'param', alors ce comporte comme un vecteur...
			if(l2$param == 1) getMethod("[", "matrix")(m, dim1, drop = drop) else getMethod("[", "matrix")(m, dim1, dim2,
					drop = drop)
		}
		else {
			#Si drop = F ,alors retourne un objet 'GLmultiPhiGroupSingle' modifie
			if(l2$param <= 1) {
				if(is(dim1, "missing"))
					dim1 <- 1:length(xgroupe)
				#validation
				tmp <- 1:length(xgroupe)
				names(tmp) <- names(xgroupe)
				lg <- tmp[dim1]
				#Nouveau Group a construire
				if(any(is.na(lg))) stop("Invalid extraction: one of the element was not found (Error in the 'GLPrivExtPHISINGLE' function)")
					#stop("Extraction invalide: un des elements n'a pas ete trouve (Erreur dans la fonction 'GLPrivExtPHISINGLE')")
				#extraction
				lg <- xgroupe[dim1]
				#Nouveau Group a construire
				#Regeneration de la liste de proband
				gind <- unlist(xgrindex, use.names = F)
				ggrou <- unlist(xgroupe, use.names = F)
				pro <- ggrou[match(1:(dim(x)[1]), gind)]
				if(any(is.na(pro)))
					stop("You can not use Drop=T for this particular object (Error in the 'GLPrivExtPHISINGLE' function)")
					#stop("Vous ne pouvez pas utiliser 'Drop=T' pour cet objet en particulier (Erreur dans la fonction 'GLPrivExtPHISINGLE')")
				#Nouvelle matrice reduite
				class(lg) <- "GLgroup"
				#Creation de l'objet
				GLPhiGroup(x, lg, proband = pro)
			}
			else stop("If Drop=T, you can only use one subscript (for wanted groups)")
				#stop("Si 'Drop=T', vous ne pouvez utiliser qu'un seul subscript (pour les groups voulus)")
		}
	}
	else stop("You can only use one or two subscript (just like a matrix)")
		#stop("Vous ne pouvez utiliser qu'un ou deux subscript (exactement comme une matrice)")
}

GLPrivOcc = function(gen, pro = 0, ancestors = 0)
{
	frequences.anc <- rep(0, length(ancestors))
	repeat {
		if(length(pro) == 0)
			break
		tmptable <- table(pro)
		tmpind <- match(as.integer(names(tmptable)), ancestors)
		tmptable <- tmptable[!is.na(tmpind)]
		tmpind <- tmpind[!is.na(tmpind)]
		frequences.anc[tmpind] <- frequences.anc[tmpind] + tmptable
		tmpind <- match(pro, gen$ind)
		mothers <- gen$mother[tmpind]
		fathers <- gen$father[tmpind]
		pro <- (c(mothers[mothers != 0], fathers[fathers != 0]))
	}
	return(frequences.anc)
}

gen.formatage... = function(...)
{
	parametres = list(...)
	retour = list()
	#Sauvegarder les noms de parametres de l'utilisateur
	#Il doit se conformer a des noms defini par le programmeur
	#ex : mother, father, sex
	#Formater les donnees
	#valeur de retour
	#une liste contenant un vecteur et une liste
	#le vecteur contiendra le nom des parametres valide seulement
	#la liste contiendra, dans le meme ordre que le vecteur, le contenu de chacun des parametres
	if(length(parametres) > 0) {
		nomParametreUtilisateur = names(parametres)
		retour[[2]] = list()
		posMere = match("mother", nomParametreUtilisateur)
		if(!is.na(posMere)) {
			retour[[1]] = c(retour[[1]], nomParametreUtilisateur[
				posMere])
			retour[[2]][[(length(retour[[2]]) + 1)]] = parametres[[
				posMere]]
		}
		posPere = match("father", nomParametreUtilisateur)
		if(!is.na(posPere)) {
			retour[[1]] = c(retour[[1]], nomParametreUtilisateur[
				posPere])
			retour[[2]][[(length(retour[[2]]) + 1)]] = parametres[[
				posPere]]
		}
		posSex = match("sex", nomParametreUtilisateur)
		if(!is.na(posSex)) {
			retour[[1]] = c(retour[[1]], nomParametreUtilisateur[
				posSex])
			retour[[2]][[(length(retour[[2]]) + 1)]] = parametres[[
				posSex]]
		}
	}
	else return(retour)
	return(retour)
}
