# 1 - gen.gc				-> garde article
# 2 - gen.gcplus			-> garde article
# 3 - gen.completeness		-> garde article
# 4 - gen.completenessVar	-> garde article
# 5 - gen.branching			-> garde article
# 6 - gen.children			-> garde article
# 7 - gen.meangendepth		-> garde article
# 8 - gen.entropyMeanVar		-> garde article
# 9 - gen.f				-> garde article
# 11 - gen.fmean			-> garde article
# 12 - gen.founder			-> garde article
# 13 - gen.half.founder		-> garde article
# 14 - gen.sibship			-> garde article
# 16 - gen.genealogy		-> garde article
# 17 - gen.lineages			-> garde article
# 17-1 gen.genout			-> garde article
# 19 - gen.implex			-> garde article
# 20 - gen.implexVar		-> garde article
# 21 - gen.max				-> garde article
# 23 - gen.min				-> garde article
# 24 - gen.mean			-> garde article
# 25 - gen.nochildren		-> garde article
# 26 - gen.nowomen			-> garde article
# 27 - gen.nomen			-> garde article
# 28 - gen.noind			-> garde article
# 32 - gen.occ				-> garde article
# 33 - gen.parent			-> garde article
# 34 - gen.phi				-> garde article
# 35 - gen.phiOver			-> garde article
# 37 - gen.phiMean			-> garde article
# 41 - gen.depth			-> garde article
# 42 - gen.pro				-> garde article
# 43 - gen.rec				-> garde article
# 44 - gen.meangendepthVar		-> garde article


#gen.gc = function(gen, pro = 0, ancestors = 0, typeCG = "IND", named = T, check = 1)
#{
#	if(length(check) != 1)
#		stop("Invalid 'check' parameter: choices are 0 or 1")
#
#	retour = gen.detectionErreur(gen = gen, pro = pro, ancestors = ancestors, print.it = print.it, named = named, typeCG = typeCG,
#							    check = c(3, 5, 11, 34, 18, 10))
#	if(retour$erreur == T)	stop(retour$messageErreur)
#	gen = retour$gen
#	pro = retour$pro
#	ancestors = retour$ancestors
#	typeCG = retour$typeCG
##	print.it = retour$print.it
#	named = retour$named
#
#	if(typeCG == "IND") {
#		if(is(pro, "GLgroup")) {
#			CG = GLPrivCG(gen = gen, pro = as.numeric(unlist(pro)), ancestors = ancestors, print.it = FALSE, named = named)
#			return(GLPrivCGgroup(CG, grppro = pro))
#		}
#		else return(GLPrivCG(gen = gen, pro = pro, ancestors = ancestors, print.it = FALSE, named = named))
#	}
#	else {
#		if(is(pro, "GLgroup")) {
#			CG = GLPrivCG(gen = gen, pro = as.numeric(unlist(pro)), ancestors = ancestors, print.it = FALSE, named = named)
#			CG = GLPrivCGgroup(CG, grppro = pro)
#			if(typeCG == "MEAN")
#				return(GLPrivCGmoyen(CG = CG, named = named))
#			if(typeCG == "CUMUL")
#				stop("CUMUL is not available per group")
#			if(typeCG == "TOTAL")
#				return(GLPrivCGtotal(CG = CG, named = named))
#			if(typeCG == "PRODUCT")
#				return(GLPrivCGproduit(CG = CG, named = named))
#		}
#		else {
#			CG = GLPrivCG(gen = gen, pro = as.numeric(unlist(pro)), ancestors = ancestors, print.it = FALSE, named = 
#				named)
#			if(typeCG == "MEAN")
#				return(GLPrivCGmoyen(CG = CG, named = named))
#			if(typeCG == "CUMUL")
#				return(GLPrivCGcumul(CG = CG, named = named))
#			if(typeCG == "TOTAL")
#				return(GLPrivCGtotal(CG = CG, named = named))
#			if(typeCG == "PRODUCT")
#				return(GLPrivCGproduit(CG = CG, named = named))
#		}
#	}
#}


gen.gc = function(gen, pro = 0, ancestors = 0, vctProb = c(0.5, 0.5, 0.5, 0.5), typeCG = "IND") #, check = 1)#named = T, 
{
	#if(length(check) != 1)	stop("Invalid 'check' parameter: choices are 0 or 1")

	retour = gen.detectionErreur(gen = gen, pro = pro, ancestors = ancestors, print.it = FALSE, named = T, typeCG = typeCG, check = c(3, 5, 11, 34, 18, 10))
	if(retour$erreur == T)	stop(retour$messageErreur)
	gen = retour$gen
	pro = retour$pro
	ancestors = retour$ancestors
	typeCG = retour$typeCG
#	print.it = retour$print.it
	named = retour$named

	if(typeCG == "IND") {
		if(is(pro, "GLgroup")) {
			CG = GLPrivCGPLUS(gen = gen, pro = as.numeric(unlist(pro)), ancestors = ancestors, vctProb = vctProb, print.it = FALSE, named = named)
			return(GLPrivCGgroup(CG, grppro = pro))
		}
		else return(GLPrivCGPLUS(gen = gen, pro = pro, ancestors = ancestors, vctProb, print.it = FALSE, named = named))
	}
	else {
		if(is(pro, "GLgroup")) {
			CG = GLPrivCGPLUS(gen = gen, pro = as.numeric(unlist(pro)), ancestors = ancestors, vctProb = vctProb, print.it = FALSE, named = named)
			CG = GLPrivCGgroup(CG, grppro = pro)
			if(typeCG == "MEAN")
				return(GLPrivCGmoyen(CG = CG, named = named))
			if(typeCG == "CUMUL")
				stop("CUMUL is not available per group")
			if(typeCG == "TOTAL")
				return(GLPrivCGtotal(CG = CG, named = named))
			if(typeCG == "PRODUCT")
				return(GLPrivCGproduit(CG = CG, named = named))
		}
		else {
			CG = GLPrivCGPLUS(gen = gen, pro = as.numeric(unlist(pro)), ancestors = ancestors, vctProb = vctProb, print.it = FALSE, named = named)
			if(typeCG == "MEAN")
				return(GLPrivCGmoyen(CG = CG, named = named))
			if(typeCG == "CUMUL")
				return(GLPrivCGcumul(CG = CG, named = named))
			if(typeCG == "TOTAL")
				return(GLPrivCGtotal(CG = CG, named = named))
			if(typeCG == "PRODUCT")
				return(GLPrivCGproduit(CG = CG, named = named))
		}
	}
}


gen.completeness = function(gen, pro = 0, genNo = -1, type = "MEAN", ...)#, check = 1)#named = T, 
{
	#Validation des parametres
	#if(length(check) != 1) stop("Invalid 'check' parameter: choices are 0 or 1")
	if( type != "IND" )# | type != "MOYSUJETS"
		if(is(gen, "vector"))
			if(length(list(...)) != 2)
				stop("Invalid '...' parameter : 'father' and 'mother' parameter names are obligatory")
	#if(check == 1) {
		retour <- gen.detectionErreur(gen = gen, pro = pro, genNo = genNo, typecomp = type, named = T, check = c(1, 5, 16, 171, 10))
		if(retour$erreur == T)
			stop(retour$messageErreur)
		gen <- retour$gen
		pro <- retour$pro
		genNo <- retour$genNo
		type <- retour$typecomp
		named <- retour$named
	#}
	#Calcule de la completude par sujet
	if(type == "IND"){ # | type == "MOYSUJETS") {
		tableau = sapply(pro, function(x, gen, genNo, named)
		{
			GLPriv.completeness3V(gen$ind, gen$father, gen$mother, pro = x, genNo = genNo, named = named)
		}
		, gen = gen, genNo = genNo, named = named)
		if(is.null(dim(tableau))) tableau <- t(as.matrix(tableau))
#		if(type == "MOYSUJETS")
#			tableau <- data.frame(apply(tableau, 1, mean))
		#Fait la moyenne 
		#if(named == T)
		  if(type == "IND")
			dimnames(tableau)[[2]] <- as.character(paste("Ind", as.character(pro)))
		dimnames(tableau)[[1]] <- as.character(genNo)
		#Rajout du numero de generation en lignes
		return(data.frame(tableau))
	}
	else if(type == "MEAN") {
		#Si c'est MEAN, calcul de la completude avec tous les sujets a la fois
		return(GLPriv.completeness3V(gen$ind, gen$father, gen$mother, pro = pro, genNo = genNo, named = named))
	}
}


gen.completenessVar = function(gen, pro = 0, genNo = -1, ...) #, check = 1, ...)#named = T, 
{
	#Validations des parametres 
	#if(length(check) != 1) stop("Invalid 'check' parameter: choices are 0 or 1")
#	if(bcorrFactor == T)
#		if(sum(N) == 0) #Le facteur de correction doit avoir une valeur numerique N taille de la population
#			stop("Correction factor must have a numerical population size value N")
	if(is(gen, "vector"))
		if(length(list(...)) != 2)
			stop("Invalid '...' parameter : 'father' and 'mother' parameter names are obligatory")

	retour = gen.detectionErreur(gen = gen, pro = pro, genNo = genNo, named = T, check = c(1, 5, 16, 10))
	if(retour$erreur == T)	stop(retour$messageErreur)
	gen = retour$gen
	pro = retour$pro
	genNo = retour$genNo
	named = retour$named

	#Selon le type de donnees, le facteur de correction sera modifie en consequence
#	if(typeCorpus == "ECH") corrFactor = length(pro)/(length(pro) - 1) else if(typeCorpus == "POP")
#		corrFactor = 1
#	if(corrFactor == T)
#		corrFactor = (corrFactor * (N - length(pro)))/N

	corrFactor = 1
	#Calcule la variance de l'indice de completude
	tab = sapply(pro, function(x, gen, genNo, named)
		GLPriv.completeness3V(gen$ind, gen$father, gen$mother, pro = x, genNo = genNo, named = named),
		gen = gen, genNo = genNo, named = named)
	
	if(is.null(dim(tab))) tab <- t(as.matrix(tab))
	tab = data.frame(apply(tab, 1, var) * corrFactor)
	dimnames(tab)[[1]] <- as.character(genNo)
	dimnames(tab)[[2]] <- "completeness.var"
	return(tab)
}


gen.branching = function(gen, pro = 0, ancestors = gen.founder(gen), bflag = 0)#, check = 1)
{
	if(sum(as.numeric(pro)) == 0)
		pro = gen.pro(gen)
	if(bflag == 0) {
		pro.segment = gen.pro(gen)
		ancestors = gen.founder(gen.branching(gen, pro.segment, ancestors, bflag = 1))
	}
	#if(length(check) != 1) stop("Invalid 'check' parameter: choices are 0 or 1")
		#stop("Param\350tre 'check' invalide: les choix disponibles sont 0 et 1")
	#if(check == 1) {
		retour = gen.detectionErreur(gen = gen, pro = pro, ancestors = ancestors, check = c(3, 36, 37))
		if(retour$erreur)
			stop(retour$messageErreur)
		gen = retour$gen
		pro = retour$pro
		ancestors = retour$ancestors
	
	#}
	#Structure necessaire pour emmagasiner le resultat la fonction de la dll
	#print(paste("taille alloue:",length(gen@.Data)))
	tmpgen <- integer(length(gen@.Data))
	tmpNelem <- integer(1)

	#print(".C(SPLUSebranche,.. commence")
	.Call("SPLUSebranche", gen@.Data, pro, length(pro), ancestors, length(ancestors), tmpgen, tmpNelem, specialsok = T)
	#print(".C(SPLUSebranche,.. fait:")
	#print(paste(length(tmpgen),tmpgen[1],tmpgen[2],tmpgen[3] ))
	#print(paste(length(gen@.Data),gen@.Data[1],gen@.Data[2],gen@.Data[3]))
	length(tmpgen) <- tmpNelem
	tmpNelem <- length(tmpgen)
	#print(length(tmpgen))
	ebranche = new("GLgen", .Data = tmpgen, Date = date())
	#print("1")
	ebranche.asc = gen.genout(ebranche)
	sexeAbsent=F
	if(length(setdiff(unique(ebranche.asc[,"sex"]), c(1,2,"H","F")))>0) 
	{
	  diff = setdiff(unique(ebranche.asc[,"sex"]), c(1,2,"H","F"))
	  ebranche.asc=data.frame(ind=ebranche.asc$ind,father=ebranche.asc$father,mother=ebranche.asc$mother) #*****
	  sexeAbsent=T
	  #warning(paste("la colonne \"sexe\" contient des valeurs non valide:",diff,"\n  Elle ne sera pas consideree pour le reste des calculs."))
	  warning(paste("The \"sex\" column contains invalid values:",diff,
					"\nThe column won't be considered for further calculations."))
	}
	#print("2")
	#print(ebranche.asc[1,])
	pro.ebranche = gen.pro(ebranche)
	#print("3")
	pro.enTrop = setdiff(pro.ebranche, pro)
	#print(paste(length(pro.ebranche),length(pro)))
	#print(pro.enTrop)
	if(sum(as.numeric(pro.enTrop)) != 0) {
		ebranche.asc = ebranche.asc[(!(ebranche.asc$ind %in% pro.enTrop)),  ]
		#ebranche.asc=data.frame(ind=ebranche.asc$ind,father=ebranche.asc$father,mother=ebranche.asc$mother) #*****
		ebranche = gen.genealogy(ebranche.asc)
	#print(ebranche.asc)
		pro.ebranche = gen.pro(ebranche)
	}
	#print("4")
	fond.ebranche = gen.founder(ebranche)
	#print("5")
	pro.quiSontFond = pro.ebranche[pro.ebranche %in% fond.ebranche]
	#print(paste("6", dim(ebranche.asc)))
	ebranche.asc = ebranche.asc[(!(ebranche.asc$ind %in% pro.quiSontFond)),  ]
	#print(paste("7", dim(ebranche.asc)))
	#ebranche.asc=data.frame(ind=ebranche.asc$ind,father=ebranche.asc$father,mother=ebranche.asc$mother)#*****
	if(dim(ebranche.asc)[1]==0) stop("No branching possible, all probands are founders.")
	else gen = gen.genealogy(ebranche.asc)
	#print("8")
	gen.validationAsc(gen)
	#print("9")
	return(gen)
}


gen.children = function(gen, individuals, ...)#, check = 1)
{
	#if(length(check) != 1) stop("Invalid 'check' parameter: choices are 0 or 1")
		#stop("Param\350tre 'check' invalide: les choix disponibles sont 0 et 1")
	if(is(gen, "vector"))
		if(length(list(...)) != 2)
			stop("Invalid '...' parameter : 'father' and 'mother' parameter names are obligatory")
			#stop("Param\350tre '...' invalide : indication du nom des param\350tres 'pere' et 'mere' est obligatoire")
	#if(check == 1) {
		retour = gen.detectionErreur(gen = gen, individuals = individuals, check = c(1, 13), ...)
		if(retour$erreur == T)
			stop(retour$messageErreur)
		gen = retour$gen
		individuals = retour$individuals
	#}
	PositionEnfantDesMeres <- match(gen$mother, individuals)
	PositionEnfantDesPeres <- match(gen$father, individuals)
	EnfantDesMere <- gen$ind[(1:length(PositionEnfantDesMeres))[!is.na(PositionEnfantDesMeres)]]
	EnfantDesPere <- gen$ind[(1:length(PositionEnfantDesPeres))[!is.na(PositionEnfantDesPeres)]]
	Enfants <- unique(c(EnfantDesMere, EnfantDesPere))
	return(Enfants)
}


gen.meangendepth = function(gen, pro = 0, type = "MEAN", ...)#, check = 1)#named = T, 
{
	#Validations des parametres
	#if(length(check) != 1) stop("Invalid 'check' parameter: choices are 0 or 1")
		#stop("Param\350tre 'check' invalide: les choix disponibles sont 0 et 1")
	if(is(gen, "vector"))
		if(length(list(...)) != 2)
			stop("Invalid '...' parameter : 'father' and 'mother' parameter names are obligatory")
			#stop("Param\350tre '...' invalide : indication du nom des param\350tres 'pere' et 'mere' est obligatoire")
	#if(check == 1) {
		retour <- gen.detectionErreur(gen = gen, pro = pro, typecomp = type, check = c(1, 5, 17))
		if(retour$erreur == T)
			stop(retour$messageErreur)
		gen <- retour$gen
		pro <- retour$pro
		type <- retour$typecomp
	#}
	if(type == "IND") {# | type == "MOYSUJETS") {
		tableau <- sapply(pro, function(x, gen)
		GLPriv.entropie3V(gen$ind, gen$father, gen$mother, pro = x), gen = gen)
		tableau <- data.frame(tableau)
#		if(type == "MOYSUJETS") {
#			tableau = data.frame(apply(tableau, 2, mean))
#			dimnames(tableau)[[1]] <- "Mean.Exp.Gen.Depth"
#		}
		#if(named == T)
			if(type == "IND")
				dimnames(tableau)[[1]] <- as.character(paste("Ind", as.character(pro)))
		dimnames(tableau)[[2]] <- "Exp.Gen.Depth"
		return(tableau)
	}
	else if(type == "MEAN")
		return(GLPriv.entropie3V(gen$ind, gen$father, gen$mother, pro = pro))
}


#gen.entropyMeanVar = function(gen, pro = 0, check = 1, ...) #typeCorpus = "ECH", bfacteurCorr = F, N = NULL,
#{
#	#Validations des parametres
#	if(length(check) != 1) stop("Invalid 'check' parameter: choices are 0 or 1")
##	if(bfacteurCorr == T)
#		if(sum(N) == 0)
##			stop("Correction factor must have a numerical population size value N")
#	if(is(gen, "vector"))
#		if(length(list(...)) != 2)
#			stop("Invalid '...' parameter : 'father' and 'mother' parameter names are obligatory")
#
#	retour = gen.detectionErreur(gen = gen, pro = pro, check = c(1, 5))
#	if(retour$erreur == T)	stop(retour$messageErreur)
#	gen = retour$gen
#	pro = retour$pro
#
#	tableau = sapply(pro, function(x, gen)
#	GLPriv.entropie3V(gen$ind, gen$father, gen$mother, pro = x), gen = gen)
##	if(typeCorpus == "ECH")
##		facteurCorr = length(pro)/(length(pro) - 1)
##	else if(typeCorpus == "POP")
##		facteurCorr = 1
##	if(bfacteurCorr == T)
##		facteurCorr = (facteurCorr * (N - length(pro)))/N
#	
#	facteurCorr = 1
#	tableau = data.frame(tableau)
#	tableau = data.frame(apply(tableau, 2, var) * facteurCorr)
#	dimnames(tableau)[[1]] <- "Mean.Exp.Gen.Depth.Var"
#	dimnames(tableau)[[2]] <- "Exp.Gen.Depth"
#	return(tableau)
#}


#gen.f = function(gen, pro = 0, nbgenerations = 0, named = T, check = 1)
#{
#	if(length(check) != 1)	stop("Invalid 'check' parameter: choices are 0 or 1")
#	retour = gen.detectionErreur(gen = gen, pro = pro, nbgenerations = nbgenerations, print.it = FALSE, named = named, check = c(3, 5, 19, 18, 10))
#	if(retour$erreur == T)	stop(retour$messageErreur)
#	gen = retour$gen
#	pro = retour$pro
#	nbgenerations = retour$nbgenerations
##	print.it = retour$print.it
#	named = retour$named
#
#	#Structure necessaire pour emmagasiner le resultat la fonction de la dll
#	tmp <- double(length(pro))
#
#	#Call de la fonction en C
#	.Call("SPLUSF", gen@.Data, pro, length(pro), nbgenerations, tmp, FALSE, specialsok = T)
#	names(tmp) <- pro
##	if(print.it) {
##		base <- c(deparse(substitute(gen)), deparse(substitute(pro)), nbgenerations)
##		header.txt <- paste("\n\t***   Calls : gen.F (", base[1], ",", base[2], ",", base[3], ")  ***\n\n")
##		cat(header.txt)
##	}
#	return(invisible(tmp))
#}


#gen.fmean = function(vectF, named = T, check = 1)
#{
#	if(length(check) != 1)
#		stop("Invalid 'check' parameter: choices are 0 or 1")
#		#stop("Param\350tre 'check' invalide: les choix disponibles sont 0 et 1")
#	#if(check == 1) {
#		retour = gen.detectionErreur(vectF = vectF, named = named, check = c(33, 10))
#		if(retour$erreur == T)
#			stop(retour$messageErreur)
#		vectF = retour$vectF
#		named = retour$named
#	#}
#	#Test pour accelerer la procedure
#	return(GLapplyF(vectF, mean, named = named))
#}


gen.founder = function(gen, ...)#, check = 1)
{
	#if(length(check) != 1) stop("Invalid 'check' parameter: choices are 0 or 1")
		#stop("Param\350tre 'check' invalide: les choix disponibles sont 0 et 1")
	if(is(gen, "vector"))
		if(length(list(...)) != 2)
			stop("Invalid '...' parameter : 'father' and 'mother' parameter names are obligatory")
			#stop("Param\350tre '...' invalide : indication du nom des param\350tres 'pere' et 'mere' est obligatoire")
	#if(check == 1) {
		retour = gen.detectionErreur(gen = gen, ..., check = 1)
		if(retour$erreur == T)
			return(retour$messageErreur)
		gen = retour$gen
	#}
	return(gen$ind[gen$father == 0 & gen$mother == 0])
}


gen.half.founder = function(gen, ...)#, check = 1)
{
	#if(length(check) != 1) stop("Invalid 'check' parameter: choices are 0 or 1")
		#stop("Param\350tre 'check' invalide: les choix disponibles sont 0 et 1")
	if(is(gen, "vector"))
		if(length(list(...)) != 2)
			stop("Invalid '...' parameter : 'father' and 'mother' parameter names are obligatory")
			#stop("Param\350tre '...' invalide : indication du nom des param\350tres 'pere' et 'mere' est obligatoire")
	#if(check == 1) {
		retour = gen.detectionErreur(gen = gen, ..., check = 1)
		if(retour$erreur == T)
			return(retour$messageErreur)
		gen = retour$gen
	#}
	return(gen$ind[(gen$father != 0 & gen$mother == 0) | (gen$father == 0 & gen$mother != 0)])
}


gen.sibship = function(gen, individuals, halfSibling = T, ...)#, check = 1)
{
	#if(length(check) != 1) stop("Invalid 'check' parameter: choices are 0 or 1")
		#stop("Param\350tre 'check' invalide: les choix disponibles sont 0 et 1")
	if(is(gen, "vector"))
		if(length(list(...)) != 2)
			stop("Invalid '...' parameter : 'father' and 'mother' parameter names are obligatory")
			#stop("Param\350tre '...' invalide : indication du nom des param\350tres 'pere' et 'mere' est obligatoire")
	#if(check == 1) {
	retour = gen.detectionErreur(gen = gen, individuals = individuals, halfSibling = halfSibling, check = c(1, 13, 14), ...)
	if(retour$erreur == T) stop(retour$messageErreur)
	gen = retour$gen
	individuals = retour$individuals
	halfSibling = retour$halfSibling
	#}
	if(halfSibling == T) {
		PositionProband = match(individuals, gen$ind)
		#Trouve les meres et les peres des probands
		Meres <- gen$mother[PositionProband]
		Peres <- gen$father[PositionProband]
		MaskMere <- Meres != 0
		MaskPere <- Peres != 0
		Meres <- (Meres/MaskMere)[!is.na(Meres/MaskMere)]
		Peres <- (Peres/MaskPere)[!is.na(Peres/MaskPere)]
		#Trouve tous les enfants de ces individuals
		sibshipMo <- gen.children(gen, individuals = Meres)#, check = 0)
		sibshipFa <- gen.children(gen, individuals = Peres)#, check = 0)
		#Vecteur contenant tous les enfants incluant les probands
		sibshipAndProband <- unique(c(sibshipMo, sibshipFa))
		#maintenant on enleve les probands
		temp <- match(sibshipAndProband, individuals)
		sibship <- sibshipAndProband[(1:length(temp))[is.na(temp)]]
		return(sibship)
	}
	else {
		PositionProband = match(individuals, gen$ind)
		#Trouve les meres et les peres des probands
		Meres <- gen$mother[PositionProband]
		Peres <- gen$father[PositionProband]
		MaskMere <- Meres != 0
		MaskPere <- Peres != 0
		Meres <- (Meres/MaskMere)[!is.na(Meres/MaskMere)]
		Peres <- (Peres/MaskPere)[!is.na(Peres/MaskPere)]
		temp1 <- match(gen$mother, Meres)
		temp2 <- match(gen$father, Peres)
		PositionsibshipAndProband <- temp1 * temp2
		#La sibship incluant les probands
		sibshipSameFaMoAndProband <- gen$ind[(1:length(PositionsibshipAndProband))[!is.na(PositionsibshipAndProband)]]
		#maintenant enlevons les probands
		temp <- match(sibshipSameFaMoAndProband, individuals)
		sibship <- sibshipSameFaMoAndProband[(1:length(temp))[is.na(temp)]]
		return(sibship)
	}
}


gen.f = function(gen, pro, depthmin= (gen.depth(gen)-1), depthmax= (gen.depth(gen)-1)) #, check = 1)#named = T, 
{
	#if(length(check) != 1)	stop("Invalid 'check' parameter: choices are 0 or 1")
	if(missing(pro))		pro = gen.pro(gen)

	retour = gen.detectionErreur(gen = gen, pro = pro, depthmin = depthmin, depthmax = depthmax, print.it = FALSE, named = T, 
							check = c(3, 5, 20, 18, 10))
	if(retour$erreur == T) stop(retour$messageErreur)
	gen = retour$gen
	pro = retour$pro
	depthmin = retour$depthmin
	depthmax = retour$depthmax
	named = retour$named

	#Structure necessaire pour emmagasiner le resultat la fonction de la dll
	ecart <- as.integer(depthmax) - as.integer(depthmin) + 1
	tmp <- double(length(pro) * ecart)

	#Call de la fonction en C
	.Call("SPLUSFS", gen@.Data, pro, length(pro), depthmin, depthmax, tmp, FALSE, specialsok = T)
	#Construction de la matrice de retour
	dim(tmp) <- c(length(pro), ecart)
	dimnames(tmp) <- list(pro, NULL)
	tmp = drop(tmp)
	return(invisible(GLmulti(tmp, depth = as.integer(depthmin:depthmax))))
}


gen.genealogy = function(ped, autoComplete=FALSE, ...)#, check = 1)
{
	if(class(ped) != "GLgen") {
	 if(dim(ped)[2]==4 && sum(colnames(ped)==c("X1","X2","X3","X4"))==4) {
	  print("No column names given. Assuming <ind>, <father>, <mother> and <sex>")
	  colnames(ped) <- c("ind", "father", "mother", "sex")
	 }
	 if(sum(c("ind","father","mother","sex") %in% colnames(ped)) < 4){
	  stop(paste(paste(c("ind","father","mother","sex")[grep(F,c("ind","father","mother","sex") %in% colnames(ped))]),
	 		"not in table columns.",collapse=""))
	 }
	 if(autoComplete & !all(is.element(ped[ped[,"father"]!=0,"father"], ped[,"ind"]))) {
	 	pereManquant <- unique(ped[grep(F, is.element(ped[,"father"], ped[,"ind"])),"father"])
	 	pereManquant <- pereManquant[-grep("^0$",pereManquant)]
	 	ajout <- matrix(c(pereManquant, rep(0, (2*length(pereManquant))), rep(1,length(pereManquant))), byrow=F, ncol=4)
	 	colnames(ajout) <- colnames(ped)
	 	ped <- rbind(ped, ajout)
	 }
	 if(autoComplete & !all(is.element(ped[ped[,"mother"]!=0,"mother"], ped[,"ind"]))) {
	 	mereManquante <- unique(ped[grep(F, is.element(ped[,"mother"], ped[,"ind"])),"mother"])
	 	mereManquante <- mereManquante[-grep("^0$",mereManquante)]
	 	ajout <- matrix(c(mereManquante, rep(0, (2*length(mereManquante))), rep(2,length(mereManquante))), byrow=F, ncol=4)
	 	colnames(ajout) <- colnames(ped)
	 	ped <- rbind(ped, ajout)
	 }
	}
	#if(length(check) != 1)	stop("Invalid 'check' parameter: choices are 0 or 1")
	retour = gen.detectionErreur(gen = ped, check = 1, ...)
	if(retour$erreur == T)	stop(retour$messageErreur)
	gen = retour$gen

	tmp2 <- NULL
	if(!is.null(gen$sex)) {
		tmp <- factor(gen$sex, levels = c("H", "h", 1, "F", "f", 2))
		tmp2 <- as(tmp, "integer")
		tmp2[tmp2 == 2 | tmp2 == 3] <- 1
		tmp2[tmp2 == 4 | tmp2 == 5 | tmp2 == 6] <- 2
	}
	n <- .Call("SPLUSCALLCreerObjetGenealogie", gen$ind, gen$father, gen$mother, tmp2)
	
	#Creation de l'objet Genealogie
	return(new("GLgen", .Data = n, Date = date()))
}


gen.lineages = function(ped, pro = 0, maternal = T, ...)#, check = 1
{
	#Creation d'un objet GLgen avec toutes les ascendances
	gen = gen.genealogy(ped, ...) #check = check,
	#Validation des parametres gen et proband

	retour = gen.detectionErreur(gen = gen, pro = pro, check = c(3, 36))
	if(retour$erreur == T)	stop(retour$messageErreur)
	gen = retour$gen
	pro = retour$pro

	#Si des sujets ne sont pas forces, par defaut les individuals n'ayant pas d'enfants sont selectionnes
	if(sum(pro == 0)) data.ind = gen.pro(gen) else data.ind = pro
	#Si c'est des lignees maternelles, les tous les peres sont mis a 0, sinon c'est les meres 
	if(maternal == T) {
		ped$father = rep(0, length(ped$father))
#		output = "M"
	}
	else {
		ped$mother = rep(0, length(ped$mother))
#		output = "F"
	}
	#On cree un objet GLgen avec les meres ou les peres a 0
	genMouP = gen.genealogy(ped, ...) #, check = check
	lig.parent.lst = c(data.ind)
	#Pour toutes les depths, on prend les parents a partir des sujets
	for(i in 1:gen.depth(gen)) {
		data.ind = unlist(gen.parent(genMouP, data.ind))
		lig.parent.lst = c(lig.parent.lst, data.ind)
	}
	#Du resultat, on extrait les individuals de la table d'ascendances qui sont presents
	gen = gen.genealogy(ped[(ped$ind %in% lig.parent.lst),  ], ...) #, check = check
	#Retourne l'objet GLgen de lignees
	return(gen)
}


gen.genout = function(gen, sorted = F)#, check = 1)
{
	#if(length(check) != 1) stop("Invalid 'check' parameter: choices are 0 or 1")
		#stop("Param\350tre 'check' invalide: les choix disponibles sont 0 et 1")
	#if(check == 1) {
		retour = gen.detectionErreur(gen = gen, sorted = sorted, check = c(3, 4))
		if(retour$erreur == T) stop(retour$messageErreur)
		gen = retour$gen
		sorted = retour$sorted
	#}
	#Structure necessaire pour emmagasiner le resultat la fonction de la dll 
	#print(paste(" ? ",gen@.Data[9]))
	taille <- gen.noind(gen)
	v <- list(ind = integer(taille), father = integer(taille), mother = integer(taille), sex = integer(taille))

	#extern "C" void SPLUSOutgen
	#(long* genealogie, long* plRetIndividu,long* plRetPere,long* plRetMere,long* mustsort)
	#param <- list(Data=gen@.Data, ind=v$ind, father=v$father, mother=v$mother, sex=v$sex, sorted=sorted)
	#param = .Call("SPLUSOutgen", param, NAOK = T)
	param = .Call("SPLUSOutgen", gen@.Data, v$ind, v$father, v$mother, v$sex, sorted)
	v <- list(ind = param$ind, father = param$father, mother = param$mother, sex = param$sex)
	#Si le numero du sex (0 ou 1 )des individuals est present, on les change pour "H" ou "F"
	#if(v$sex[1] == -1) v <- v[1:3]
	#else 			v[[4]] <- factor(v[[4]], labels = c("H", "F"))
	return(invisible(data.frame(v)))
}

gen.implex = function(gen, pro = 0, genNo = -1, type = "MEAN", onlyNewAnc = F, ...)#, check = 1 named = T, 
{
	#Validations des parametres 
	#if(length(check) != 1) stop("Invalid 'check' parameter: choices are 0 or 1")
		#stop("Param\350tre 'check' invalide: les choix disponibles sont 0 et 1")
	if(is(gen, "vector"))
		if(length(list(...)) != 2)
			stop("Invalid '...' parameter : 'father' and 'mother' parameter names are obligatory")
			#stop("Param\350tre '...' invalide : indication du nom des param\350tres 'pere' et 'mere' est obligatoire")
	#if(check == 1) {
		retour <- gen.detectionErreur(gen = gen, pro = pro, genNo = genNo, typecomp = type, named = T, check = c(1,	5, 16, 17, 10))
		if(retour$erreur == T)
			stop(retour$messageErreur)
		gen <- retour$gen
		pro <- retour$pro
		genNo <- retour$genNo
		named <- retour$named
		type <- retour$typecomp
	#}
	#Les ancetres se repetent sur plusieurs generations
	#Si on veut les ancetres distincts par generation nouveaux ou pas la fonctionnalite utilisee sera differente
	if(onlyNewAnc == F) fctApp <- GLPriv.implex3V else fctApp <- gen.implex3V
	#Les ancetres ne sont comptes qu'a leur 1ere apparition
	#Selon le type du calcul  
	#Calcule de l'implex par sujet 
	if(type == "IND" | type == "MEAN") {
		tableau = sapply(pro, function(x, gen, genNo, fctApp, named)
		{
			fctApp(gen$ind, gen$father, gen$mother, pro = x, genNo = genNo, named = named)
		}
		, gen = gen, genNo = genNo, fctApp = fctApp, named = named)
		if(is.null(dim(tableau))) tableau <- t(as.matrix(tableau))
		#Selon le resultat, on applique au tableau une operation de moyenne ou pas
		if(type == "MEAN")	tableau = data.frame(apply(tableau, 1, mean))
		#if(named == T)
		#dimnames(tableau)[[2]] <- "implex"
		names(tableau) <- "implex"
		if(type == "IND")	dimnames(tableau)[[2]] <- as.character(paste("Ind", as.character(pro)))
		dimnames(tableau)[[1]] <- as.character(genNo)
		return(data.frame(tableau))
	}
	else if(type == "ALL")
		return(fctApp(gen$ind, gen$father, gen$mother, pro = pro, genNo = genNo, named = named))
}


gen.implexVar = function(gen, pro = 0, onlyNewAnc = F, genNo = -1, ...)# check = 1,named = T, 
{
	#Validation des parametres
	#if(length(check) != 1) stop("Invalid 'check' parameter: choices are 0 or 1")
#	if(bfacteurCorr == T)
#		if(sum(N) == 0)
#			stop("Correction factor must have a numerical population size value N")
	if(is(gen, "vector"))
		if(length(list(...)) != 2)
			stop("Invalid '...' parameter : 'father' and 'mother' parameter names are obligatory")

	retour = gen.detectionErreur(gen = gen, pro = pro, genNo = genNo, named = T, check = c(1, 5, 16, 10))
	if(retour$erreur == T)	stop(retour$messageErreur)
	gen = retour$gen
	pro = retour$pro
	genNo = retour$genNo
	named = retour$named

	#Si on veut les ancetres distincts par generation nouveaux ou pas la fonctionnalite utilisee sera differente
	if(onlyNewAnc == F) fctApp <- GLPriv.implex3V else fctApp <- gen.implex3V
	#Selon le type de donnees, le facteur de correction sera modifie en consequence
#	if(typeCorpus == "ECH") facteurCorr = length(pro)/(length(pro) - 1) else if(typeCorpus == "POP")
#		facteurCorr = 1
#	if(bfacteurCorr == T)
#		facteurCorr = (facteurCorr * (N - length(pro)))/N

	facteurCorr = 1
	tableau = sapply(pro, function(x, gen, fctApp, genNo, named) {
				fctApp(gen$ind, gen$father, gen$mother, pro = x, genNo = genNo, named = named)
			}, gen = gen, fctApp = fctApp, genNo = genNo, named = named)
	if(is.null(dim(tableau))) tableau <- t(as.matrix(tableau))
	tableau = data.frame(apply(tableau, 1, var) * facteurCorr)
	dimnames(tableau)[[1]] <- as.character(genNo)
	dimnames(tableau)[[2]] <- "implex.var"
	return(tableau)
}


#gen.max = function(gen, individuals, named = T, check = 1)
#{
#	#On appel la fonction qui permet d'avoir
#	#le numero de generation de tout les individuals
#	dfData.numgen = gen.generation(gen, as.integer(individuals))
#	dfResult = as.data.frame(as.numeric(names(dfData.numgen))) #named.index.rowcol( dfData.numgen, "numeric")
#	dfResult[, 2] = dfData.numgen
#	dimnames(dfResult)[[2]] <- c("ind", "numgen")
#	return(dfResult)
#}

gen.max = function(gen, individuals)#, check = 1) #, ancestors=0)named = T, 
{
	#if(length(check) != 1)	stop("Invalid 'check' parameter: choices are 0 or 1")

	retour = gen.detectionErreur(gen = gen, individuals = individuals, ancestors = 0, named = T, check = c(3, 13, 10))
	if(retour$erreur == T)	stop(retour$messageErreur)
	gen			= retour$gen
	individuals = retour$individuals
	named		= retour$named

	#Structure necessaire pour emmagasiner le resultat la fonction de la dll
	nPro = length(individuals)
	ret = integer(nPro)
	#extern "C"  void  SPLUSnumeroGen(long* Genealogie, long* lpProband, NProband, retour) 
	.Call("SPLUSnumeroGen", gen@.Data, as.integer(individuals), nPro, ret)
	#if(named)
	names(ret) <- individuals
	return(ret)
}

gen.min = function(gen, individuals)#, check = 1) #, ancestors=0)named = T, 
{
	#if(length(check) != 1) stop("Invalid 'check' parameter: choices are 0 or 1")
		#stop("Param\350tre 'check' invalide: les choix disponibles sont 0 et 1")
	#if(check == 1) {
		retour = gen.detectionErreur(gen = gen, individuals = individuals, ancestors = 0, named = T, check = c(3,13,10))
		if(retour$erreur == T)
			stop(retour$messageErreur)
		gen = retour$gen
		individuals = retour$individuals
		named = retour$named
	#}
	#Structure necessaire pour emmagasiner le resultat la fonction de la dll
	nPro = length(individuals)
	ret = integer(nPro)
	#extern "C"  void  SPLUSnumeroGen(long* Genealogie, long* lpProband, NProband, retour) 
	.Call("SPLUSnumGenMin", gen@.Data, as.integer(individuals), nPro, ret)
	#if(named)
		names(ret) <- individuals
	return(ret)
}

gen.mean = function(gen, individuals)#, check = 1) #, ancestors=0)named = T, 
{
	#if(length(check) != 1) stop("Invalid 'check' parameter: choices are 0 or 1")
		#stop("Param\350tre 'check' invalide: les choix disponibles sont 0 et 1")
	#if(check == 1) {
		retour = gen.detectionErreur(gen = gen, individuals = individuals, ancestors = 0, named = T, check = c(3,13,10))
		if(retour$erreur == T)
			stop(retour$messageErreur)
		gen = retour$gen
		individuals = retour$individuals
		named = retour$named
	#}
	#Structure necessaire pour emmagasiner le resultat la fonction de la dll
	nPro = length(individuals)
	ret = double(nPro)
	#extern "C"  void  SPLUSnumeroGen(long* Genealogie, long* lpProband, NProband, retour) 
	.Call("SPLUSnumGenMoy", gen@.Data, as.integer(individuals), nPro, ret)
	#if(named)
		names(ret) <- individuals
	return(ret)
}

gen.nochildren = function(gen, individuals)#, check = 1)#named = T, 
{
	#if(length(check) != 1) stop("Invalid 'check' parameter: choices are 0 or 1")
		#stop("Param\350tre 'check' invalide: les choix disponibles sont 0 et 1")
	#if(check == 1) {
		retour = gen.detectionErreur(gen = gen, individuals = individuals, named = T, check = c(3, 13, 10))
		if(retour$erreur == T)
			stop(retour$messageErreur)
		gen = retour$gen
		individuals = retour$individuals
		named = retour$named
	#}
	#Structure necessaire pour emmagasiner le resultat la fonction de la dll		
	ret <- integer(length(individuals))
	#extern "C" void SPLUSChild(long* Genealogie, long* plProband,long* lNProband, long* retour)
	.Call("SPLUSChild", gen@.Data, individuals, length(individuals), ret, specialsok = T)
	#if(named)
		names(ret) <- individuals
	return(ret)
}

gen.nowomen = function(gen)#, check = 1)
{
	#if(length(check) != 1)stop("Invalid 'check' parameter: choices are 0 or 1")
		#stop("Param\350tre 'check' invalide: les choix disponibles sont 0 et 1")
	#if(check == 1) {
		retour = gen.detectionErreur(gen = gen, check = 3)
		if(retour$erreur == T)
			stop(retour$messageErreur)
		gen = retour$gen
	#}
	if(gen@.Data[12] == -1) return(NA)
	return(gen@.Data[9] - gen@.Data[12])
}

gen.nomen = function(gen)#, check = 1)
{
	#if(length(check) != 1) stop("Invalid 'check' parameter: choices are 0 or 1")
		#stop("Param\350tre 'check' invalide: les choix disponibles sont 0 et 1")
	#if(check == 1) {
		retour = gen.detectionErreur(gen = gen, check = 3)
		if(retour$erreur == T)
			stop(retour$messageErreur)
		gen = retour$gen
	#}
	if(gen@.Data[12] == -1) return(NA)
	return(gen@.Data[12])
}

gen.noind = function(gen)#, check = 1)
{
	#if(length(check) != 1) stop("Invalid 'check' parameter: choices are 0 or 1")
		#stop("Param\350tre 'check' invalide: les choix disponibles sont 0 et 1")
	#if(check == 1) {
		retour = gen.detectionErreur(gen = gen, check = c(3))
		if(retour$erreur == T)
			stop(retour$messageErreur)
		gen = retour$gen
	#}
	return(gen@.Data[9])
}

gen.occ = function(gen, pro = 0, ancestors = 0, typeOcc = "IND", ...) # check = 1,
{
	#if(length(check) != 1) stop("Invalid 'check' parameter: choices are 0 or 1")
		#stop("Param\350tre 'check' invalide: les choix disponibles sont 0 et 1")
	if(is(gen, "vector"))
		if(length(list(...)) != 2)
			stop("Invalid '...' parameter : 'father' and 'mother' parameter names are obligatory")
			#stop("Param\350tre '...' invalide : indication du nom des param\350tres 'pere' et 'mere' est obligatoire")
	#if(check == 1) {
		retour = gen.detectionErreur(gen, pro = pro, ancestors = ancestors, check = c(1, 5, 11), ...)
		if(retour$erreur == T)
			stop(retour$messageErreur)
		gen = retour$gen
		pro = retour$pro
		ancestors = retour$ancestors
	#}
	#Les probands sont consideres individuellement
	#Les probands sont divises en groupe
	if(is(pro, "GLgroup")) {
		occurences <- matrix(0, nrow = length(ancestors), ncol = length(pro))
		for(i in 1:length(pro))
			occurences[, i] <- GLPrivOcc(gen, pro = pro[[i]], ancestors = ancestors)
		dimnames(occurences) <- list(ancestors, names(pro))
		return(occurences)
	}
	else {
		occurences <- matrix(0, nrow = length(ancestors), ncol = length(pro))
		for(i in 1:length(pro))
			occurences[, i] <- GLPrivOcc(gen, pro = pro[i], ancestors = ancestors)
		dimnames(occurences) <- list(ancestors, pro)
		if(typeOcc == "IND")
			return(occurences)
		else if(typeOcc == "TOTAL") {
			dfResult.occtot = data.sum(as.data.frame(occurences))
			dimnames(dfResult.occtot)[[1]] <- dimnames(occurences)[[1]]
			dimnames(dfResult.occtot)[[2]] <- c("nb.occ")
			return(dfResult.occtot)
		}
		else
			print("Please choose between \"IND\" and \"TOTAL\" for the variable typeOcc.")
	}
}

gen.parent = function(gen, individuals, output = "FaMo", ...)#, check = 1
{
	#if(length(check) != 1) stop("Invalid 'check' parameter: choices are 0 or 1")
		#stop("Param\350tre 'check' invalide: les choix disponibles sont 0 et 1")
	if(is(gen, "vector"))
		if(length(list(...)) != 2)
			stop("Invalid '...' parameter : 'father' and 'mother' parameter names are obligatory")
			#stop("Param\350tre '...' invalide : indication du nom des param\350tres 'pere' et 'mere' est obligatoire")
	#if(check == 1) {
		retour = gen.detectionErreur(gen = gen, individuals = individuals, output = output, check = c(1, 13, 15), ...)
		if(retour$erreur == T)
			stop(retour$messageErreur)
		gen = retour$gen
		individuals = retour$individuals
		output = retour$output
	#}
	PositionProband = match(individuals, gen$ind)
	Meres <- gen$mother[PositionProband]
	Peres <- gen$father[PositionProband]
	Meres <- Meres[!is.na(Meres)]
	Peres <- Peres[!is.na(Peres)]
	Meres <- unique(Meres)
	Peres <- unique(Peres)
	if(output == "FaMo")
		return(list(Fathers=Peres[Peres > 0], Mothers=Meres[Meres > 0]))
	else if(output == "Fa")
		return(Peres[Peres > 0])
	else if(output == "Mo")
		return(Meres[Meres > 0])
}

#gen.phi = function(gen, pro = 0, nbgenerations = 0, print.it = F, named = T, check = 1)
#{
#	if(length(check) != 1)
#		stop("Invalid 'check' parameter: choices are 0 or 1")
#		#stop("Param\350tre 'check' invalide: les choix disponibles sont 0 et 1")
#	#if(check == 1) {
#		retour = gen.detectionErreur(gen = gen, pro = pro, nbgenerations = nbgenerations, print.it = print.it, named = 
#			named, check = c(3, 5, 19, 18, 10))
#		if(retour$erreur == T)
#			stop(retour$messageErreur)
#		if(length(retour$pro) < 2)
#			stop("Invalid 'pro' parameter: must be a numerical vector of at least 2 proband")
#			#stop("Param\350tre 'prop' invalide: doit \352tre un vecteur num\351rique de 2 proposants minimum")
#		gen = retour$gen
#		pro = retour$pro
#		nbgenerations = retour$nbgenerations
#		print.it = retour$print.it
#		named = retour$named
#	#}
#	#Structure necessaire pour emmagasiner le resultat la fonction de la dll
#	tmp <- double(length(pro) * length(pro))
#	#extern "C" void SPLUSPhiMatrix(long* Genealogie,long* proband, long *NProband,long *Niveau,double* pdRetour, long *printit)
#	#Call de la fonction en C
#	.Call("SPLUSPhiMatrix", gen@.Data, pro, length(pro), as.integer(nbgenerations), tmp, print.it, specialsok = T)
#	dim(tmp) <- c(length(pro), length(pro))
#	#if(named)
#		dimnames(tmp) <- list(pro, pro)
#	if(print.it) {
#		base <- c(deparse(substitute(gen)), deparse(substitute(pro)), nbgenerations)
#		header.txt <- paste("***   Calls : gen.phi (", base[1], ",", base[2], ",", base[3], ")  ***")
#		cat(header.txt, "\n")
#	}
#	return(invisible(tmp))
#}

gen.phiOver = function(phiMatrix, threshold)
{
	if(!is.matrix(phiMatrix))
		return("erreur on doit avoir une matrice")
	n = dim(phiMatrix)[1]
	phiMatrix[phiMatrix >= 0.5] = 0
	phiMatrix[lower.tri(phiMatrix)] = 0
	ind = dimnames(phiMatrix)[[1]]
	indices = matrix(rep(1:n, each = n), n, n)
	ran = indices[phiMatrix >= threshold]
	col = t(indices)[phiMatrix >= threshold]
	if(is.null(ind))
		ind = 1:n
	else ind = as.numeric(ind)
	data.frame(line = ran, column = col, pro1 = ind[ran], pro2 = ind[col], kinship = phiMatrix[phiMatrix >= threshold])
}

gen.phiMean = function(phiMatrix)#, check = 1)#named = T, 
{
	#if(length(check) != 1) stop("Invalid 'check' parameter: choices are 0 or 1")
		#stop("Param\350tre 'check' invalide: les choix disponibles sont 0 et 1")
	#if(check == 1) {
		retour = gen.detectionErreur(matricephi = phiMatrix, named = T, check = c(28, 10))
		if(retour$erreur == T)
			stop(retour$messageErreur)
		phiMatrix = retour$matricephi
		named = retour$named
	#}
	#Test pour accelerer la procedure
	if(class(phiMatrix) == "matrix") mean(phiMatrix[phiMatrix < 0.5]) else GLapplyPhi(phiMatrix, function(x)
		mean(x[x < 0.5]), named = named)
}

#gen.phiMT = function(gen, pro = 0, nbgenerations = 0, print.it = F, named = T, check = 1)
#{
#	if(length(check) != 1)
#		stop("Invalid 'check' parameter: choices are 0 or 1")
#		#stop("Param\350tre 'check' invalide: les choix disponibles sont 0 et 1")
#	#if(check == 1) {
#		retour = gen.detectionErreur(gen = gen, pro = pro, nbgenerations = nbgenerations, print.it = print.it, named = named, check = c(3, 5, 19, 18, 10))
#		if(retour$erreur == T)
#			return(retour$messageErreur)
#		gen = retour$gen
#		pro = retour$pro
#		nbgenerations = retour$nbgenerations
#		print.it = retour$print.it
#		named = retour$named
#	#}
#	#Structure necessaire pour emmagasiner le resultat la fonction de la dll
#	tmp <- double(length(pro) * length(pro))
#	#extern "C" void SPLUSPhiMatrixMT(long* Genealogie,long* proband,long *NProband,long *Niveau,double* pdRetour, long *printit)
#	.Call("SPLUSPhiMatrixMT", gen@.Data, pro, length(pro), as.integer(nbgenerations), tmp, print.it, specialsok = T)
#	dim(tmp) <- c(length(pro), length(pro))
#	#if(named)
#		dimnames(tmp) <- list(pro, pro)
#	if(print.it) {
#		base <- c(deparse(substitute(gen)), deparse(substitute(pro)), nbgenerations)
#		header.txt <- paste("***   Calls : gen.phiMT (", base[1], ",", base[2], ",", base[3], ")  ***")
#		cat(header.txt, "\n")
#	}
#	return(invisible(tmp))
#}

gen.phi = function(gen, pro, depthmin = (gen.depth(gen)-1), depthmax = (gen.depth(gen)-1), MT = F)#, check = 1)#named = T, 
{
	#if(length(check) != 1)	stop("Invalid 'check' parameter: choices are 0 or 1")
	if(missing(pro))		pro = gen.pro(gen)
	if( depthmin<0 | depthmin>(gen.depth(gen)-1) | depthmax<0 | depthmax>(gen.depth(gen)-1) )
		stop("depthmin and depthmax must be between 0 and (gen.depth(gen)-1)")
	
	retour = gen.detectionErreur( gen=gen, pro=pro, depthmin=depthmin, depthmax=depthmax, print.it=F, named=T, check=c(3,5,20,18,10))
	if(retour$erreur == T)	stop(retour$messageErreur)

	gen		= retour$gen
	pro		= retour$pro
	depthmin	= retour$depthmin
	depthmax	= retour$depthmax
#	print.it	= retour$print.it
	named	= retour$named

	#a faire un peu plus tard
	if(MT) {
	  ecart <- as.integer(depthmax) - as.integer(depthmin) + 1
	  np <- length(pro)
	  npp <- length(pro) * length(pro)
	  #Structure necessaire pour emmagasiner le resultat la fonction de la dll
	  rmatrix <- double(ecart * npp)
	  moyenne <- double(ecart)
	  .Call("SPLUSPhisMT", gen@.Data, pro, length(pro), as.integer(depthmin), as.integer(depthmax), moyenne, rmatrix, FALSE, specialsok=T)
	}
	else {
#	  depthmaxtmp = depthmax
#	  depthmintmp = depthmin
	  liste = list()
	  j = 1
	  for(i in depthmin:depthmax) {
		depthmintmp = i
		depthmaxtmp = i
		ecart <- as.integer(depthmaxtmp) - as.integer(depthmintmp) + 1
		np <- length(pro)
		#Structure necessaire pour emmagasiner le resultat la fonction de la dll
		npp <- length(pro) * length(pro)
		rmatrix <- double(ecart * npp)
		moyenne <- double(ecart)
		print.it=F
		.Call("SPLUSPhis", gen@.Data, pro, length(pro), depthmintmp, depthmaxtmp, moyenne, rmatrix, print.it, specialsok = T)
		dim(rmatrix) <- c(np, np, ecart)

		dimnames(rmatrix) <- list(pro, pro, NULL)
		rmatrix <- drop(rmatrix)
#		if(print.it) {
#			base <- c(deparse(substitute(gen)), deparse(substitute(pro)), depthmintmp, depthmaxtmp)
#			header.txt <- paste("***   Calls : gen.phis (", base[1], ",", base[2], ",", base[3], ",", base[4], ")  ***")
#			cat(header.txt, "\n")
#		}
		liste[[j]] = rmatrix
		j = j + 1
	  }
	  sortie.lst = c()
	  for(i in 1:length(liste))	sortie.lst = c(sortie.lst, liste[[i]])
	  ecart <- as.integer(depthmax) - as.integer(depthmin) + 1
	  np <- length(pro)
	
	  #Structure necessaire pour emmagasiner le resultat la fonction de la dll
	  npp <- length(pro) * length(pro)
	  rmatrix <- double(ecart * npp)
	  rmatrix <- sortie.lst
	}
	dim(rmatrix) <- c(np, np, ecart)

	dimnames(rmatrix) <- list(pro, pro, NULL)
	rmatrix <- drop(rmatrix)
	return(invisible(GLmulti(rmatrix, depth = as.integer(depthmin:depthmax))))
}
# print.it = F,
#gen.phis = function(gen, depthmin, depthmax, pro, named = T, check = 1)
#{
#	if(length(check) != 1)	stop("Invalid 'check' parameter: choices are 0 or 1")
#	if(missing(pro))		pro = gen.pro(gen)
#
#	retour = gen.detectionErreur(gen=gen,pro=pro,depthmin=depthmin,depthmax=depthmax,print.it=FALSE,named=named,check=c(3,5,20,18,10))
#	if(retour$erreur == T)	stop(retour$messageErreur)
#	gen		 = retour$gen
#	pro		 = retour$pro
#	depthmin = retour$depthmin
#	depthmax = retour$depthmax
##	print.it = retour$print.it
#	named	 = retour$named
#
#	#a faire un peu plus tard
#	depthmaxtmp = depthmax
#	depthmintmp = depthmin
#	liste = list()
#	j = 1
#	for(i in depthmintmp:depthmaxtmp) {
#		depthmin = i
#		depthmax = i
#		ecart <- as.integer(depthmax) - as.integer(depthmin) + 1
#		np <- length(pro)
#		#Structure necessaire pour emmagasiner le resultat la fonction de la dll
#		npp <- length(pro) * length(pro)
#		rmatrix <- double(ecart * npp)
#		moyenne <- double(ecart)
#		.Call("SPLUSPhis", gen@.Data, pro, length(pro), depthmin, depthmax, moyenne, rmatrix, print.it, specialsok = T)
#		dim(rmatrix) <- c(np, np, ecart)
#		dimnames(rmatrix) <- list(pro, pro, NULL)
#		rmatrix <- drop(rmatrix)
#		if(print.it) {
#			base <- c(deparse(substitute(gen)), deparse(substitute(pro)), depthmin, depthmax)
#			header.txt <- paste("***   Calls : gen.phis (", base[1], ",", base[2], ",", base[3], ",", base[4], ")  ***")
#			cat(header.txt, "\n")
#		}
#		liste[[j]] = rmatrix
#		j = j + 1
#	}
#	sortie.lst = c()
#	for(i in 1:length(liste))
#		sortie.lst = c(sortie.lst, liste[[i]])
#	ecart <- as.integer(depthmaxtmp) - as.integer(depthmintmp) + 1
#	np <- length(pro)
#	#Structure necessaire pour emmagasiner le resultat la fonction de la dll
#	npp <- length(pro) * length(pro)
#	rmatrix <- double(ecart * npp)
#	rmatrix <- sortie.lst
#	dim(rmatrix) <- c(np, np, ecart)
#	#if(named)
#		dimnames(rmatrix) <- list(pro, pro, NULL)
#	rmatrix <- drop(rmatrix)
#	return(invisible(GLmulti(rmatrix, depth = as.integer(depthmintmp:depthmaxtmp))))
#}

#gen.phisMT = function(gen, depthmin, depthmax, pro, print.it = F, named = T, check = 1)
#{
#	if(length(check) != 1)
#		stop("Invalid 'check' parameter: choices are 0 or 1")
#		#stop("Param\350tre 'check' invalide: les choix disponibles sont 0 et 1")
#	if(missing(pro))
#		pro = gen.pro(gen)
#	#if(check == 1) {
#		retour = gen.detectionErreur(gen = gen, pro = pro, depthmin = depthmin, depthmax = depthmax, print.it = print.it,
#			named = named, check = c(3, 5, 20, 18, 10))
#		if(retour$erreur == T)
#			stop(retour$messageErreur)
#		gen = retour$gen
#		pro = retour$pro
#		depthmin = retour$depthmin
#		depthmax = retour$depthmax
#		print.it = retour$print.it
#		named = retour$named
#	#}
#	#a faire un peu plus tard
#	ecart <- as.integer(depthmax) - as.integer(depthmin) + 1
#	np <- length(pro)
#	npp <- length(pro) * length(pro)
#	#Structure necessaire pour emmagasiner le resultat la fonction de la dll
#	rmatrix <- double(ecart * npp)
#	moyenne <- double(ecart)
#	#extern "C" void SPLUSPhis(long* Genealogie,long* proband, long *NProband,long *NiveauMin,long *NiveauMax,double* pdRetour, double *MatrixArray, long *printit)
#	.Call("SPLUSPhisMT", gen@.Data, pro, length(pro), as.integer(depthmin), as.integer(depthmax), moyenne, rmatrix, print.it, specialsok = T)
#	dim(rmatrix) <- c(np, np, ecart)
#	#if(named)
#		dimnames(rmatrix) <- list(pro, pro, NULL)
#	rmatrix <- drop(rmatrix)
#	if(print.it) {
#		base <- c(deparse(substitute(gen)), deparse(substitute(pro)), depthmin, depthmax)
#		header.txt <- paste("***   Calls : gen.phis (", base[1], ",", base[2], ",", base[3], ",", base[4], ")  ***")
#		cat(header.txt, "\n")
#	}
#	return(invisible(GLmulti(rmatrix, depth = as.integer(depthmin:depthmax))))
#}

gen.depth = function(gen)
{
	return(depth(gen))
}

gen.pro = function(gen, ...) #, check = 1
{
	#if(length(check) != 1) stop("Invalid 'check' parameter: choices are 0 or 1")
		#stop("Param\350tre 'check' invalide: les choix disponibles sont 0 et 1")
	if(is(gen, "vector"))
		if(length(list(...)) != 2)
			stop("Invalid '...' parameter : 'father' and 'mother' parameter names are obligatory")
			#stop("Param\350tre '...' invalide : indication du nom des param\350tres 'pere' et 'mere' est obligatoire")
	#print("genPro : 1e verifications faites")
	#if(check == 1) {
		retour = gen.detectionErreur(gen = gen, check = 1, ...)
		#return(1)
		if(retour$erreur)
			return(retour$messageErreur)
		gen = retour$gen
		#print(paste("genPro",retour$erreur))
	#}
	#print(paste("gen.pro post",length(gen$ind)))
	#print(paste("gen.pro post",length(gen$father)))
	#print(paste("gen.pro post",length(gen$mother)))
	#print(paste("gen.pro post",length(gen$sex)))
	return(sort(gen$ind[is.na(match(gen$ind, c(gen$father, gen$mother)))]))
}

gen.rec = function(gen, pro = 0, ancestors = 0, ...) #, check = 1
{
	#if(length(check) != 1) stop("Invalid 'check' parameter: choices are 0 or 1")
		#stop("Param\350tre 'check' invalide: les choix disponibles sont 0 et 1")
	if(is(gen, "vector"))
		if(length(list(...)) != 2)
			stop("Invalid '...' parameter : 'father' and 'mother' parameter names are obligatory")
			#stop("Param\350tre '...' invalide : indication du nom des param\350tres 'pere' et 'mere' est obligatoire")
	#if(check == 1) {
		retour = gen.detectionErreur(gen = gen, pro = pro, ancestors = ancestors, check = c(1, 5, 11), ...)
		if(retour$erreur == T)
			stop(retour$messageErreur)
		gen = gen.genealogy(retour$gen)#, check = 0)
		pro = retour$pro
		ancestors = retour$ancestors
	#}
	if(is(pro, "GLgroup")) {
		nombreAncetre <- length(ancestors)
		nombreGroupe <- length(pro)
		rec <- matrix(0, nrow = nombreAncetre, ncol = nombreGroupe)
		for(i in 1:nombreGroupe) {
			contr <- t(gen.gc(gen, pro[[i]], ancestors))
			rec[, i] <- (contr > 0) %*% rep(1, dim(contr)[2])
		}
		dimnames(rec) <- list(ancestors, names(pro))
		return(rec)
	}
	else {
		contr <- t(gen.gc(gen, pro, ancestors))
		recouv <- (contr > 0) %*% rep(1, dim(contr)[2])
		return(recouv)
	}
}

gen.meangendepthVar = function(gen, pro = 0, type = "MEAN", ...)#, check = 1, named = T, 
{
	#Validation des parametres
	#if(length(check) != 1) stop("Invalid 'check' parameter: choices are 0 or 1")
	if(is(gen, "vector"))
		if(length(list(...)) != 2)
			stop("Invalid '...' parameter : 'father' and 'mother' parameter names are obligatory")

	retour <- gen.detectionErreur(gen = gen, pro = pro, typecomp = type, check = c(1, 5, 17))
	if(retour$erreur == T)	stop(retour$messageErreur)
	gen <- retour$gen
	pro <- retour$pro
	type <- retour$typecomp

	if(type == "IND") {# | type == "MOYSUJETS") {
		tableau <- sapply(pro, function(x, gen, pro, genNo, T)
		GLPriv.variance3V(gen$ind, gen$father, gen$mother, pro = x), gen = gen, pro = pro)
		tableau <- data.frame(tableau)
#		if(type == "MOYSUJETS") {
#			tableau <- data.frame(apply(tableau, 2, mean))
#			dimnames(tableau)[[1]] <- "Mean.Exp.Gen.Depth"
#		}
		if(type == "IND")
			dimnames(tableau)[[1]] <- as.character(paste("Ind", as.character(pro)))
		dimnames(tableau)[[2]] <- "Mean.Gen.Depth"
		return(tableau)
	}
	else if(type == "MEAN")
		return(GLPriv.variance3V(gen$ind, gen$father, gen$mother, pro = pro))
}

#gen.entropyVar2 = function(gen, pro = 0, typeCorpus = "ECH", bfacteurCorr = F, N = NULL, check = 1, ...)
#{
#	#Validation des parametres
#	if(length(check) != 1) stop("Invalid 'check' parameter: choices are 0 or 1")
#		#stop("Param\350tre 'check' invalide: les choix disponibles sont 0 et 1")
#	if(bfacteurCorr == T)
#		if(sum(N) == 0)
#			stop("Correction factor must have a numerical population size value N")
#			#stop("Le facteur de correction doit avoir une valeur num\351rique N taille de la population")
#	if(is(gen, "vector"))
#		if(length(list(...)) != 2)
#			stop("Invalid '...' parameter : 'father' and 'mother' parameter names are obligatory")
#			#stop("Param\350tre '...' invalide : indication du nom des param\350tres 'pere' et 'mere' est obligatoire")
#	#if(check == 1) {
#		retour = gen.detectionErreur(gen = gen, pro = pro, check = c(1, 5))
#		if(retour$erreur == T)
#			stop(retour$messageErreur)
#		gen = retour$gen
#		pro = retour$pro
#	#}
#	tableau = sapply(pro, function(x, gen)
#	GLPriv.variance3V(gen$ind, gen$father, gen$mother, pro = x), gen = gen)
#	if(typeCorpus == "ECH")
#		facteurCorr = length(pro)/(length(pro) - 1)
#	else if(typeCorpus == "POP")
#		facteurCorr = 1
#	if(bfacteurCorr == T)
#		facteurCorr = (facteurCorr * (N - length(pro)))/N
#	tableau = data.frame(tableau)
#	tableau = data.frame(apply(tableau, 2, var) * facteurCorr)
#	dimnames(tableau)[[1]] <- "Prof.varEntropie.var"
#	dimnames(tableau)[[2]] <- "Prof.varEntropie.var"
#	return(tableau)
#}
