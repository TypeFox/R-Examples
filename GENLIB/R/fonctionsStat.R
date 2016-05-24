# 1 - gen.fCI					-> garde article
# 2 - gen.phiCI				-> garde article


gen.fCI = function(vectF, prob = c(0.025, 0.05, 0.95, 0.975), b = 5000, print.it = F)#, check = 1)#named = T, 
{
	#if(length(check) != 1) stop("Invalid 'check' parameter: choices are 0 or 1")
		#stop("Param\350tre 'check' invalide: les choix disponibles sont 0 et 1")
	#if(check == 1) {
		retour = gen.detectionErreur(vectF = vectF, prob = prob, b = b, print.it = print.it, named = T, check = c(33, 21, 22, 18, 10))
		if(retour$erreur == T)
			stop(retour$messageErreur)
		vectF = retour$vectF
		prob = retour$prob
		b = retour$b
		print.it = retour$print.it
		named = retour$named
	#}
	FUN <- function(x, pprogress, b, prob)
	{
		drop(t(bcanon(x = x, nboot = b, theta = mean, alpha = prob)$confpoints[,2]))
		#drop(bcanon(x=x, nboot=b, theta=mean, alpha=prob))$confpoints
		#drop(boot.ci(boot(x, mean, R=b), prob)$bca)
		#limits.bca(bootstrap(x, mean, B = b, trace = pprogress), prob))
	}
	if(!is.null(group(vectF))) {
		#Validation si tous les groups sont de meme taille
		tai = sapply(group(vectF), length)
		if(any(tai <= 1))
			stop("Invalid parameter: groups in matricephi must all contain at least 2 probands")
			#stop("Parametre invalide: Les groups contenus dans matricephi doivent tous inclure au moins 2 probands")
	}
	#Test pour accelerer la procedure
	nam <- sapply(prob * 100, function(x)	paste(x, "%", collapse = "1", sep = ""))
	
	return(GLapplyF(vectF, FUN, FunReturnLength = length(prob), named = named, namesVector = nam, pprogress = print.it, b = b, prob = prob))
}

gen.phiCI = function(phiMatrix, prob = c(0.025, 0.05, 0.95, 0.975), b = 5000, print.it = F)#, check = 1)#named = T, 
{
	#if(length(check) != 1) stop("Invalid 'check' parameter: choices are 0 or 1")
		#stop("Param\350tre 'check' invalide: les choix disponibles sont 0 et 1")
	#if(check == 1) {
		retour = gen.detectionErreur(matricephi = phiMatrix, prob = prob, b = b, print.it = print.it, named = T, 
								check = c(29, 21, 22, 18, 10))
		if(retour$erreur == T)
			stop(retour$messageErreur)
		matricephi = retour$matricephi
		prob = retour$prob
		b = retour$b
		print.it = retour$print.it
		named = retour$named
	#}
	#Si les groups sont petits et s'ils ne sont pas de meme taille, une erreur se produit 
	##Il faut verifier avec Louis
	if(!is.null(group(matricephi))) {
		#Validation si tous les groups sont de meme taille
		tai = sapply(group(matricephi), length)
		if(any(tai <= 2))
			stop("Invalid 'matricephi' parameter: groups in 'matricephi' must all contain at least 3 probands")
			#stop("Param\350tre 'matricephi' invalide: Les groups contenus dans 'matricephi' doivent tous inclure au moins 3 probands")
		if(any(tai != tai[1]))
			stop("Invalid 'matricephi' parameter: all groups in 'matricephi' must be the same size")
			#stop("Param\350tre 'matricephi' invalide: Les groups contenus dans 'matricephi' doivent tous etre de la meme taille")
	}
	FUN <- function(x, pprogress, b, prob)
	{
		tmp <- bootstrap(x=1:dim(x)[1], nboot = 5000, theta = function(a, r) { mean(r[a, a][r[a, a] < 0.5]) }, r = x )
		quantile(tmp[[1]], c(0.025, 0.05, 0.95, 0.975))
		#tmp <- boot( 1:dim(x)[1], function(a, r){ mean(r[a, a][r[a, a] < 0.5]) }, R = b, args.stat=list(r = x)) #, trace = pprogress) ,assign.frame1=T
		#limits.emp(tmp, prob = prob)[1,  ]
	}
	#fin Glapplu  param=pprogress=print.it,b=b,prob=prob) 
	#Test pour accelerer la procedure
	if(class(matricephi) == "matrix") {
		ret <- FUN(matricephi, pprogress = print.it, b = b, prob = prob)
		if(!named)
			class(ret) <- "numeric"
		return(ret)
	} else {
		nam <- sapply(prob * 100, function(x)
		paste(x, "%", collapse = "1", sep = ""))
		ret <- GLapplyPhi(matricephi, FUN, FunReturnLength = length(prob), named = named, namesVector = nam, pprogress = 
			print.it, b = b, prob = prob)
		return(ret)
	}
}


