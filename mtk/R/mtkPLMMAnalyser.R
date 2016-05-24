# Mexico Toolkit
# Author(s) : R. Faivre, INRA-MIA Toulouse 
# Repository: https://mulcyber.toulouse.inra.fr/projects/baomexico/
# This file was generated automatically by the tool: mtk.AnalyserAddons() built by Juhui WANG, INRA-MIA, JOUY, FRANCE.


#' A sub-class of the class \code{\linkS4class{mtkAnalyser}} used to perform the sensitivity  analysis
#' with the "PLMM" method defined in the "newplmm.R" file.
#' For more details, see the help of the function "plmm()" defined in the "newplmm.R" file.
#' @title The mtkPLMMAnalyser class
#' @exportClass mtkPLMMAnalyser

setClass("mtkPLMMAnalyser",
		contains=c("mtkAnalyser")
	)

#' The constructor.
#' @param mtkParameters a vector of [\code{\linkS4class{mtkParameter}}] representing the parameters necessary to run the Analyser.
#' @param listParameters a named list defining the parameters necessary to run the Analyser. It gives non object-oriented way to specify the parameters.
#' @return an object of class \code{\linkS4class{mtkPLMMAnalyser}}
#' @examples mtkPLMMAnalyser()
#' @export mtkPLMMAnalyser
#' @title The constructor

mtkPLMMAnalyser<- function(mtkParameters=NULL, listParameters=NULL) {
	p<-mtkParameters
		if(!is.null(listParameters))
			p <- make.mtkParameterList(listParameters)
					
	res <- new("mtkPLMMAnalyser",service="PLMM", parameters=p)
		return(res)
		}


#' Performs  sensitivity analysis  with the method  "PLMM" defined in the "newplmm.R" file. 
#' @title The run method
#' @param this an object of class \code{\linkS4class{mtkPLMMAnalyser}}
#' @param context an object of class \code{\linkS4class{mtkExpWorkflow}}
#' @return invisible()
#' @exportMethod run
	
setMethod(f="run", signature=c(this="mtkPLMMAnalyser",
			context="mtkExpWorkflow"),
		definition=function(this, context){
			if(this@state) return(invisible())
			nameThis<-deparse(substitute(this))
			
			X <- context@processesVector$design@result@main 
			Y <- context@processesVector$evaluate@result@main
			parameters<-getParameters(this)
			parameters<-c(context@processesVector$design@result@information,parameters)
			##!!
			##!! Pre-processing the input data, the processing of the method to implement follows:
			##!!
	
	
					analysisOutput<-eval(do.call("plmm",c(list(X=X,Y=Y), parameters)))

			##!!
			##!!  post-processing the output of the method:
			##!!
	
		this@result <- mtkPLMMAnalyserResult(main=analysisOutput$main, information=analysisOutput$information)
		this@state<-TRUE
			
			assign(nameThis, this, envir=parent.frame())
			return(invisible())
		})
	
				####################### 
				## THE ORIGINAL CODE ## 
				####################### 


plmm <- function(X, Y, degree.pol = 1, rawX = FALSE, numY = 1, listeX = NULL, all = FALSE, which = "best", lang = "en", digits = options()$digits, colors = c("red", "orange","blue"), legend.loc = NULL, ...)
{
	if(ncol(Y)>1) { nomY <- names(Y)[numY] ; Y <- data.frame(Y[,numY]) ; names(Y) <- nomY }
	if(!is.null(listeX)) X <- X[ , listeX]
	DATA <- cbind(X,Y)
	resultats <- plmm.mtk(DATA,degpol = degree.pol, raw = rawX)
	information <- list(AnalysisMethod ="PLMM", nomsX = names(X),nomY = resultats$best$nomY, degree.pol = degree.pol, rawX = rawX, numY = numY, listeX = listeX, all = all, which = which, lang = lang, digits = digits, colors = colors, legend.loc = legend.loc)
	resultat <- list(main = resultats, information = information)
	class(resultat) <- c(class(resultat), "plmm")
	resultat
}
#################################################################################
#################################################################################
plmm.mtk <- function (data, degpol = 1, raw = FALSE) 
{
	posY <- ncol(data)
	nomY <- names(data)[posY]
        X <- data[, -posY]
        Y <- data[, posY]
        nomsX <- names(X)
#################################################################################
#################################################################################
	GenerePlanNiveau <- function(degre = 2, facteurs = LETTERS[1:5]) {
		addFacteur <- function(mat, new, degre) {
			rownames(mat) <- NULL
			mat2 <- cbind(mat, new = rep(0:degre, rep(nrow(mat), degre + 1)))
			colnames(mat2)[ncol(mat2)] <- new
			mat2[apply(mat2, 1, sum) < c(degre + 1), ]
		}
 		suite <- list()
		for (i in 1:2) suite[[i]] <- seq(0, degre)
		names(suite) <- facteurs[1:2]
		mat0 <- expand.grid(suite)
		mat1 <- mat0[apply(mat0, 1, sum) < c(degre + 1), ]
		if (length(facteurs) > 2) {
			for (j in 3:length(facteurs)) {
				mat2 <- addFacteur(mat1, facteurs[j], degre)
				mat1 <- mat2
			}
		}
		rownames(mat1) <- apply(mat1, 1, paste, collapse = "")
		mat1
	}
#################################################################################
	GenereXplmm <- function(plan, croisement, degre = 3, RAW = FALSE) {
## Scaling the design to prevent from computational side effetcs
		plan <- scale(plan)
		grosX <- array(apply(plan, 2, poly, degree = degre, raw = RAW), dim = c(nrow(plan), degre, ncol(plan)))
		composeX <- function(v, nx = nrow(plan), val = grosX) {
			if (sum(v) == 0) 
				rep(1, nx)
			else {
				qui <- seq(1, length(v))[v > 0]
				quoi <- v[v > 0]
				reste <- val[, quoi[1], qui[1]]
				if (length(qui) > 1) {
					for (j in 2:length(qui)) reste <- reste * val[, quoi[j], qui[j]]
				}
				reste
			}
		}
		apply(croisement, 1, composeX)
	}
#################################################################################
	calc.plm <- function(Y, X, gpn, ssbase) {
		resultat <- matrix(ncol = 4, nrow = ncol(gpn))
		dimnames(resultat) <- list(names(gpn), c("Main", "Sans", "Total","Partiel"))
#
## utilisation de as.matrix pour eviter l'utilisation d'un vecteur au lieu d'une matrice
		for (i in 1:ncol(gpn)) {
			luiseul <- !(apply(gpn[, -i, drop = FALSE], 1, sum) >  0)
			euxseuls <- !(gpn[, i] > 0)
			formuleM <- lm.fit(as.matrix(X[, luiseul]), Y)
			formuleS <- lm.fit(as.matrix(X[, euxseuls]), Y)
## Par le R2 partiel
			formulePX <- lm.fit(as.matrix(X[, euxseuls]), as.matrix(X[, luiseul]))
			formulePX$residuals <- as.matrix(formulePX$residuals)
			formulePX$residuals[,1] = 1

			formuleP <- lm.fit(formulePX$residuals, formuleS$residuals)

			resultat[i, 1] <- 1 - sum(formuleM$residuals^2)/ssbase
			resultat[i, 2] <- 1 - sum(formuleS$residuals^2)/ssbase

			resultat[i, 4] <- 1 - sum(formuleP$residuals^2)/sum(formuleS$residuals^2)
		}
		resultat[, 3] <- modmax - resultat[, 2]
#
		return(resultat)
	}
#################################################################################
	calc.plmadj <- function(sortie = full) {
#
		r2adjusted <- function(val, n, p) 1 - (1-val)*(n-1)/(n-p-1)
		Cpk <- function(p = 3, k = 3) gamma(k+p+1)/gamma(k+1)/gamma(p+1)
#
		size  <- sortie$size
		dlib <- sortie$dlib
		nbf <- length(sortie$nomsX)
		modmax <- 1 - sortie$residual
		se.base <- sortie$se.base
		X.gpn <- sortie$gpn
#     
## Les informations se trouvent dans la structure X.gpn
# dlib.base contient le nombre de parametres estimes par modele (constante exclue)
		dlib.base <- matrix(NA, ncol = nbf, nrow = 2)
		for (i in 1:ncol(X.gpn)) {
			luiseul <- !(apply(X.gpn[, -i, drop = FALSE], 1, sum) >  0)
			euxseuls <- !(X.gpn[, i] > 0)
			dlib.base[1,i] <- sum(luiseul) -1
			dlib.base[2,i] <- sum(euxseuls) -1 
		}
#

######### Pour eviter une division par 0
	if(dlib>0) r2adj <- r2adjusted(modmax, size, size - dlib -1) else r2adj <- modmax
######### 

		main <- r2adjusted(sortie$main, size, dlib.base[1,])
		sans <- r2adjusted(sortie$sans, size - dlib.base[2,], dlib.base[1,])
		total <- r2adj - sans
## if(positif) total <- replace(total,total<0, 0)
		confusion <- total - main
		partiel <- r2adjusted(sortie$partiel, size - dlib.base[2,] -1, dlib.base[1,]) 
## if(positif) partiel <- replace(partiel,partiel<0, 0)
		specific <- apply(cbind(main, total), 1, min)
##
		output <- list(nomsX = sortie$nomsX, nomY = sortie$nomY, residual = 1 - modmax, main = main, 
			sans = sans, total = total, partiel = partiel, degree = sortie$degree, 
			dlib = dlib, rse = sortie$rse, se.base = sortie$se.base, size = size, gpn = X.gpn, 
			raw = sortie$raw, best = sortie$best, adjusted = TRUE, call = sortie$call)
#
		return(output)
	}
#################################################################################
	plmm.model <- function(nomY,nomX,degmod = 4) {
# Si la dernière colonne du dataframe est la variable à expliquer
#    nomY=names(dataframe)[ncol(dataframe)] 
#    nomX=names(dataframe)[-ncol(dataframe)]
		nom <- paste(nomY, "~ polym(", nomX[1])
		for(i in nomX[-1]) nom <- paste(nom, ",", i)
		nom <- paste(nom, ", degree =",degmod,")")
		nom
	}
#################################################################################
#################################################################################

	xx.gpn <- GenerePlanNiveau(degre = degpol, facteurs = nomsX)
	xx.X <- GenereXplmm(plan = X, croisement = xx.gpn, degre = degpol, RAW = raw)

	size <- length(Y)

	se.base <- sqrt(var(Y))

	expand.X <- as.data.frame(xx.X[,-1])
##    complet <- lm.fit(xx.X, Y)
	complet <- lm(Y ~ ., data = expand.X)
	dlib <- complet$df.residual
	sss <- sum(complet$residuals^2)
######### Pour eviter une division par 0
	if(dlib>0) rss <- sss/dlib else rss <- sss
#########
	rse <- sqrt(rss)
	ssbase <- (size - 1) * var(Y)
	modmax <- 1 - sss/ssbase
	call <- formula(plmm.model(nomY, nomsX, degmod = degpol))

	resultat <- calc.plm(Y, xx.X, xx.gpn, ssbase)

	full <- list(nomsX = nomsX, nomY = nomY, residual = 1 - modmax, main = resultat[, 1], 
		sans = resultat[, 2], total = resultat[, 3], partiel = resultat[, 4], degree = degpol, 
		dlib = dlib, rse = rse, se.base = sqrt(var(Y)), size = size, gpn = xx.gpn, raw = raw, 
		best = FALSE, adjusted = FALSE, call = call)
## Fin de la methode standard
## Debut du calcul avec le R2 ajuste

	full.R2adjusted <- calc.plmadj(full)

## Fin du calcul avec le R2 ajuste

## Recuperation d'une partie de stepwise pour avoir acces au dataframe dans l'appel a stepAIC
	require(MASS)
	rhs <- paste(c("~", deparse(formula(complet)[[3]])), collapse = "")
	rhs <- gsub(" ", "", rhs)
	mod <- update(complet, . ~ 1)
	lower <- ~1
	upper <- eval(parse(text = rhs))

	suite <- stepAIC(mod, scope = list(lower = lower, upper = upper), direction = "both", 
		k = log(length(Y)), trace=0)
# remplace suite = stepwise(complet, direction="forward/backward", trace=0)
	dlib <- suite$df.residual
	sss <- sum(suite$residuals^2)
	rss <- sss/dlib
	rse <- sqrt(rss)
#    ssbase <- (size - 1) * var(Y)
	modmax <- 1 - sss/ssbase
	call <- suite$call

	select.X <- pmatch(names(complet$coef), names(suite$coeff))
	xs.gpn <- xx.gpn[ !is.na(select.X),]
	xs.X <- xx.X[ , !is.na(select.X)]

	resultat <- calc.plm (Y, xs.X, xs.gpn, ssbase)

	best <- list(nomsX = names(xx.gpn), nomY = full$nomY, residual = 1 - modmax, main = resultat[, 
		1], sans = resultat[, 2], total = resultat[, 3], partiel = resultat[, 4], degree = degpol, 
		dlib = dlib, rse = rse, se.base = full$se.base, size = size, gpn = xs.gpn, raw = raw, 
		best = TRUE, adjusted = FALSE, call = call)

## Debut du calcul avec le R2 ajuste sur le modele best
	best.R2adjusted <- calc.plmadj(best)
## Fin du calcul avec le R2 ajuste  sur le modele best


	class(full) <- c(class(full), "plmm")
	class(full.R2adjusted) <- c(class(full.R2adjusted), "plmm")
	class(best) <- c(class(best), "plmm")
	class(best.R2adjusted) <- c(class(best.R2adjusted), "plmm")
	sortie <- list( best = best , full = full, best.adjustedR2 = best.R2adjusted, full.adjustedR2 = 		full.R2adjusted)
	class(sortie) <- c(class(sortie), "plmm")
	return(sortie)
}
####################################################################################################
####################################################################################################
summary.plmm <- function (result, all = result$information$all, which = result$information$which, lang = result$information$lang, digits = result$information$digits, ...) 
{
#    Le contenu est object = result$main
	liste <- which
	if(all) liste <- c("best","full","best.adjustedR2","full.adjustedR2")
#
	for(jj in liste) {
		object <- result$main[[jj]]
		add.text1 <- ""
## english version
		if (lang == "en") {
#main
			if(object$adjusted)  {
				cat(paste("Percentage of factor contributions (adjusted multiple R-squared):", 						format(100 * (1 - object$residual), digits = digits)) ) 
			}
			else {
			cat(paste("Percentage of factor contributions (multiple R-squared):", 
            			format(100 * (1 - object$residual), digits = digits)) )
			}
#model
			if(object$best) add.text1 <- " (stepwise selection) "
			cat(paste("\n\n Polynomial Linear MetaModel of degree", object$degree, add.text1, "\n\n"))
			if(substring(deparse(object$call)[1],1,3) == "lm(") {
				bidon <- deparse(object$call)
				bidon[1] <- paste(object$nomY, substring(bidon[1], 15, nchar(bidon[1]) ))
	        		if(object$raw) bidon[length(bidon)] <- paste(substring(bidon[length(bidon)], 1,
					nchar(bidon[length(bidon)]) -3), "raw X") else
				bidon[length(bidon)] <- paste(substring(bidon[length(bidon)], 1,
					nchar(bidon[length(bidon)]) -3),"X")
				cat(paste("",bidon[1]))
				if(length(bidon)>1) {
					for(i in 2:length(bidon)) cat(paste("\n",bidon[i]))
					}
				}
			else {
				cat(paste("", deparse(object$call)))
			}
#stats
			cat(paste("\n\nInitial  standard error of", format(object$se.base, digits = digits), 
				"with", object$size - 1, "degrees of freedom"))
			cat(paste("\nResidual standard error of", format(object$rse, digits = digits), 
				"with", object$dlib, "degrees of freedom", "\n\n"))
			alone <- apply(cbind(object$main, object$total), 1, min)
			total <- apply(cbind(object$main, object$total), 1, max)
			confusion <- (object$total - object$main)
#summary
			if(object$adjusted) cat("Indices based on adjusted R-squared\n\n")
			print(round(100 * cbind("First Order SI" = alone, "Total SI" = total, 
					Interaction = confusion), digits=digits))
		}
## version française
		if (lang == "fr") {
#titre
			if(object$adjusted)  cat(
				paste("Pourcentage des contributions des facteurs (R2 multiple ajust\u00e9) :", 
					format(100 * (1 - object$residual), digits = digits)) ) else
				cat(paste("Pourcentage des contributions des facteurs (R2 multiple) :", 
					format(100 * (1 - object$residual), digits = digits)) )
#modele
			if(object$best) add.text1 <- " (s\u00e9lection pas \u00e0 pas) "
			cat(paste("\n\n M\u00e9tamod\u00e8le Lin\u00e9aire Polynomial de degr\u00e9", object$degree, add.text1, "\n\n"))
			if(substring(deparse(object$call)[1],1,3) == "lm(") {
				bidon <- deparse(object$call)
				bidon[1] <- paste(object$nomY, substring(bidon[1], 15, nchar(bidon[1]) ))
				if(object$raw) bidon[length(bidon)] <- paste(substring(bidon[length(bidon)], 1,
						nchar(bidon[length(bidon)]) -3), "raw X") else
					bidon[length(bidon)] <- paste(substring(bidon[length(bidon)], 1,
						nchar(bidon[length(bidon)]) -3),"X")
				cat(paste("",bidon[1]))
				if(length(bidon)>1) {
					for(i in 2:length(bidon)) cat(paste("\n",bidon[i]))
					}

			}
			else {
				cat(paste("", deparse(object$call)))
			}
#stats
			cat(paste("\n\nErreur standard initiale   de", format(object$se.base, digits = digits), 
				"avec", object$size - 1, "degr\u00e9s de libert\u00e9"))
			cat(paste("\nErreur standard r\u00e9siduelle de", format(object$rse, digits = digits), 
				"avec", object$dlib, "degr\u00e9s de libert\u00e9", "\n\n"))
			alone <- apply(cbind(object$main, object$total), 1, min)
			total <- apply(cbind(object$main, object$total), 1, max)
			confusion <- (object$total - object$main)
#resume
			if(object$adjusted) cat("Indices bas\u00e9s sur les R2 ajust\u00e9s\n\n")
			print(round(100 * cbind("Indice de 1er ordre" = alone, "Indice total" = total, 
				Interaction = confusion), digits = digits))
		}   
		if(length(liste) > 1 & match(jj,liste) <4) cat(paste("\n",paste(rep(".",60),collapse = ""),"\n\n"))
# Fin de la boucle sur la liste de résumés
		}
	invisible(100 * cbind(Seul = object$main, Specifique = alone, Total = total, Interaction = confusion))
}
####################################################################################################
####################################################################################################
plot.plmm <- function (result, ylim = c(0, 1), colors = result$information$colors, all = result$information$all, which = result$information$which,  lang = result$information$lang, legend.loc = result$information$legend.loc, ...) 
{
#    Le contenu est object = result$main
	liste.base <- c("best","full","best.adjustedR2","full.adjustedR2")
	liste <- which
	oldpar <- par()$mfrow
	par(mfrow = c(1,1))
	if(all) { liste <- liste.base ; par(mfrow = c(2,2)) }
	degmax = result$main[[1]]$degree
	titles <- paste(c("Best metamodel", "Full polynomial metamodel"),
		paste("(degree ",degmax,")",sep="") )
	title.add <- c(titles, paste(titles,"(adjusted R2)", sep="\n"))
	titres <- paste(c("Meilleur m\u00e9tamod\u00e8le", "M\u00e9tamod\u00e8le polynomial complet"), 
		paste("(degr\u00e9 ",degmax,")",sep="") )

	titre.add <- c(titres, paste(titres,"(R2 ajust\u00e9)", sep="\n"))
#
	for(jj in liste) {
		object <- result$main[[jj]]
		S <- rbind(object$main, 1 - object$residual - object$sans - object$main)
		colnames(S) <- object$nomsX
		barplot(S, ylim = ylim, col = colors, ...)
		for (i in 1:length(object$nomsX)) segments(0.1 + 1.2 * (i - 1), 
			object$main[i], 1.2 * i + 0.1, object$main[i], lwd = 2, col = colors[3])
		abline(h = 1 - object$residual,lty = 2)
		if (lang == "en") {
			title (main = title.add[match(jj,liste.base)])
			if(!is.null(legend.loc)) legend(legend.loc, c("main effect", "interactions"), fill = colors)
		}
		if (lang == "fr") {
			title (main = titre.add[match(jj,liste.base)])
			if(!is.null(legend.loc)) legend(legend.loc, c("effet principal", "interactions"), fill = colors)
		}
	}
	par(mfrow = oldpar)
	invisible()
}
####################################################################################################
####################################################################################################

# Mexico Toolkit
# Author(s) : R. Faivre, INRA-MIA Toulouse 
# Repository: https://mulcyber.toulouse.inra.fr/projects/baomexico/
# This file was generated automatically by the tool: mtk.AnalyserAddons() built by Juhui WANG, INRA-MIA, JOUY, France.


#' A sub-class of the class \code{\linkS4class{mtkAnalyserResult}} used to hold the results of the sensitivity  analysis
#' with the "PLMM" method defined in the "newplmm.R" file.
#' For more details, see the help of the function "plmm()" defined in the "newplmm.R" file.
#' @title The mtkPLMMAnalyserResult class
#' @exportClass mtkPLMMAnalyserResult

setClass("mtkPLMMAnalyserResult",
	
				contains=c("mtkAnalyserResult")
		)

#' The constructor.
#'  @param main a data-frame to hold the main results produced by the Analyser.
#'  @param information a named list to provide supplementary information about the analysis process and its results.
				
#' @return an object of class \code{\linkS4class{mtkPLMMAnalyserResult}}
#' @examples mtkmtkPLMMAnalyserResult()
#' @export mtkmtkPLMMAnalyserResult
#' @title The constructor

mtkPLMMAnalyserResult <- function(main, information=NULL) {
	res <- new("mtkPLMMAnalyserResult", main=main, information=information)
				return(res)
				}
#' Shows a summary of the results  produced by the Analyser.
#' For more details, see the help of the function "summary.plmm()" defined in the "newplmm.R" file.
#' @title The summary method
#' @param object an object of class \code{\linkS4class{mtkPLMMAnalyserResult}}
#' @return invisible()
#' @exportMethod summary

setMethod(f="summary", "mtkPLMMAnalyserResult",
	definition=function(object,...){
summary.plmm(list(main=object@main,information=object@information),...)
	})

#' Plots the results  produced by the Analyser.
#' For more details, see the help of the function "plot.plmm()" defined in the "newplmm.R" file.
#' @title The plot method
#' @param x an object of class \code{\linkS4class{mtkPLMMAnalyserResult}}

#' @return invisible()
#' @exportMethod plot

setMethod(f="plot", "mtkPLMMAnalyserResult",
	definition=function(x,y,...){
if(!missing(y)) plot.plmm(list(main=x@main,information=x@information),y, ...)
else plot.plmm(list(main=x@main,information=x@information), ...)
	})
