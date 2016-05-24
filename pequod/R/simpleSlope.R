# object=mod5
## moderated regression with mean centering
#library(car)
#data(Ginzberg)

simple.slope <- function(object, pred, mod1, mod2 = "none", 
	coded = "none") {
	regr <- object$Stepfin
	macov <- vcov(regr)
	nam <- dimnames(macov)[[1]]
	
	jj <- strsplit(nam, ".XX.")
	
	
	ord <- "none" %in% mod2
	
	n <- 0
	k <- 0
	ljj <- length(jj)
	for (i in (1:ljj)) {
		n <- n + 1
		k[n] <- length(jj[[i]])
	}
	
	
	fI <- length(jj[k == 1])
	#numero elementi ordine 1 + intercetta
	fII <- length(jj[k == 2])
	#numero elementi ordine 2
	fIII <- length(jj[k == 3])
	#numero elementi ordine 3(NON NECESSARIA DA ELIMINARE)
	
	
	if (ord == TRUE) {
		pos_pred <- grep(paste0("\\b",pred,"\\b"), jj[1:fI])
		#posizione predittore + intercetta
		pos_mod <- grep(paste0("\\b",mod1,"\\b"), jj[1:fI])
		#posizione moderatore + intercetta
		
		
		pos_pred_mod1 <- fI + grep(paste0("\\b",mod1,"\\b"), jj[(fI + 1):(fI + fII)])
		#posizione moderatore + intercetta
		pos_pred_mod2 <- fI + grep(paste0("\\b",pred,"\\b"), jj[(fI + 1):(fI + fII)])
		#posizione moderatore + intercetta 
		posint2 <- sum((pos_pred_mod1 %in% pos_pred_mod2) * pos_pred_mod1)
		#posiz finale
		
		nomX <- names(regr$coef)[pos_pred]
		nomZ <- names(regr$coef)[pos_mod]
		nomY <- names(regr$model)[1]
		
		
		if (nomX %in% coded) {
			#sgamare il coding e sostituire i punti
			tax <- table(regr$model[nomX])
			tax <- as.numeric(dimnames(tax)[[1]])
			
			X_1L <- min(tax)
			X_1H <- max(tax)
			
		}
		else {
			#punti del predittore
			X_sd <- sd(regr$model[[pos_pred]], na.rm = T)
			X_1L <- mean(regr$model[[pos_pred]], na.rm = T) - 
				X_sd
			X_1H <- mean(regr$model[[pos_pred]], na.rm = T) + 
				X_sd
		}
		
		if (nomZ %in% coded) {
			#sgamare il coding e sostituire i punti
			taz <- table(regr$model[nomZ])
			taz <- as.numeric(dimnames(taz)[[1]])
			
			M_1L <- min(taz)
			M_1H <- max(taz)
			
		}
		else {
			#punti del moderatore
			M_sd <- sd(regr$model[[pos_mod]], na.rm = T)
			M_1L <- mean(regr$model[[pos_mod]], na.rm = T) - 
				M_sd
			M_1H <- mean(regr$model[[pos_mod]], na.rm = T) + 
				M_sd
		}
		
		x <- c(X_1L, X_1H)
		z <- c(M_1L, M_1H)
		
		
		
		#loop per le matrici
		n <- 0
		i <- 0
		k <- 0
		matr <- 0
		slope <- 0
		slope_se <- 0
		t_num <- 0
		for (i in (1:length(x))) {
			for (k in 1:(length(z))) {
				n <- n + 1
				# matrice del y in funzione di x e z
				matr[n] <- (regr$coef[[pos_pred]] + (regr$coef[[posint2]] * 
				z[k])) * x[i] + (regr$coef[[1]] + (regr$coef[[pos_mod]] * 
				z[k]))
				# slope di x in funzione di z
				slope[k] <- (regr$coef[[pos_pred]] + (regr$coef[[posint2]] * 
				z[k]))
				# errore standard delle slope
				slope_se[k] <- sqrt(macov[pos_pred, pos_pred] + 
				2 * z[k] * macov[pos_pred, posint2] + macov[posint2, 
				posint2] * (z[k]^2))
			}
		}
		
		# numeratore test t
		t_num <- slope
		
		#output matrice
		hmatr <- matrix(matr, ncol = 2)
		
		#output matrice
		if (nomZ %in% coded) {
			leg1 <- paste("Low", nomZ, "(", M_1L, ")")
			leg3 <- paste("High", nomZ, "(", M_1H, ")")
		}
		else {
			leg1 <- paste("Low", nomZ, "(-1 SD)")
			leg3 <- paste("High", nomZ, "(+1 SD)")
		}
		
		if (nomX %in% coded) {
			pr1 <- paste("Low", nomX, "(", X_1L, ")")
			pr3 <- paste("High", nomX, "(", X_1H, ")")
		}
		else {
			pr1 <- paste("Low", nomX, "(-1 SD)")
			pr3 <- paste("High", nomX, "(+1 SD)")
		}
		
		pr <- c(pr1, pr3)
		mo <- c(leg1, leg3)
		pmatr <- matrix(matr, ncol = 2, dimnames = list(mo, pr))
		t_value <- t_num/slope_se
		dimm <- dim(regr$model)
		df_t <- dimm[1] - dimm[2]
		p_value <- 2 * (1 - pt(abs(t_value), df = df_t))
		
		resu <- cbind(slope, slope_se, t_value, p_value)
		colnames(resu) <- c("simple slope", "standard error", 
			"t-value", "p.value")
		rownames(resu) <- mo
		
		
		
		#Baur & Curran
		a <- ((1.96)^2) * (macov[posint2, posint2]) - ((regr$coef[[posint2]])^2)
		b <- ((1.96)^2) * (macov[pos_pred, posint2]) - ((regr$coef[[posint2]]) * 
			(regr$coef[[pos_pred]]))
		co <- ((1.96)^2) * (macov[pos_pred, pos_pred]) - ((regr$coef[[pos_pred]])^2)
		
		
		
		del <- (b^2) - a * co
		x1 <- (-b - sqrt(del))/a
		x2 <- (-b + sqrt(del))/a
		
		int <- cbind(min(c(x1, x2)), max(c(x1, x2)))
		colnames(int) <- c("lower CI", "upper CI")
		rownames(int) <- nomZ
		
		orde <- 2
		
		
		pippo <- list(nomY = nomY, nomX = nomX, X_1L = X_1L, 
			X_1H = X_1H, orde = orde, Points = pmatr, simple_slope = resu, 
			conf95 = int, Df = df_t)
	}
	
	
	
	
	if (ord == FALSE) {
		
		# ordine 1
		pos_pred <- grep(paste0("\\b",pred,"\\b"), jj[1:fI])
		#posizione predittore + intercetta
		pos_mod <- grep(paste0("\\b",mod1,"\\b"), jj[1:fI])
		#posizione moderatore + intercetta
		pos_mod1 <- grep(paste0("\\b",mod2,"\\b"), jj[1:fI])
		#posizione moderatore 2 + intercetta
		
		
		# ordine 2
		pos_II_mod1 <- fI + grep(paste0("\\b",mod1,"\\b"), jj[(fI + 1):(fI + fII)])
		#posizione moderatore 1 + intercetta
		pos_II_mod2 <- fI + grep(paste0("\\b",mod2,"\\b"), jj[(fI + 1):(fI + fII)])
		#posizione moderatore + intercetta
		pos_II_pred <- fI + grep(paste0("\\b",pred,"\\b"), jj[(fI + 1):(fI + fII)])
		#posizione moderatore + intercetta
		
		pos_pred_mod1 <- intersect(pos_II_pred, pos_II_mod1)
		pos_pred_mod2 <- intersect(pos_II_pred, pos_II_mod2)
		pos_mod1_mod2 <- intersect(pos_II_mod1, pos_II_mod2)
		
		# ordine 3
		pos_III_mod1 <- fI + fII + grep(paste0("\\b",mod1,"\\b"), jj[(fI + fII + 
			1):(fI + fII + fIII)])
		#posizione moderatore 1 + intercetta
		pos_III_mod2 <- fI + fII + grep(paste0("\\b",mod2,"\\b"), jj[(fI + fII + 
			1):(fI + fII + fIII)])
		#posizione moderatore + intercetta
		pos_III_pred <- fI + fII + grep(paste0("\\b",pred,"\\b"), jj[(fI + fII + 
			1):(fI + fII + fIII)])
		#posizione moderatore + intercetta
		
		pos_pred_mod1_mod2 <- intersect(intersect(pos_III_pred, 
			pos_III_mod1), pos_III_mod2)
		
		
		
		
		
		nomX <- names(regr$coef)[pos_pred]
		nomZ <- names(regr$coef)[pos_mod]
		nomW <- names(regr$coef)[pos_mod1]
		nomY <- names(regr$model)[1]
		
		
		
		
		if (nomX %in% coded) {
			#sgamare il coding e sostituire i punti
			tax <- table(regr$model[nomX])
			tax <- as.numeric(dimnames(tax)[[1]])
			
			X_1L <- min(tax)
			X_1H <- max(tax)
		}
		else {
			#punti del predittore
			X_sd <- sd(regr$model[[pos_pred]], na.rm = T)
			X_1L <- mean(regr$model[[pos_pred]], na.rm = T) - 
				X_sd
			X_1H <- mean(regr$model[[pos_pred]], na.rm = T) + 
				X_sd
		}
		
		if (nomZ %in% coded) {
			#sgamare il coding e sostituire i punti
			taz <- table(regr$model[nomZ])
			taz <- as.numeric(dimnames(taz)[[1]])
			M1_1L <- min(taz)
			M1_1H <- max(taz)
		}
		else {
			#punti del predittore
			M1_sd <- sd(regr$model[[pos_mod]], na.rm = T)
			M1_1L <- mean(regr$model[[pos_mod]], na.rm = T) - 
				M1_sd
			M1_1H <- mean(regr$model[[pos_mod]], na.rm = T) + 
				M1_sd
		}
		
		if (nomW %in% coded) {
			#sgamare il coding e sostituire i punti
			taw <- table(regr$model[nomW])
			taw <- as.numeric(dimnames(taw)[[1]])
			M2_1L <- min(taw)
			M2_1H <- max(taw)
		}
		else {
			#punti del predittore
			M2_sd <- sd(regr$model[[pos_mod1]], na.rm = T)
			M2_1L <- mean(regr$model[[pos_mod1]], na.rm = T) - 
				M2_sd
			M2_1H <- mean(regr$model[[pos_mod1]], na.rm = T) + 
				M2_sd
		}
		
		
		
		
		#matrice dei punti
		x <- c(X_1L, X_1H)
		z <- c(M1_1L, M1_1H)
		w <- c(M2_1L, M2_1H)
		
		
		
		#loop per le matrici
		n <- 0
		i <- 0
		k <- 0
		f <- 0
		matr <- 0
		slope <- 0
		slope_se <- 0
		t_num <- 0
		
		for (i in (1:length(x))) {
			for (k in 1:(length(z))) {
				for (f in 1:(length(w))) {
				n <- n + 1
				# matrice del y in funzione di x e z
				matr[n] <- (regr$coef[[pos_pred]] + (regr$coef[[pos_pred_mod1]] * 
					z[k]) + (regr$coef[[pos_pred_mod2]] * w[f]) + 
					(regr$coef[[pos_pred_mod1_mod2]] * w[f] * 
					z[k])) * x[i] + regr$coef[[1]] + (regr$coef[[pos_mod]] * 
					z[k]) + (regr$coef[[pos_mod1]] * w[f]) + 
					(regr$coef[[pos_mod1_mod2]] * z[k] * w[f])
				slope[n] <- (regr$coef[[pos_pred]] + (regr$coef[[pos_pred_mod1]] * 
					z[k]) + (regr$coef[[pos_pred_mod2]] * w[f]) + 
					(regr$coef[[pos_pred_mod1_mod2]] * w[f] * 
					z[k]))
				slope_se[n] <- sqrt(macov[pos_pred, pos_pred] + 
					macov[pos_pred_mod1, pos_pred_mod1] * (z[k]^2) + 
					macov[pos_pred_mod2, pos_pred_mod2] * (w[f]^2) + 
					macov[pos_pred_mod1_mod2, pos_pred_mod1_mod2] * 
					(z[k]^2) * (w[f]^2) + 2 * macov[pos_pred, 
					pos_pred_mod1] * z[k] + 2 * macov[pos_pred, 
					pos_pred_mod2] * w[f] + 2 * macov[pos_pred, 
					pos_pred_mod1_mod2] * z[k] * w[f] + 2 * macov[pos_pred_mod1, 
					pos_pred_mod2] * w[f] * z[k] + 2 * macov[pos_pred_mod1, 
					pos_pred_mod1_mod2] * (z[k]^2) * (w[f]) + 
					2 * macov[pos_pred_mod2, pos_pred_mod1_mod2] * 
					z[k] * (w[f]^2))
				}
			}
		}
		
		
		#punti simple slope differences
		b4 <- regr$coef[[pos_pred_mod1]]
		b5 <- regr$coef[[pos_pred_mod2]]
		b7 <- regr$coef[[pos_pred_mod1_mod2]]
		s44 <- macov[pos_pred_mod1, pos_pred_mod1]
		s55 <- macov[pos_pred_mod2, pos_pred_mod2]
		s77 <- macov[pos_pred_mod1_mod2, pos_pred_mod1_mod2]
		s45 <- macov[pos_pred_mod1, pos_pred_mod2]
		s47 <- macov[pos_pred_mod1, pos_pred_mod1_mod2]
		s57 <- macov[pos_pred_mod2, pos_pred_mod1_mod2]
		ZH <- M1_1H
		ZL <- M1_1L
		WH <- M2_1H
		WL <- M2_1L
		
		d12 <- (b5 + b7 * ZH)/sqrt(s55 + (ZH^2) * s77 + 2 * ZH * 
			s57)
		
		d13 <- (b4 + b7 * WH)/sqrt(s44 + (WH^2) * s77 + 2 * WH * 
			s47)
		d24 <- (b4 + b7 * WL)/sqrt(s44 + (WL^2) * s77 - 2 * WL * 
			s47)
		d34 <- (b5 + b7 * ZL)/sqrt(s55 + (ZL^2) * s77 - 2 * ZL * 
			s57)
		d14 <- (b4 * (ZH - ZL) + b5 * (WH - WL) + b7 * (ZH * 
			WH - ZL * WL))/sqrt(((ZH - ZL)^2) * s44 + ((WH - 
			WL)^2) * s55 + ((ZH * WH - ZL * WL)^2) * s77 + 2 * 
			((ZH - ZL) * (WH - WL) * s45 + (ZH - ZL) * (ZH * 
				WH - ZL * WL) * s47 + (WH - WL) * (ZH * WH - 
				ZL * WL) * s57))
		d23 <- (b4 * (ZH - ZL) + b5 * (WL - WH) + b7 * (ZH * 
			WL - ZL * WH))/sqrt(((ZH - ZL)^2) * s44 + ((WL - 
			WH)^2) * s55 + ((ZH * WL - ZL * WH)^2) * s77 + 2 * 
			((ZH - ZL) * (WL - WH) * s45 + (ZH - ZL) * (ZH * 
				WL - ZL * WH) * s47 + (WL - WH) * (ZH * WL - 
				ZL * WH) * s57))
		
		#output matrice
		if (nomZ %in% coded) {
			leg1 <- paste("Low", nomZ, "(", M1_1L, ")")
			leg3 <- paste("High", nomZ, "(", M1_1H, ")")
		}
		else {
			leg1 <- paste("Low", nomZ, "(-1 SD)")
			leg3 <- paste("High", nomZ, "(+1 SD)")
		}
		
		if (nomW %in% coded) {
			legm1 <- paste("Low", nomW, "(", M2_1L, ")")
			legm3 <- paste("High", nomW, "(", M2_1H, ")")
		}
		else {
			legm1 <- paste("Low", nomW, "(-1 SD)")
			legm3 <- paste("High", nomW, "(+1 SD)")
		}
		
		if (nomX %in% coded) {
			pr1 <- paste("Low", nomX, "(", X_1L, ")")
			pr3 <- paste("High", nomX, "(", X_1H, ")")
		}
		else {
			pr1 <- paste("Low", nomX, "(-1 SD)")
			pr3 <- paste("High", nomX, "(+1 SD)")
		}
		
		r1 <- paste(leg1, ",", legm1, "[1]")
		r2 <- paste(leg1, ",", legm3, "[2]")
		r3 <- paste(leg3, ",", legm1, "[3]")
		r4 <- paste(leg3, ",", legm3, "[4]")
		
		
		
		pr <- c(pr1, pr3)
		mo <- c(r1, r2, r3, r4)
		pmatr <- matrix(matr, ncol = 2, dimnames = list(mo, pr))
		pslope <- matrix(slope[1:4], ncol = 1)
		pslope_se <- matrix(slope_se[1:4], ncol = 1)
		
		
		t_value <- pslope/pslope_se
		dimm <- dim(regr$model)
		df_t <- dimm[1] - dimm[2]
		p_value <- 2 * (1 - pt(abs(t_value), df = df_t))
		
		resu <- cbind(pslope, pslope_se, t_value, p_value)
		rownames(resu) <- mo
		colnames(resu) <- c("simple slope", "standard error", 
			"t-value", "p.value")
		
		
		
		#diff matr
		di12 <- paste(r4, " vs. ", r3)
		di13 <- paste(r4, " vs. ", r2)
		di24 <- paste(r3, " vs. ", r1)
		di34 <- paste(r2, " vs. ", r1)
		di14 <- paste(r4, " vs. ", r1)
		di23 <- paste(r3, " vs. ", r2)
		
		nomo <- c(di12, di13, di14, di24, di23, di34)
		
		difslo <- c(d12, d13, d14, d24, d23, d34)
		
		p_value1 <- 2 * (1 - pt(abs(difslo), df = df_t))
		
		b_value <- p_value1 * 6
		maxx <- c(1, 1, 1, 1, 1, 1)
		bvalue1 <- ifelse(b_value >= maxx, maxx, b_value)
		madif <- cbind(difslo, p_value1, bvalue1)
		rownames(madif) <- nomo
		colnames(madif) <- c("t-value", "p.value", "    Bonferroni.p")
		
		
		
		orde <- 3
		
		pippo = list(nomY = nomY, nomX = nomX, X_1L = X_1L, X_1H = X_1H, 
			orde = orde, Points = pmatr, simple_slope = resu, 
			delta_slope = madif, Df = df_t)
		
	}
	
	return(pippo)
	
}


#########################################################
#                OBJECT simple.slope                    #
#########################################################


simpleSlope <- function(object, pred, mod1, mod2, 
	coded, ...) UseMethod("simpleSlope")

## simpleSlope default ################### 
simpleSlope.default <- function(object, pred, mod1, 
	mod2 = "none", coded = "none", ...) {
	object <- as(object, "lmres")
	pred <- as.character(pred)
	mod1 <- as.character(mod1)
	mod2 <- as.character(mod2)
	coded <- as.character(coded)
	est1 <- simple.slope(object, pred, mod1, mod2, coded)
	est1$call <- match.call()
	class(est1) <- "simpleSlope"
	est1
}

## print simpleSlope  ################### \t
print.simpleSlope <- function(x, ...) {
	cat("Simple Slope:\n")
	print(x$simple_slope)
	
}


#                OBJECT SUMMARY LMRES                   #
#########################################################

summary.simpleSlope <- function(object, ...) {
	
	
	## summary default ##
	
	if (object$orde == 2) {
		orde <- object$orde
		
		#Points
		Points <- object$Points
		Points <- round(Points, 4)
		
		
		#Simple slope
		simple_slope <- object$simple_slope
		simple_slope <- round(simple_slope, 4)
		
		
		#Bauer & Curran
		conf95 <- object$conf95
		conf95 <- round(conf95, 4)
		
		#output list
		sout <- list(nomY = object$nomY, orde = orde, Points = Points, 
			simple_slope = simple_slope, conf95 = conf95, Df = object$Df)
	}
	
	if (object$orde == 3) {
		orde <- object$orde
		
		#Points
		Points <- object$Points
		Points <- round(Points, 4)
		
		
		#Simple slope
		simple_slope <- object$simple_slope
		simple_slope <- round(simple_slope, 4)
		
		delta_slope <- object$delta_slope
		delta_slope <- round(delta_slope, 4)
		
		#output list
		sout <- list(nomY = object$nomY, orde = orde, Points = Points, 
			simple_slope = simple_slope, delta_slope = delta_slope, 
			Df = object$Df)
		
	}
	
	class(sout) <- "summary.simpleSlope"
	sout
}



# print summary.simpleSlope ################### 

print.summary.simpleSlope <- function(x, ...) {
	cat("\n")
	cat("** Estimated points of", x$nomY, " **\n")
	cat("\n")
	printCoefmat(x$Points, justify = "centre")
	cat("\n")
	cat("\n")
	cat("\n")
	cat("** Simple Slopes analysis ( df=", x$Df, ") **\n")
	cat("\n")
	printCoefmat(x$simple_slope, P.values = TRUE, has.Pvalue = TRUE, 
		justify = "centre", digits = 3)
	cat("\n")
	cat("\n")
	cat("\n")
	
	if (x$orde == 2) {
		
		cat("** Bauer & Curran 95% CI **\n")
		cat("\n")
		printCoefmat(x$conf95, justify = "centre")
	}
	
	if (x$orde == 3) {
		
		cat("** Slope Difference Test (( df=", x$Df, "); Dawson & Richter, 2006) **\n")
		
		printCoefmat(x$delta_slope, P.values = TRUE, has.Pvalue = TRUE, 
			justify = "centre")
	}
	
}








