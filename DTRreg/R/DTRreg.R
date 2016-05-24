# Doubly-robust dynamic treatment regimen estimation via dynamic weighted ordinary least squares (dWOLS), assuming linear blip functions and logistic regression for treatment model
#' @export
DTRreg <- function(outcome, blip.mod, treat.mod, tf.mod, data=NULL, method = "gest", weight = "default", var.estim="none", B = 200, M = 0, truncate = 0, verbose = "FALSE", interrupt = "FALSE", treat.range = NULL, missing = "default", interactive = FALSE, treat.mod.man = NULL, type = "DTR") {
	# if interactive mode chosen, build standardized input from user input
	if (interactive==TRUE & var.estim != "bs.quiet") {
		cat("DTR estimation interactive mode!\n")
		outcome <- get(readline("Enter outcome variable: "))
		K <- readline("Enter number of stages: ")
		for (j in 1:K) {
			var <- c()
			new.var <- ""
			cat("Enter each stage",j,"blip covariate (type STOP when done)\n")
			while (new.var != "STOP") {
				new.var <- readline(": ")
				if (new.var != "STOP") {var <- append(var,new.var)}
			}
			# convert input to formula form
			if (length(var) == 0) {var <- "1"}
			if (length(var) == 1) {blip.mod[[j]] <- paste("~",var)}
			else {blip.mod[[j]] <- paste("~",paste(var,collapse="+"))}
			cat("Stage",j,"blip model:",blip.mod[[j]],"\n")
			blip.mod[[j]] <- as.formula(blip.mod[[j]])
		}
		for (j in 1:K) {
			var <- c()
			new.var <- ""
			cat("Enter each stage",j,"treatment-free covariate (type STOP when done)\n")
			while (new.var != "STOP") {
				new.var <- readline(": ")
				if (new.var != "STOP") {var <- append(var,new.var)}
			}
			if (length(var) == 0) {var <- "1"}
			if (length(var) == 1) {tf.mod[[j]] <- paste("~",var)}
			else {tf.mod[[j]] <- paste("~",paste(var,collapse="+"))}
			cat("Stage",j,"treatment-free model:",tf.mod[[j]],"\n")
			tf.mod[[j]] <- as.formula(tf.mod[[j]])

		}
		for (j in 1:K) {
			var <- c()
			new.var <- ""
			cat("Enter stage",j,"treatment variable\n")
			out <- readline(": ")
			cat("Enter each stage",j,"treatment model covariate (type STOP when done)\n")
			while (new.var != "STOP") {
				new.var <- readline(": ")
				if (new.var != "STOP") {var <- append(var,new.var)}
			}
			if (length(var) == 0) {var <- "1"}
			if (length(var) == 1) {treat.mod[[j]] <- paste(out,"~",var)}
			else {treat.mod[[j]] <- paste(out,"~",paste(var,collapse="+"))}
			cat("Stage",j,"treatment model:",treat.mod[[j]],"\n")
			treat.mod[[j]] <- as.formula(treat.mod[[j]])
		}
		method <- readline("Dynamic WOLS (dwols), G-estimation (gest), or Q-learning (qlearn)?\n> ")
		var.estim <- readline("Variance estimation: none, bootstrap or sandwich?\n> ")
	}
	obj <- list()
	# checking input validity
	try(match.arg(method,c("dwols","gest","qlearn")))
	try(match.arg(var.estim,c("none","bootstrap","sandwich","bs.quiet")))
	# if q-learning selected, then set method to dwols and define qlearn to 1 for later
	qlearn <- 0
	if (method == "qlearn") {
		method <- "dwols"
		qlearn <- 1
	}
	if (var.estim == "sandwich" & method == "dwols") {stop("Sandwich variance estimation only available with G-estimation.\n")}
	if (length(treat.mod) != length(blip.mod) | length(treat.mod) != length(tf.mod)) {stop("-treat.mod-, -blip.mod- and -tf.mod- must all have same length.")}
	if (truncate <0 | truncate > 0.5) {stop("-truncate- must lie between 0 and 0.5.")}
	# infer stages from treat.mod
	K <- length(treat.mod)
	obj$K <- K
	# if no data set, create one from the models
	if (is.null(data)) {
		# check for any intercept-only 
		n <- length(outcome)
		data <- cbind(outcome,do.call("cbind",lapply(blip.mod,get_all_vars_DTR,n)),do.call("cbind",lapply(treat.mod,get_all_vars_DTR,n)),do.call("cbind",lapply(tf.mod,get_all_vars_DTR,n)))
		data <- data[!duplicated(lapply(data,summary))]
		obj$outcome <- outcome
	} else {
		obj$outcome <- eval(substitute(outcome),envir=data)
	}
	# make sure data is a data frame (can't be a matrix)
	data <- data.frame(data)
	obj$data <- data

	obj$blip.mod <- blip.mod
	obj$treat.mod <- treat.mod
	obj$tf.mod <- tf.mod
	Y <- obj$outcome
	if (any(is.na(Y))) {
		stop("Missing values in the outcome are not allowed.")
	}
	obj$obs.Y <- Y
	obj$blip.list <- obj$psi <- obj$covmat <- obj$beta <- list()
	drop <- c()
	N <- length(Y)
	keep <- c(1:N)
	# work in stages starting from final stage (K)
	for (j in K:1) {
		# store current outcome for use if missingness occurs
		Y.orig <- Y
		# get data: Hpsi = blip, Hbeta = treatment-free, Halpha = treatment
		A <- model.response(model.frame(treat.mod[[j]],data, na.action='na.pass'))
		Hpsi <- model.matrix(blip.mod[[j]],model.frame(blip.mod[[j]],data, na.action='na.pass'))
		Hbeta <- model.matrix(tf.mod[[j]],model.frame(tf.mod[[j]],data, na.action='na.pass'))
		Halpha <- model.matrix(treat.mod[[j]],model.frame(treat.mod[[j]],data, na.action='na.pass'))
		# list of treatment variables
		treat.var <- all.vars(treat.mod[[j]])[1]
		# store blip variable names for output
		obj$blip.list[[j]] <- colnames(Hpsi)
		# missing data: identify if missing in Hpsi, Hbeta, Halpha or A; drop and store sample size
		if (missing == "default" | missing == "ipcw") {
			# pi = vector weights for IPCW
			pi <- rep(1,N)
			# drop adds to current drop so future missingness is accommodated
			new.drop <- which(apply(is.na(cbind(A,Hpsi,Hbeta,Halpha)),1,any))
			drop <- unique(new.drop,drop)
			# if missing, later than stage 1, and requested, calculate IPCW
			if (length(new.drop) > 0 & missing == "ipcw" & j > 1) {
				# this works on 'keep' from previous stage: full population who *could* be missing now
				Miss <- rep(1,N)
				Miss[new.drop] <- 0
				Miss <- Miss[keep]
				# build full history from all previous stages
				H.list <- unique(c(unlist(sapply(blip.mod[(j-1):1],all.vars)),unlist(sapply(tf.mod[(j-1):1],all.vars)),unlist(sapply(treat.mod[(j-1):1],all.vars))))
				H <- as.matrix(data[,which(colnames(data) %in% H.list)])
				H <- H[keep,]
				pi[which(!apply(is.na(H),1,any))] <- fitted(glm(Miss~H))
				obj$ipcw[[j]] <- pi
			}
			keep <- keep[!(keep %in% drop)]
			pi <- pi[keep]
			A <- A[keep]
			Hpsi <- Hpsi[keep,]
			Hbeta <- Hbeta[keep,]
			Halpha <- Halpha[keep,]
			Y <- Y[keep]
		}
		# force matrices in case intercept-only models
		Hpsi <- as.matrix(Hpsi)
		Hbeta <- as.matrix(Hbeta)
		Halpha <- as.matrix(Halpha)
		obj$n[[j]] <- length(keep)
		n <- length(Y)
		# check if treatments are binary, if not, and dWOLS selected, abort
		A.bin <- as.numeric(length(unique(A)) <= 2)
		if (A.bin == 0 & method %in% c("dwols","dWOLS")) {
			stop("Non-binary treatment only suitable for G-estimation analysis.\n")
		}
		# if binary and not 0/1, recode (and print warning)
		if (A.bin & !all(unique(A) %in% c(0,1))) {
			A <- ifelse(A==min(A),0,1)
			# suppress output if bootstrapping
			if (var.estim!="bs.quiet") {
				cat("Warning: stage",j,"treatment recoded to 0/1.\n")
			}
		}
		# if no treat.range but outcome non-binary, set treat range to min/max of A
		if (A.bin == 0 & is.null(treat.range)) {
			treat.range <- c(min(A),max(A))
			if (var.estim!="bs.quiet") {
				cat("Warning: no treatment range specified, [min,max] of observed treatment used.\n")
			}
		}
		obj$treat.range[[j]] <- treat.range
		# treatment model: if binary logistic, if non-binary linear
		if (A.bin == 1) {
			obj$cts[[j]] <- "bin"
			# note that if user-specified we take specified P(A = 1) instead
			if (class(treat.mod.man) == "list") {
				if (class(treat.mod.man[[j]]) != "formula") {
					Ahat <- treat.mod.man[[j]]
				}
			} else {
				alpha <- glm(A~Halpha,binomial)
				obj$treat.mod.fitted[[j]] <- alpha
				Ahat <- fitted(alpha)
			}
		} else {
			alpha <- lm(treat.mod[[j]],data)
			obj$treat.mod.fitted[[j]] <- alpha
			Ahat <- fitted(alpha)
			if (treat.var %in% colnames(Hpsi)){
				# continuous with quadratic blip
				obj$cts[[j]] <- "cts.q"
			} else {
				# continuous with linear blip
				obj$cts[[j]] <- "cts.l"
			}
		}
		# if dWOLS
		if (method == "dwols") {
			# weights
			if (is.function(weight) == FALSE) {
				w <- abs(A - Ahat)
			} else {
				# weight function specified for pi = 1, therefore
				w <- weight(Ahat)
				w[A == 0] <- w[A==0]*Ahat[A==0]/(1-Ahat[A==0])
			}
			# if q-learning selected, set all weights to 1
			if (qlearn == 1) {w <- rep(1,n)}
			# update with IPCW
			w <- w*pi
			# blip parameter estimates extracted via dimensions of Hbeta and Hpsi
			est <- solve(t(cbind(Hbeta,A*Hpsi))%*%cbind(w*Hbeta,w*A*Hpsi)) %*% t(cbind(Hbeta,A*Hpsi)) %*% (w*Y)
			psi <- est[(dim(Hbeta)[2]+1):(dim(Hbeta)[2]+dim(Hpsi)[2])]
			beta <- est[1:dim(Hbeta)[2]]
		}
		# g-estimation
		if (method == "gest") {
			Hd <- cbind(Hbeta,A*Hpsi)
			if (A.bin == 1) {
				# binary treatment
				Hw <- cbind(Hbeta,Hpsi*(A - Ahat))
			} else {
				# continuous treatment, requires some fiddly cleanup
				# ***INPUT CLEANUP STARTS
				# separate Hpsi into A/not-A components (if applicable)
				p1.list <- p2.list <- c()
				for (p.var in colnames(Hpsi)) {
					# note there may be interaction terms, so have to separate on :
					if (treat.var %in% strsplit(p.var,split=":")[[1]]) {
					p2.list <- append(p2.list,which(colnames(Hpsi) == p.var))} 						else {
					p1.list <- append(p1.list,which(colnames(Hpsi) == p.var))
					}
				}
				Hpsi1 <- as.matrix(Hpsi[,p1.list])
				# if we have A terms need to extract what they're multiplied with
				# 3 options: A by itself, or an A interaction as :A or A:
				psi2.list <- c()
				for (p.var in colnames(Hpsi)[p2.list]) {
					# if p.var is treat.var, we do nothing (assumed intercept term)
					if (p.var != treat.var) {
						# if variable starts with A: then replace with nothing
						if (strsplit(p.var,split=":")[[1]][1] == treat.var) {psi2.list <- append(psi2.list,paste(strsplit(p.var,split=":")[[1]][-1],collapse=":"))}
						# if variable ends with :A then replace with nothing
						if (tail(strsplit(p.var,split=":")[[1]],n=1) == treat.var) {psi2.list <- append(psi2.list,paste(head(strsplit(p.var,split=":")[[1]],n=-1),collapse=":"))}
						# if variable has :A: in the middle, replace this with :
						if (strsplit(p.var,split=paste(":",treat.var,":",sep=""))[[1]][1] != p.var) {psi2.list <- append(psi2.list,paste(strsplit(p.var,split=paste(":",treat.var,":",sep=""))[[1]],collapse=":"))}
					}
				}
				# turn psi2.list into a model so we can get a new data frame
				# intercept always implied
				if (length(psi2.list) > 0) {psi2.mod <- paste("~1+",paste(psi2.list,collapse="+")); obj$b.list[[j]] <- psi2.list} else {psi2.mod <- "~1"}
				if (length(p1.list) > 1) {psi1.mod <- paste("~1+",paste(colnames(Hpsi)[p1.list[-1]],collapse="+"))} else {psi1.mod <- "~1"}
				# store both components of the blip model for prediction function
				obj$blip.list.cts[[j]] <- c(psi1.mod,psi2.mod)
				Hpsi2 <- model.matrix(as.formula(psi2.mod),data)
				# ***INPUT CLEANUP ENDS
				A2 <- A^2
				A2hat <- (summary(alpha)$sigma)^2 + Ahat^2
				# if no A^2 terms, just use Hpsi
				ifelse(length(p2.list) > 0,Hw <- cbind(Hbeta,Hpsi1*(A-Ahat),Hpsi2*(A2-A2hat)),Hw <- cbind(Hbeta,Hpsi*(A-Ahat)))
			}


######

# w-matrix for analysis
#W3 <- (diag(A3-A3hat) - (A3-A3hat)*H3.beta %*% solve(t(H3.beta) %*% H3.beta) %*% t(H3.beta))
# estimate parameters
#psi3 <- solve(t(H3.psi) %*% W3 %*% (A3*H3.psi)) %*% t(H3.psi) %*% W3 %*% Y3
#		W <- (diag(A-Ahat) - (A-Ahat)*Hbeta %*% solve(t(Hbeta) %*% Hbeta) %*% t(Hbeta))
#		psi <- solve(t(Hpsi) %*% W %*% (A*Hpsi)) %*% t(Hpsi) %*% W %*% (Y*pi)
#		beta <- solve(t(Hbeta) %*% Hbeta) %*% t(Hbeta) %*% (Y*pi - A*Hpsi%*%psi)
######
			# get estimates
			est <- solve(t(Hw)%*%Hd) %*% t(Hw) %*% (Y*pi)
			# blip parameters
			psi <- est[(dim(Hbeta)[2]+1):(dim(Hbeta)[2]+dim(Hpsi)[2])]
			# non-blip parameters
			beta <- est[1:dim(Hbeta)[2]]
		}
		# use estimates to identify optimal treatments
		# if treatment is binary:
		if (A.bin == 1) {opt <- as.numeric(Hpsi %*% psi > 0)}
		# if treatment is continuous
		if (A.bin == 0) {
			if (obj$cts[[j]] == "cts.l") {
				# if no treatment term within blip, then max/min based on sign
				opt <- as.numeric(Hpsi %*% psi > 0)*max(treat.range) + as.numeric(Hpsi %*% psi <= 0)*min(treat.range)
			} else {
				# otherwise optimal is found by maximizing wrt A, need to check second derivative
				# if sign(Hpsi2 %*% psi2) < 0, then blip is maxed by setting first derivative to zero
				# otherwise evaluate at max/min or treatment range and pick the larger
				# dif is blip(minA) - blip(maxA); so if negative take max, otherwise min
				sec.deriv <- Hpsi2 %*% psi[(dim(Hpsi1)[2]+1):length(psi)]
				dif <- sign(min(treat.range)*(Hpsi1 %*% psi[1:(dim(Hpsi1)[2])] + min(treat.range)*Hpsi2 %*% psi[(dim(Hpsi1)[2]+1):length(psi)]) - max(treat.range)*(Hpsi1 %*% psi[1:(dim(Hpsi1)[2])] + max(treat.range)*Hpsi2 %*% psi[(dim(Hpsi1)[2]+1):length(psi)]))
				opt <- as.numeric(sec.deriv < 0)*(-(Hpsi1 %*% psi[1:(dim(Hpsi1)[2])])/(2*Hpsi2 %*% psi[(dim(Hpsi1)[2]+1):length(psi)])) + as.numeric(sec.deriv >= 0)*(as.numeric(dif >= 0)*min(treat.range) + as.numeric(dif < 0)*max(treat.range))
			}
			# make sure opt is between limits
			opt <- pmin(pmax(opt,min(treat.range)),max(treat.range))
		}
		# ensure opt is in vector form
		opt <- as.vector(opt)
		# store estimates, optimal treatments and blip parameters
		obj$beta[[j]] <- beta
		obj$psi[[j]] <- psi
		obj$opt.treat[[j]] <- opt
		# Sandwich variance estimator done stage-by-stage
		if (var.estim == "sandwich") {
			D <- A - Ahat
			# HD: what precedes the (Y-gamma-beta) term in the est eq.
			if (method == "gest") {HD <- Hpsi*(A-Ahat)}
			if (method == "dwols") {HD <- A*Hpsi*abs(A-Ahat)}
			# residuals are the (Y-gamma-beta) residuals
			residuals <- as.vector(Y - ((A*Hpsi)%*%psi) - Hbeta%*%beta)
			# estimating equation = (Y - gamma - E[Y - gamma|H])*(A - Ahat)*Hpsi
			U <- residuals * HD
			# dU/dpsi:
			M <- (1/n) * -crossprod(HD, A * Hpsi)
			# if binary
			if (A.bin == 1) {
				score.alpha <- D * Halpha
				score.varsig <- residuals * Hbeta
				U.adj <- U
				# dT.dalpha = Ahat(1-Ahat) Halpha Halpha
				dT.dalpha <- (1/n)*(-crossprod(Ahat*(1-Ahat)*Halpha, Halpha))
				dU.dalpha <- (1/n)*(-crossprod(Ahat*(1-Ahat)*residuals*Hpsi,Halpha))
				U.adj <- U.adj - t(dU.dalpha%*%solve(dT.dalpha) %*% t(score.alpha))
				# now with respect to beta
				dT.dvarsig <- (1/n) * (-crossprod(Hbeta))
				dU.dvarsig <- (1/n) * (-crossprod(HD, Hbeta))
				U.adj <- U.adj-t(dU.dvarsig%*%solve(dT.dvarsig)%*%t(score.varsig))
				if (j == K) {U.store <- M.store <- deriv.regret.store <- list()}
				U.store[[j]] <- matrix(nrow=N,ncol=dim(U)[2])
				deriv.regret.store[[j]] <- matrix(nrow=N,ncol=dim(Hpsi)[2])
				U.store[[j]][keep,] <- U
				M.store[[j]] <- M
				# derivative of regrets wrt psi
				deriv.regret.store[[j]][keep,] <- (opt - A) * Hpsi
				if (j < K) {
					for (l in (j+1):K) {
						dU2.dpsi2 <- -M.store[[l]]
						dU1.dpsi2<-(1/n)*(crossprod(HD,deriv.regret.store[[l]][keep,]))
						U.adj <- U.adj - t(dU1.dpsi2 %*% solve(dU2.dpsi2) %*% t(U.store[[l]][keep,]))
					}
				}
			}
			# if linear, then score = (1/s^2) * (A - Ahat)*Halpha 
			if (A.bin == 0) {
				score.alpha <- D * Halpha
				score.varsig <- residuals * Hbeta
				U.adj <- U
				# dT.dalpha = -1/s^2 t(Halpha) * Halpha
				dT.dalpha <- (1/n) * (-crossprod(Halpha,Halpha))
				# dU/dalpha = (y - gamma - G)*Hpsi*[...], if no A^2 term, [...] is -Halpha
				# A^2 term then also -2*Halpha*alpha-hat - 2/(n-2)*sum(Halpha*(Ahat-A))
				# N.B. this slight ugliness is because of the var(A-hat|...) term
				if (obj$cts[[j]] == "cts.l") {
					dU.dalpha <- (1/n) * (-crossprod(residuals * Hpsi, Halpha))
				} else {
					Hpsi.cts <- cbind(Hpsi1,Hpsi2*2*(Ahat + (1/(dim(Hbeta)[1]-2))*sum(Ahat-A)))
					dU.dalpha <- (1/n) * (-crossprod(residuals * Hpsi.cts,Halpha))
				}
				U.adj <- U.adj - t(dU.dalpha %*% solve(dT.dalpha) %*% t(score.alpha))#
				dT.dvarsig <- (1/n) * (-crossprod(Hbeta))
				dU.dvarsig <- (1/n) * (-crossprod(HD, Hbeta))
				U.adj <- U.adj - t(dU.dvarsig %*% solve(dT.dvarsig) %*% t(score.varsig))
				# storage for earlier intervals
				if (j == K) {U.store <- M.store <- deriv.regret.store <- list()}
				U.store[[j]] <- U
				M.store[[j]] <- M
				# derivative of regrets wrt psi, depends on whether A term in reg
				if (obj$cts[[j]] == "cts.l") {deriv.regret.store[[j]] <- (opt - A) * Hpsi}
				if (obj$cts[[j]] == "cts.q") {deriv.regret.store[[j]] <- opt*(cbind(Hpsi1,opt*Hpsi2))-A*Hpsi}
				if (j < K) {
					for (l in (j+1):K) {
		    				dU2.dpsi2 <- -M.store[[l]]
						dU1.dpsi2 <- (1/n) * (crossprod(HD, deriv.regret.store[[l]]))
						U.adj <- U.adj - t(dU1.dpsi2 %*% solve(dU2.dpsi2) %*% t(U.store[[l]]))
					}
				}
			}
			IF <- t(solve(M, t(U.adj)))
			Sigma <- (1/n) * var(IF)
			obj$covmat[[j]] <- Sigma
		}
		# if type = "DTR": update pseudo-Y by adding regret, if SNMM (i.e. not DTR), subtract blip, note that if missing we take Y = Y.orig
		if (obj$cts[[j]] %in% c("cts.l","bin")) {
			# for linear blip
			Y.orig[keep] <- Y.orig[keep] + (as.numeric(type == "DTR")*(opt) - A)*(Hpsi%*%psi)
			Y.orig[drop] <- Y.orig[drop]
			obj$regret[[j]] <- (opt - A)*(Hpsi%*%psi)
		} else {
			# for quadratic blip
			Y.orig[keep] <- Y.orig[keep] + as.numeric(type == "DTR")*(opt)*(Hpsi1 %*% psi[1:length(p1.list)] + as.numeric(type == "DTR")*(opt)*Hpsi2 %*% psi[(length(p1.list)+1):(length(psi))]) - A*(Hpsi%*%psi)
			Y.orig[drop] <- Y.orig[drop]
			obj$regret[[j]] <- opt*(Hpsi1 %*% psi[1:length(p1.list)] + opt*Hpsi2 %*% psi[(length(p1.list)+1):(length(psi))]) - A*(Hpsi%*%psi)
		}
		Y <- Y.orig
		# fitted values at each stage
		Y.fit <- Hbeta %*% beta + A*(Hpsi %*% psi)
		k <- j
		while (k < K) {
			k <- k + 1
			Y.fit <- Y.fit - obj$regret[[k]][keep]
		}
		obj$fitted[[j]] <- Y.fit
	}
	# optimal outcome according to blip estimates is current pseudo-Y
	obj$opt.Y <- Y
	obj$type <- type
	# standard errors
	if (var.estim=="bootstrap") {
		if (M == 0) {M <- nrow(data)}
		# note time for progress reports
		ptm <- proc.time()
		psi.boot <- list()
		# continue indicator
		cont <- "u"
		for (i in 1:B) {
			data.boot <- data[sample(1:M, replace=TRUE),]
			psi.boot[[i]] <- DTRreg(outcome, blip.mod, treat.mod, tf.mod, data.boot, method, var.estim="bs.quiet", verbose=verbose,treat.mod.man=treat.mod.man)$psi
			# if verbose, display ETA and give option to abort
			if (verbose == T & i >= 10) {
				# only display if projected > 30s, then only display every 30s
				eta <- (B-i)*((proc.time() - ptm)[3]/i)
				# if very long (> 10 mins) give option to abort but if ignored continue
				if (i > 50 & eta > 600 & cont != "y" & interrupt == T) {
					cont <- readline(paste("Estimated run time",round(eta/60),"minutes. Continue? y/n: "))
					if (cont == "n" | cont == "no" | cont == "NO") {stop("Aborted.\n")} else {cont <- "y"}
				}
				if (i == 10) {last <- eta+31}
				if (eta > 30 & eta < (last-30)) {
					cat("Approximately",round(eta),"seconds remaining.\n")
					last <- eta
				}
			}
		}
		psi.boot <- do.call(function(...) mapply(rbind, ..., SIMPLIFY=FALSE), psi.boot)
		# if requested, truncate estimates based on specified percentile
		if (truncate > 0) {
			for (j in 1:K) {
				psi.boot[[j]] <- apply(psi.boot[[j]],2,DTRreg.trunc,truncate)
			}
		}
		covmat <- lapply(psi.boot, var)
		obj$covmat <- covmat
		obj$psi.boot <- psi.boot
	}
	# if variance estimated (and we're not currently bootstrapping, estimate non-regularity
	if (!(var.estim %in% c("none","bs.quiet")) & A.bin == 1) {
		for (j in K:1) {
			Hpsi <- model.matrix(blip.mod[[j]],data)
			# lower/upper limits on psi
			psi.l <- obj$psi[[j]] - 1.96*sqrt(diag(obj$covmat[[j]]))
			psi.u <- obj$psi[[j]] + 1.96*sqrt(diag(obj$covmat[[j]]))
			# look at max/min value of blip based on sign of covariates and max/min of parameter CIs
			Hpsi.p <- apply(Hpsi,c(1,2),FUN=function(x) max(x,0))
			Hpsi.n <- apply(Hpsi,c(1,2),FUN=function(x) min(x,0))
			# lower blip is positive values * lower CI + negative values * upper CI
			blip.l <- Hpsi.p %*% psi.l + Hpsi.n %*% psi.u
			blip.u <- Hpsi.p %*% psi.u + Hpsi.n %*% psi.l
			obj$nonreg[j] <- sum(blip.l < 0 & blip.u > 0)/dim(Hpsi)[1]
		}
	}
	# return
	class(obj) <- "DTRreg"
	obj
}

# printing
#' @export
print.DTRreg <- function(x,...) {
	K <- x$K
	# longest blip variable name lengths, if applicable, get length of blip excluding treatment terms
	max.length <- 0
	p.length <- c()
	for (j in 1:K) {
		max.length <- max(max.length,nchar(x$blip.list[[j]]))
		if (j <= length(x$blip.list.cts)) {if (is.null(x$blip.list.cts[[j]]) == F) {p.length[j] <- length(strsplit(x$blip.list.cts[[j]][1],split="\\+")[[1]])} else {p.length[j] <- 0}}
	}
	if (x$type == "DTR") {cat("DTR estimation over",K,"stages:\n\n")}
	cat("Blip parameter estimates")
	if (is.null(x$covmat[1])) {cat(" (standard errors not requested)")}
	cat("\n")
	# first blank space should be the length of longest variable + 1
	blank.paste <- function(n) {if (n > 0) {for (i in 1:n) {cat(" ")}}}
	blank.paste(max.length)
	# if no standard errors, just print estimates
	if (is.null(x$covmat[1][[1]])) {
		cat(" Estimate\n")
		for (j in 1:K) {
			cat("Stage ",j," (n = ",x$n[[j]],")\n",sep="")
			for (p in 1:length(x$psi[[j]])) {
				blank.paste(max.length-nchar(x$blip.list[[j]][p]))
				cat(x$blip.list[[j]][p])
				blank.paste(3 - as.numeric(x$psi[[j]][p] < 0))
				cat(paste(sprintf("%.4f",x$psi[[j]][p]),"\n",sep=""))
			}
			cat("\n")
		}
	}
	# otherwise, print standard errors and CIs
	if (is.null(x$covmat[1][[1]]) == F) {
		cat(" Estimate Std. Error    95% Conf. Int\n")
		for (j in 1:K) {
			cat("Stage ",j," (n = ",x$n[[j]],")\n",sep="")
			for (p in 1:length(x$psi[[j]])) {
				blank.paste(max.length-nchar(x$blip.list[[j]][p]))
				cat(x$blip.list[[j]][p])
				blank.paste(3 - as.numeric(x$psi[[j]][p] < 0))
				cat(paste(sprintf("%.4f",x$psi[[j]][p]),"     ",sprintf("%.4f",sqrt(x$covmat[[j]][p,p])),sep=""))
				blank.paste(2 - as.numeric(x$psi[[j]][p] < 0))
				cat(paste("[",sprintf("%.4f",x$psi[[j]][p]-1.96*sqrt(x$covmat[[j]][p,p])),",",sprintf("%.4f",x$psi[[j]][p]+1.96*sqrt(x$covmat[[j]][p,p])),"]\n",sep=""))
			}
			cat("\n")
		}
	}
	# print warnings if non-regularity concerns
	if (x$type == "DTR") {
	if (is.null(x$covmat[1][[1]]) == F) {
		for (j in 1:K) {
			if (x$cts[[j]] == "bin") {
			if (x$nonreg[j] > 0.05) {cat("Warning: possible non-regularity at stage ",j," (prop = ",round(x$nonreg[j],3),")\n",sep="")}
		}
		}
		cat("\n")
	}
	# print DTRregs explicitly
	if (!("cts.q" %in% x$cts)) {cat("Recommended dynamic treatment regimen:\n")}
	# if binary just treat/don't treat
	for (j in 1:K) {
		if (x$cts[[j]] == "bin") {
			cat("Stage ",j,": treat if ",sprintf("%.4f",x$psi[[j]][1]),sep="")
			if (length(x$psi[[j]]) > 1) {
				for (p in 2:length(x$psi[[j]])) {
					if(x$psi[[j]][p] >= 0) {cat(" + ")} else {cat(" - ")}
					cat(sprintf("%.4f",abs(x$psi[[j]][p])),x$blip.list[[j]][p])
				}
			}
		cat(" > 0\n")
		} else {
			# set to max/min of range depending on sign of blip
			if (x$cts[[j]] == "cts.l") {
				cat("Stage ",j,": maximum treatment if ",sprintf("%.4f",x$psi[[j]][1]),sep="")
				if (length(x$psi[[j]]) > 1) {
					for (p in 2:length(x$psi[[j]])) {
						if(x$psi[[j]][p] >= 0) {cat(" + ")} else {cat(" - ")}
						cat(sprintf("%.4f",abs(x$psi[[j]][p])),x$blip.list[[j]][p])
					}
				}
			cat(" > 0, otherwise minimum treatment.\n")
			}
		}
	}
	}
}
#' @export
summary.DTRreg <- function(object,...) {
	print.DTRreg(object)
}

# optimal outcome prediction for new data ASSUMING CORRECT MODELS
#' @export
predict.DTRreg <- function(object, newdata, treat.range=NULL, ...) {
	x <- object
	Hpsi <- model.matrix(x$blip.mod[[1]],newdata)
	Hbeta <- model.matrix(x$tf.mod[[1]],newdata)
	if (is.null(treat.range) & x$cts[[1]] != "bin") {treat.range <- x$treat.range[[1]]}
	# options: binary, cts with range, or cts without range
	# binary:
	if (x$cts[[1]] == "bin") {
		opt <- as.numeric(Hpsi%*%x$psi[[1]] > 0)
		obj <- Hbeta%*%x$beta[[1]] + opt*(Hpsi%*%x$psi[[1]])
	}
	# cts with no A^2 term in blip
	if (x$cts[[1]] == "cts.l") {
		opt <- as.numeric(Hpsi%*%x$psi[[1]] > 0)*max(treat.range) + as.numeric(Hpsi%*%x$psi[[1]] <= 0)*min(treat.range)
		obj <- Hbeta%*%x$beta[[1]] + opt*(Hpsi%*%x$psi[[1]])
	}
	# cts with A^2 term in blip
	if (x$cts[[1]] == "cts.q") {
		Hpsi1 <- model.matrix(as.formula(x$blip.list.cts[[1]][1]), newdata)
		Hpsi2 <- model.matrix(as.formula(x$blip.list.cts[[1]][2]), newdata)
		psi <- x$psi[[1]]
		sec.deriv <- Hpsi2 %*% psi[(dim(Hpsi1)[2]+1):length(psi)]

		dif <- sign(min(treat.range)*(Hpsi1 %*% psi[1:(dim(Hpsi1)[2])] + min(treat.range)*Hpsi2 %*% psi[(dim(Hpsi1)[2]+1):length(psi)]) - max(treat.range)*(Hpsi1 %*% psi[1:(dim(Hpsi1)[2])] + max(treat.range)*Hpsi2 %*% psi[(dim(Hpsi1)[2]+1):length(psi)]))
		opt <- as.numeric(sec.deriv >= 0)*(-(Hpsi1 %*% psi[1:(dim(Hpsi1)[2])])/(2*Hpsi2 %*% psi[(dim(Hpsi1)[2]+1):length(psi)])) + as.numeric(sec.deriv < 0)*(as.numeric(dif >= 0)*min(treat.range) + as.numeric(dif > 0)*max(treat.range))

		opt <- pmin(pmax(opt,min(treat.range)),max(treat.range))
		obj <- Hbeta%*%x$beta[[1]] + opt*(Hpsi1%*%psi[1:dim(Hpsi1)[2]] + Hpsi2%*%psi[(dim(Hpsi1)[2]+1):length(psi)])
	}
	return(obj)
}
#' @export
coef.DTRreg <- function(object,...) {
	x <- object
	K <- x$K
	obj <- x$psi
	for (j in 1:K) {
		names(obj)[j] <- paste("stage",j,sep="")
		names(obj[[j]]) <- x$blip.list[[j]]
	}
	return(obj)
}
# confidence interval
#' @export
confint.DTRreg <- function(object, parm = NULL, level = 0.95, ...) {
	x <- object
	K <- x$K
	obj <- x$psi
	if (is.null(x$covmat[1][[1]]) == F) {
		for (j in 1:K) {
			names(obj)[j] <- paste("stage",j,sep="")
			obj[[j]] <- matrix(nrow=length(x$psi[[j]]),ncol=2)
			rownames(obj[[j]]) <- x$blip.list[[j]]
			colnames(obj[[j]]) <- c(paste(100*((1 - level)/2),"%"), paste(100*((1 + level)/2),"%"))
			for (p in 1:length(x$psi[[j]])) {
				obj[[j]][p,] <- c(x$psi[[j]][p] + qnorm((1-level)/2) * sqrt(x$covmat[[j]][p,p]), x$psi[[j]][p] + qnorm((1+level)/2)*sqrt(x$covmat[[j]][p,p]))
			}
		}
	return(obj)
	} else {
		cat("confint() only available when DTRreg called with a variance estimation option.\n")
	}
}
# diagnostics
#' @export
plot.DTRreg <- function(x, ...) {
	# original Y is in x$obs.Y
	# fitted Y (one for each stage) is in x$fitted
	cat("Y residuals (for treatment-free and blip models):\n")
	for (j in 1:x$K) {
		advance <- readline("Hit <Return> to see next plot:")
		plot(x$fitted[[j]],x$obs.Y - x$fitted[[j]],xlab="Fitted values",ylab="Residuals",main=paste("Stage",j,"Residuals vs Fitted"))
		abline(h=0,lty=2,col=4)
		points(lowess(x$fitted[[j]],x$obs.Y - x$fitted[[j]]),type="l",col=2,lwd=2)
		if (length(x$blip.list[[j]]) > 1) {
			for (var in x$blip.list[[j]][2:length(x$blip.list[[j]])]) {
				advance <- readline("Hit <Return> to see next plot:")
				plot(x$data[,which(colnames(x$data) == var)],x$obs.Y-x$fitted[[j]],xlab=paste(var),ylab="Residuals",main=paste("Stage",j,"Residuals vs",var))
				abline(h=0,lty=2,col=4)
				points(lowess(x$data[,which(colnames(x$data) == var)],x$obs.Y-x$fitted[[j]]),type="l",col=2,lwd=2)
			}
		}
	}
	cat("Treatment models:\n")
	for (j in 1:x$K) {
		cat("Stage",j,"\n")
		plot(x$treat.mod.fitted[[j]],sub.caption=paste("Stage",j,"treatment"))
	}
}

# truncation function for truncate option
DTRreg.trunc <- function(x,p) {
	l <- quantile(x,p); u <- quantile(x,1-p)
	x[x < l] <- l; x[x > u] <- u
	x
}

# get_all_vars modified to accommodate intercept-only components
get_all_vars_DTR <- function(x,n) {
	if (dim(get_all_vars(x))[1] == 0) {
		return(rep(1,n))
	} else {
		return(get_all_vars(x))
	}
}



