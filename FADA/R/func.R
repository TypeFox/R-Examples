bivprob = function(rho, lower, upper = -lower, mean = 0) {
	nu = 0
	low = rep(as.double((lower - mean)), 2)
	upp = rep(as.double((upper - mean)), 2)
	if (any(lower == upper)) 
		return(0)
	infin = c(2, 2)
	infin = as.integer(infin)
	low = replace(low, low == -Inf, 0)
	upp = replace(upp, upp == Inf, 0)
	rho = as.double(rho)
	prob = as.double(0)
	a = lapply(rho, function(r, low, upp) biv.nt.prob(df = Inf, lower = low, upper = upp, 
		mean = rep(0, 2), S = matrix(c(1, r, r, 1), 2, 2)), low = low, upp = upp)
	return(unlist(a))
}
Dt = function(rho) {
	threshold = 0.05
	ut = qnorm(1 - threshold/2)
	delta = unlist(lapply(rho, bivprob, lower = -ut)) - (1 - threshold)^2
	dt <- delta/(threshold * (1 - threshold))
	return(dt)
}
VarInflation <- function(data.train, Blist, maxnbfactors, dig) {
	m <- ncol(data.train)
	n <- nrow(data.train)
	vecrho <- round(seq(10^(-dig), 1, 10^(-dig)), digits = dig)
	vecdt <- unlist(lapply(vecrho, Dt))
	sampled <- sample(1:m, min(1000, m))
	sampsize <- length(sampled)
	cordata <- crossprod(data.train[, sampled, drop = FALSE])/(n - 1)
	sdt <- sapply(1:(maxnbfactors + 1), function(i) {
		B <- matrix(Blist[[i]][sampled, ], nrow = sampsize)
		sdb <- sqrt(1 - rowSums(B^2))
		matrho <- cordata - tcrossprod(B)
		matrho <- sweep(matrho, 2, FUN = "/", STATS = sdb)
		matrho <- sweep(matrho, 1, FUN = "/", STATS = sdb)
		rho <- matrho[col(matrho) > row(matrho)]
		rho[abs(rho) >= 1] <- 1
		veccor <- sort(round(abs(rho), digits = dig))
		duplic <- duplicated(veccor)
		vduplic <- sort(unique(veccor[duplic]))
		vunic <- setdiff(unique(veccor), vduplic)
		dtunic <- vecdt[is.element(vecrho, vunic)]
		dtduplic <- vecdt[is.element(vecrho, vduplic)]
		vmatch <- match(vecrho, veccor, 0)
		nboccur <- diff(c(vmatch[vmatch > 0], length(veccor) + 1))
		nboccur <- nboccur[nboccur > 1]
		tmp <- 2 * (m - 1) * (sum(dtunic) + crossprod(nboccur, dtduplic))/(sampsize * 
			(sampsize - 1))
		return(tmp)
	})
	names(sdt) <- paste(0:maxnbfactors, "factors")
	return(sdt)
}
ifa = function(Psi, B) {
	if (class(B) == "numeric") 
		B = matrix(B, ncol = 1)
	q = ncol(B)
	Phi = rep(0, length(Psi))
	Phi[abs(Psi) > 1e-05] = 1/Psi[abs(Psi) > 1e-05]
	PhiB = tcrossprod(Phi, rep(1, q))
	PhiB = PhiB * B
	G = diag(q) + t(B) %*% PhiB
	GinvtPhiB = tcrossprod(solve(G), PhiB)
	Phib2 = tcrossprod(PhiB, t(GinvtPhiB))
	iS = diag(Phi) - Phib2
	PhiB2 = crossprod(PhiB, B)
	GinvtPhiB2 = crossprod(solve(G), PhiB2)
	Phib2 = tcrossprod(PhiB, t(GinvtPhiB2))
	iSB = PhiB - Phib2
	return(list(iS = iS, iSB = iSB))
}
emfa = function(data, nbf, EM = TRUE, minerr = 1e-06, verbose = FALSE) {
	n = nrow(data)
	m = ncol(data)
	my = crossprod(rep(1, n), data)/n
	vy = crossprod(rep(1, n), data^2)/n - my^2
	vy = (n/(n - 1)) * vy
	cdata = scale(data, center = my, scale = FALSE)
	csdata = scale(data, center = my, scale = sqrt(vy))
	S = crossprod(csdata)/(n - 1)
	if (((n > m) & (m <= 200) & (m >= 3)) & (!EM)) {
		if (nbf == 0) {
			B = NULL
			Psi = rep(1, m)
		}
		if (nbf > 0) {
			fa = factanal(csdata, factors = nbf, rotation = "varimax")
			B = fa$loadings
			class(B) = "matrix"
			Psi = fa$uniquenesses
			Psi = Psi * vy
			B = matrix(rep(sqrt(vy), ncol(B)), nrow = nrow(B)) * B
			sB = scale(t(B), center = FALSE, scale = sqrt(Psi))
			G = solve(diag(nbf) + tcrossprod(sB))
			sB = scale(t(B), center = FALSE, scale = Psi)
		}
	}
	if ((n <= m) | (m > 200) | EM) {
		if (nbf == 0) {
			B = NULL
			Psi = rep(1, m)
		}
		if (nbf > 0) {
			if (verbose) 
				print(paste("Fitting EM Factor Analysis Model with", nbf, "factors"))
			eig = fast.svd((1/sqrt((n - 1))) * t(csdata))
			evectors = eig$u[, 1:nbf]
			evalues = eig$d^2
			if (nbf > 1) 
				B = evectors[, 1:nbf] * matrix(sqrt(evalues[1:nbf]), ncol = nbf, 
					nrow = m, byrow = TRUE)
			if (nbf == 1) 
				B = matrix(evectors, nrow = m, ncol = 1) * sqrt(evalues[1])
			b2 = rowSums(B^2)
			Psi = 1 - b2
			crit = 1
			while (crit > minerr) {
				inv = ifa(Psi, B)
				Cyz = crossprod(S, inv$iSB)
				Czz = crossprod(inv$iSB, Cyz) + diag(nbf) - crossprod(B, inv$iSB)
				Bnew = tcrossprod(Cyz, solve(Czz))
				Psinew = 1 - rowSums(Bnew * Cyz)
				crit = mean((Psi - Psinew)^2)
				B = Bnew
				Psi = Psinew
				if (verbose) 
					print(paste("Objective criterion in EM-FA : ", signif(crit, 6)))
			}
			Psi = Psi * vy
			B = matrix(rep(sqrt(vy), ncol(B)), nrow = nrow(B)) * B
			sB = scale(t(B), center = FALSE, scale = sqrt(Psi))
			G = solve(diag(nbf) + tcrossprod(sB))
			sB = scale(t(B), center = FALSE, scale = Psi)
		}
	}
	res = list(B = B, Psi = Psi)
	return(res)
}
nbfactors <- function(data.train, maxnbfactors = 12, diagnostic.plot, minerr = 0.001, 
	EM = TRUE, jumps.nbfactor = 0.05) {
	dig <- 2
	m <- ncol(data.train)
	n <- nrow(data.train)

	my = crossprod(rep(1, n), data.train)/n
	vy = crossprod(rep(1, n), data.train^2)/n - my^2
	vy = (n/(n - 1)) * vy
	cdata = scale(data.train, center = my, scale = FALSE)
	csdata = scale(data.train, center = my, scale = sqrt(vy))
	S = crossprod(csdata)/(n - 1)
	eig = fast.svd((1/sqrt((n - 1))) * t(csdata))

	falist <- vector(length = maxnbfactors + 1, "list")
	falist[[1]] <- list(B = matrix(0, ncol = 1, nrow = m))
	falist[-1] <- lapply(1:maxnbfactors, emfa.nbf, csdata = csdata, S = S, eig = eig, 
		vy = vy, minerr = minerr, EM = EM, verbose = FALSE)
	Blist <- lapply(falist, function(fa, m) matrix(fa$B, nrow = m), m = m)
	sdt <- VarInflation(data.train, Blist, maxnbfactors, dig)
	if (diagnostic.plot) {
		dev.new()
		plot(0:maxnbfactors, sdt, ylab = "Variance Inflation Criterion", xlab = "Number of factors", 
			bty = "l", lwd = 1.25, type = "b", pch = 16, cex.lab = 1.25, cex = 1.25, 
			cex.axis = 1.25)
	}
	if (which.min(sdt) == 1) 
		opt <- 0
	if (which.min(sdt) > 1) {
		jumps <- -diff(sdt)/sdt[-length(sdt)]
		opt <- max((1:maxnbfactors)[jumps > jumps.nbfactor])
	}
	list(criterion = sdt, optimalnbfactors = opt)
}

emfa.nbf = function(csdata, S, eig, vy, nbf, EM = TRUE, minerr = 1e-06, verbose = FALSE) {
	n <- nrow(csdata)
	m <- ncol(csdata)
	if (((n > m) & (m <= 200) & (m >= 3)) & (!EM)) {
		if (nbf == 0) {
			B = NULL
			Psi = rep(1, m)
		}
		if (nbf > 0) {
			fa = factanal(csdata, factors = nbf, rotation = "varimax")
			B = fa$loadings
			class(B) = "matrix"
			Psi = fa$uniquenesses
			Psi = Psi * vy
			B = matrix(rep(sqrt(vy), ncol(B)), nrow = nrow(B)) * B
			sB = scale(t(B), center = FALSE, scale = sqrt(Psi))
			G = solve(diag(nbf) + tcrossprod(sB))
			sB = scale(t(B), center = FALSE, scale = Psi)
		}
	}
	if ((n <= m) | (m > 200) | EM) {
		if (nbf == 0) {
			B = NULL
			Psi = rep(1, m)
		}
		if (nbf > 0) {
			if (verbose) 
				print(paste("Fitting EM Factor Analysis Model with", nbf, "factors"))
			evectors = eig$u[, 1:nbf]
			evalues = eig$d^2
			if (nbf > 1) 
				B = evectors[, 1:nbf] * matrix(sqrt(evalues[1:nbf]), ncol = nbf, 
					nrow = m, byrow = TRUE)
			if (nbf == 1) 
				B = matrix(evectors, nrow = m, ncol = 1) * sqrt(evalues[1])
			b2 = rowSums(B^2)
			Psi = 1 - b2
			crit = 1
			while (crit > minerr) {
				inv = ifa(Psi, B)
				Cyz = crossprod(S, inv$iSB)
				Czz = crossprod(inv$iSB, Cyz) + diag(nbf) - crossprod(B, inv$iSB)
				Bnew = tcrossprod(Cyz, solve(Czz))
				Psinew = 1 - rowSums(Bnew * Cyz)
				crit = mean((Psi - Psinew)^2)
				B = Bnew
				Psi = Psinew
				if (verbose) 
					print(paste("Objective criterion in EM-FA : ", signif(crit, 6)))
			}
			Psi = Psi * vy
			B = matrix(rep(sqrt(vy), ncol(B)), nrow = nrow(B)) * B
			sB = scale(t(B), center = FALSE, scale = sqrt(Psi))
			G = solve(diag(nbf) + tcrossprod(sB))
			sB = scale(t(B), center = FALSE, scale = Psi)
		}
	}
	res = list(B = B, Psi = Psi)
	return(res)
}

LassoML <- function(data.train, ...) {
	p <- ncol(data.train$x)
	n <- nrow(data.train$x)
	nbclass <- length(unique(data.train$y))
	cl <- sort(unique(data.train$y))
	if (!all(cl == c(1:nbclass))) {
		stop("Group variable must be 1,2, ...")
	}
	family <- ifelse(nbclass == 2, "binomial", "multinomial")
	cvmod <- cv.glmnet(x = as.matrix(data.train$x), y = data.train$y, family = family, 
		type.measure = "class")
	lambda.min <- cvmod$lambda.min
	mod <- glmnet(x = as.matrix(data.train$x), y = data.train$y, family = family, 
		lambda = lambda.min, ...)
	proba.train <- predict(mod, newx = as.matrix(data.train$x), type = "response")
	if (nbclass == 2) {
		proba.train <- matrix(c(1 - proba.train, proba.train), ncol = 2, byrow = FALSE)
	}
	if (nbclass > 2) {
		proba.train <- proba.train[, , 1]
	}
	return(list(proba.train = proba.train, model = mod))
}

FADA.tmp <- function(faobject, method, sda.method, alpha,...) {
	fadta <- faobject$fa.training
	fatest <- faobject$fa.testing
	groups <- faobject$groups
	p <- ncol(faobject$fa.training)
	nbclass <- length(unique(groups))
	if (method == "glmnet") {
		out <- LassoML(list(x = fadta, y = groups), ...)
		selected <- out$selected
		proba.test <- predict(out$mod, newx = as.matrix(fatest), type = "response")
		if (nbclass == 2) {
			proba.test <- matrix(c(1 - proba.test, proba.test), ncol = 2, byrow = FALSE)
		}
		predict.test <- apply(proba.test, 1, which.max)
		out <- out$model
		proba.train <- predict(out, fadta, type = "response")
	}
	if (method == "sda") {
		ranking.LDA <- sda::sda.ranking(fadta, groups, verbose = FALSE,...)
		if (sda.method == "lfdr") {
			selected <- as.numeric(ranking.LDA[ranking.LDA[, "lfdr"] < 0.8, "idx"])
		} else {
			thr <- which.max(ranking.LDA[1:round(alpha * p), "HC"])
			selected <- as.numeric(ranking.LDA[1:thr, "idx"])
		}
		out <- sda::sda(fadta[, selected, drop = FALSE], groups, verbose = FALSE,...)
		pred <- sda::predict.sda(out, fatest[, selected, drop = FALSE], verbose = FALSE)
		proba.test <- pred$posterior
		predict.test <- pred$class
		proba.train <- sda::predict.sda(out, fadta[, selected, drop = FALSE], 
			verbose = FALSE)$posterior
	}
	if (method == "sparseLDA") {
		Xc <- normalize(fadta)
		Xn <- Xc$Xc
		out <- sparseLDA::sda(Xn, factor(groups), ...)
		Xctest <- normalizetest(fatest, Xc)
		Xctest <- matrix(Xctest, nrow = nrow(fatest), byrow = FALSE)
		colnames(Xctest) <- colnames(Xn)
		pred <- sparseLDA::predict.sda(out, Xctest)
		selected <- out$varIndex
		proba.test <- pred$posterior
		predict.test <- pred$class
		proba.train <- sparseLDA::predict.sda(out, Xn)$posterior
	}
	return(list(method = method, selected = selected, proba.train = proba.train, 
		proba.test = proba.test, predict.test = predict.test, mod = out))
}

cv.FADA <- function(train.x, train.y, test.x, test.y, nbf.cv, method,sda.method,maxnbfactors, 
	min.err, EM, maxiter, alpha,...) {
	fa.train <- decorrelate.train(list(x = train.x, y = train.y), 
		nbf = nbf.cv, maxnbfactors = maxnbfactors, diagnostic.plot = FALSE, 
		min.err = min.err, verbose = FALSE, EM = EM, maxiter = maxiter,...)
	fa.test <- decorrelate.test(fa.train, list(x = test.x))
	fada <- FADA.tmp(fa.test, method, sda.method, alpha,...)
	return(mean(fada$predict.test != test.y))
}