# FREGAT (c) 2016

genotypes <- function() {#(r, reg, genodata, mode, weights, fweights, beta.par, positions, flip.genotypes, test, et cetera) {
	if (!gtype) {
		Z <- as.matrix(genodata[, reg])
	} else if (gtype == 1) {
		Z <- as.matrix(GenABEL::as.double.snp.data(genodata[, reg]))
	} else if (gtype == 2) {
		Z <- read.plink.region(bed, snvnames, idnames, reg)
		Z <- subset.SnpMatrix(Z, measured.ids)
		Z <- SnpMatrix2numeric(Z)
	} else if (gtype == 3) {
		if (getRversion() >= "3.2.0") {
			invisible(capture.output(invisible(capture.output(Z <- seqminer::readVCFToMatrixByGene(genodata, geneFile, reg, annoType), type = 'output')), type = 'message'))
		} else { Z <- seqminer::readVCFToMatrixByGene(genodata, geneFile, reg, annoType) }
		if (is.null(Z[[1]])) return(list(m0 = 0, Z = NULL, w = NULL, pos = NULL))
		Z <- t(Z[[1]])
		if (!dim(Z)[2] == 0) Z <- Z[measured.ids, ]
	}
	if (dim(Z)[2] == 0) return(list(m0 = 0, Z = NULL, w = NULL, pos = NULL))
	w <- pos <- NULL
	if (test == 'famFLM') { pos <- positions[reg]
	} else if (!is.null(weights)) w <- weights[reg]
	Z <- as.matrix(Z)
	m0 <- dim(Z)[2]

	Z[Z == -9] <- NA

	MAF <- colMeans(Z, na.rm = TRUE) / 2

	if (any(is.na(MAF))) {
		if (all(is.na(MAF))) return(list(m0 = m0, Z = NULL, w = w, pos = pos))
		index <- !is.na(MAF)
		Z <- Z[, index]
		MAF <- MAF[index]
		if (test == 'famFLM') { pos <- as.vector(pos[index])
		} else if (!is.null(weights)) w <- as.vector(w[index])
	}

	if (any(MAF > .5) & test %in% c('famBT', 'famSKAT')) {
		if (flip.genotypes) {
			v <- MAF > .5
			Z[, v] <- abs(2 - Z[, v])
			MAF <- 1 - MAF
		} else warning("Check whether minor allele homozygotes are coded as 2 at reg ", r, ", use 'flip.genotypes = TRUE' if needed")
	}

	if (mode != 'add') {
		Z <- round(Z)
		index <- sapply(1:dim(Z)[2], function(x) all(Z[, x] >= 0, na.rm = T) & all(Z[, x] <= 2, na.rm = T))
		if (sum(index) == 0) return(list(m0 = m0, Z = 'not in range', w = NULL, pos = NULL))
		Z <- Z[, index]
		MAF <- MAF[index]
		if (test == 'famFLM') { pos <- as.vector(pos[index])
		} else if (!is.null(weights)) w <- as.vector(w[index])
		if (mode == 'rec') Z[Z == 1] <- 0
		Z[Z == 2] <- 1
	}

	v <- sapply(1:dim(Z)[2], function(x) var(Z[, x], na.rm = TRUE))
	Z <- as.matrix(Z[, v > 0])
	if (test == 'famFLM') { pos <- pos[v > 0]
	} else if (!is.null(weights)) { w <- w[v > 0]
	} else if (!is.null(fweights)) MAF <- MAF[v > 0]
	if (dim(Z)[2] == 0) return(list(m0 = m0, Z = NULL, w = NULL, pos = NULL))

	#impute missing genotypes
	if (mode == 'add') {
		if (impute.method == 'mean') {
			for (z in 1:dim(Z)[2]) Z[is.na(Z[, z]), z] <- mean(Z[, z], na.rm = T)
		} else { # BLUE imputation
			for (z in 1:dim(Z)[2]) {
				gen <- as.vector(Z[, z]) # these are genotypes of all individuals for one marker
				BLUE <- sum(colInvOmega[!is.na(gen)] * gen[!is.na(gen)]) / sum(colInvOmega[!is.na(gen)])
				Z[is.na(gen), z] <- BLUE
			}
		}
	} else {
		if (impute.method == 'mean' & mode == 'dom') {
			for (z in 1:dim(Z)[2]) {
				q <- mean(Z[, z], na.rm = T) / 2
				Z[is.na(Z[, z]), z] <- 2 * q - q ^ 2
			}
		} else if (impute.method == 'mean' & mode == 'rec') {
			for (z in 1:dim(Z)[2]) {
				q <- mean(Z[, z], na.rm = T) / 2
				Z[is.na(Z[, z]), z] <- q ^ 2
			}
		} else if (impute.method == 'blue' & mode == 'dom') {
			for (z in 1:dim(Z)[2]) {
				gen <- as.vector(Z[, z])
				q <- sum(colInvOmega[!is.na(gen)] * gen[!is.na(gen)]) / sum(colInvOmega[!is.na(gen)]) / 2
				Z[is.na(gen), z] <- 2 * q - q ^ 2
			}
		} else if (impute.method == 'blue' & mode == 'rec') {
			for (z in 1:dim(Z)[2]) {
				gen <- as.vector(Z[, z])
				q <- sum(colInvOmega[!is.na(gen)] * gen[!is.na(gen)]) / sum(colInvOmega[!is.na(gen)]) / 2
				Z[is.na(gen), z] <- q ^ 2
			}
		}
	}

	if (test == 'famBT') {
		if (!is.null(fweights)) w <- fweights(MAF)
		return(list(m0 = m0, Z = Z, w = w, pos = pos))
	}

	if (omit.linear.dependent) {
		# detect and eliminate linear-dependent genetic variants
		dqr <- qr(Z) # QR-decomposition
		index <- dqr$pivot[1:dqr$rank] # indexes of dependent genetic variants
#	} else {
#		index <- !duplicated(Z, MARGIN = 2)
#	}
		Z <- as.matrix(Z[, index]) # elimination
		if (test == 'famFLM') {
			pos <- as.vector(pos[index])
		} else if (!is.null(weights)) {
			w <- as.vector(w[index])
		} else if (!is.null(fweights)) {
			MAF <- MAF[index]
			w <- fweights(MAF)
		}
	} else {
		if (!is.null(fweights)) w <- fweights(MAF)
	}
	if (test == 'famFLM') {
#		pos <- as.vector(pos[index])
		if (any(is.na(pos))) {
			pos <- NULL
		} else {
			Z <- as.matrix(Z[, order(pos)])
			pos <- pos[order(pos)]
			if (sum(duplicated(pos)) > 0) {
				for (i in 1:length(pos)) {
					pos[duplicated(pos)] <- pos[duplicated(pos)] + 1
					if (sum(duplicated(pos)) == 0) break()
				}
			}
		if (max(pos) > min(pos)) pos <- (pos - min(pos)) / (max(pos) - min(pos))
		if (max(pos) == min(pos)) pos <- .5
		}
	if (flip.genotypes & dim(Z)[2] > 1) Z <- flipper(Z)
	}
	#return(list(beta.par = beta.par, fweights = fweights, w = w,maf=MAF))
	return(list(m0 = m0, Z = Z, w = w, pos = pos))
}

flipper <- function(Z) {
	for (j in 1:(dim(Z)[2] - 1)) {
		A <- sum(abs(Z[, j] - Z[, j + 1]))
		B <- sum(abs(Z[, j] - (2 - Z[, j + 1])))
		if (A > B) Z[, j + 1] <- 2 - Z[, j + 1] 
	}
	return(Z)
}
