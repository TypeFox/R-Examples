notests <- FALSE
if(notests) q(save = "no")
stopifnot(require(FAiR))

lower <- sqrt(.Machine$double.eps)
man <- make_manifest(covmat = Harman74.cor)
set.seed(12345)
par <- runif(99, max = 1/3)
par[76:99] <- log(par[52:75])

Phi  <- diag(3)
free <- lower.tri(Phi)
Phi  <- new("parameter.cormat", x = Phi, free = free, num_free = sum(free))

fixed <- matrix(NA_real_, nrow = nrow(cormat(man)), ncol = 3)
rownames(fixed) <- rownames(cormat(man))
free <- is.na(fixed)

# basic mapping rule
beta <- new("parameter.coef.SEFA", x = fixed, free = free, num_free = sum(free),
		mapping_rule = mapping_rule)
res <- make_restrictions(manifest = man, beta = beta, Phi = Phi, 
			discrepancy = "MLE", criteria = list())
fits <- fitS4(par, res, man, lower, TRUE)
stopifnot(which(fits != -1)[1] == length(fits))
res2 <- restrictions2model(par, res, man, lower, TRUE)$restrictions
!coef(res2)

# Thurstone simple structure
formals(res@beta@mapping_rule)$row_complexity <- 2
fits <- fitS4(par, res, man, lower, TRUE)
stopifnot(which(fits != -1)[1] == length(fits))
res2 <- restrictions2model(par, res, man, lower, TRUE)$restrictions
!coef(res2)

# unspecific exclusions by row
formals(res@beta@mapping_rule)$row_complexity <- rep(1:2, 6)
fits <- fitS4(par, res, man, lower, TRUE)
stopifnot(which(fits != -1)[1] == length(fits))
res2 <- restrictions2model(par, res, man, lower, TRUE)$restrictions
!coef(res2)

# quasi Yates
formals(res@beta@mapping_rule)$row_complexity <- NA_integer_
formals(res@beta@mapping_rule)$quasi_Yates <- TRUE
fits <- fitS4(par, res, man, lower, TRUE)
stopifnot(which(fits != -1)[1] == length(fits))
res2 <- restrictions2model(par, res, man, lower, TRUE)$restrictions
!coef(res2)

# weak Thurstone
formals(res@beta@mapping_rule)$quasi_Yates <- FALSE
formals(res@beta@mapping_rule)$weak_Thurstone <- TRUE
fits <- fitS4(par, res, man, lower, TRUE)
stopifnot(which(fits != -1)[1] == length(fits))
res2 <- restrictions2model(par, res, man, lower, TRUE)$restrictions
!coef(res2)

# Butler
formals(res@beta@mapping_rule)$weak_Thurstone <- FALSE
formals(res@beta@mapping_rule)$Butler <- TRUE
fits <- fitS4(par, res, man, lower, TRUE)
stopifnot(which(fits != -1)[1] == length(fits))
res2 <- restrictions2model(par, res, man, lower, TRUE)$restrictions
!coef(res2)

# viral
formals(res@beta@mapping_rule)$Butler <- FALSE
formals(res@beta@mapping_rule)$viral <- TRUE
fits <- fitS4(par, res, man, lower, TRUE)
stopifnot(which(fits != -1)[1] == length(fits))
res2 <- restrictions2model(par, res, man, lower, TRUE)$restrictions
!coef(res2)

# high communality
formals(res@beta@mapping_rule)$viral <- FALSE
formals(res@beta@mapping_rule)$communality <- TRUE
fits <- fitS4(par, res, man, lower, TRUE)
stopifnot(which(fits != -1)[1] == length(fits))
res2 <- restrictions2model(par, res, man, lower, TRUE)$restrictions
!coef(res2)

## Inequality restrictions
formals(res@beta@mapping_rule)$communality <- FALSE

# ranks_rows_1st
ranks <- matrix(3, nrow = 24, ncol = 3, dimnames = dimnames(fixed))
ranks[1,1] <- 2
res <- make_restrictions(man, beta = beta, Phi = Phi, 
			criteria = list("ranks_rows_1st"),
			methodArgs = list(row_ranks = ranks))
fits <- fitS4(par, res, man, lower, TRUE)
stopifnot(which(fits != -1)[1] == length(fits))
res2 <- restrictions2model(par, res, man, lower, TRUE)$restrictions
(coef(res2) * (coef(res2) %*% cormat(res2)))[1,]

# ranks_cols_1st
ranks <- matrix(24, nrow = 24, ncol = 3, dimnames = dimnames(fixed))
ranks[1,1] <- 23
res <- make_restrictions(man, beta = beta, Phi = Phi, 
			criteria = list("ranks_cols_1st"),
			methodArgs = list(col_ranks = ranks))
fits <- fitS4(par, res, man, lower, TRUE)
stopifnot(which(fits != -1)[1] == length(fits))
res2 <- restrictions2model(par, res, man, lower, TRUE)$restrictions
(FC <- coef(res2) * (coef(res2) %*% cormat(res2)))[,1]

# indicators_1st
indicators <- apply(FC, 2, which.max)
res <- make_restrictions(man, beta = beta, Phi = Phi, 
			criteria = list("indicators_1st"),
			methodArgs = list(indicators = indicators))
fits <- fitS4(par, res, man, lower, TRUE)
stopifnot(which(fits != -1)[1] == length(fits))

# evRF_1st
res <- make_restrictions(man, beta = beta, Phi = Phi, 
			criteria = list("evRF_1st"),
			methodArgs = list())
fits <- fitS4(par, res, man, lower, TRUE)
stopifnot(which(fits != -1)[1] == length(fits))
res2 <- restrictions2model(par, res, man, lower, TRUE)$restrictions
PF  <- cormat(res2)
RF  <- cov2cor(chol2inv(chol(PF)))
C <- fitted(res2, reduced = FALSE, standardized = TRUE)
det(C)^(1/ncol(C))
det(RF)^(1/ncol(RF))

# evPF_1st
res <- make_restrictions(man, beta = beta, Phi = Phi, 
			criteria = list("evPF_1st"),
			methodArgs = list())
fits <- fitS4(par, res, man, lower, TRUE)
stopifnot(which(fits != -1)[1] == length(fits))
res2 <- restrictions2model(par, res, man, lower, TRUE)$restrictions
PF <- cormat(res2)
C <- fitted(res2, reduced = FALSE, standardized = TRUE)
det(C)^(1/ncol(C))
det(PF)^(1/ncol(PF))

# h2_over_FC_1st
res <- make_restrictions(man, beta = beta, Phi = Phi, 
			criteria = list("h2_over_FC_1st"),
			methodArgs = list())
fits <- fitS4(par, res, man, lower, TRUE)
stopifnot(which(fits != -1)[1] == length(fits))
res2 <- restrictions2model(par, res, man, lower, TRUE)$restrictions
(FC <- coef(res2) * (coef(res2) %*% cormat(res2)))
cbind(FC, NA, rowSums(FC))

# no_neg_suppressors_1st
res <- make_restrictions(man, beta = beta, Phi = Phi, 
			criteria = list("no_neg_suppressors_1st"),
			methodArgs = list(FC_threshold = 0))
fits <- fitS4(par, res, man, lower, TRUE)
stopifnot(which(fits != -1)[1] == length(fits))
res2 <- restrictions2model(par, res, man, lower, TRUE)$restrictions
(FC <- coef(res2) * (coef(res2) %*% cormat(res2)))
FC

# gv_1st
res <- make_restrictions(man, beta = beta, Phi = Phi,
			criteria = list("gv_1st"),
			methodArgs = list())
fits <- fitS4(par, res, man, lower, TRUE)
stopifnot(which(fits != -1)[1] == length(fits))
res2 <- restrictions2model(par, res, man, lower, TRUE)$restrictions
PF  <- cormat(res2)
RF  <- cov2cor(chol2inv(chol(PF)))
det(PF)
det(RF)

# distinguishability_1st
res <- make_restrictions(man, beta = beta, Phi = Phi,
			criteria = list("distinguishability_1st"),
			methodArgs = list())
fits <- fitS4(par, res, man, lower, TRUE)
stopifnot(which(fits != -1)[1] == length(fits))
res2 <- restrictions2model(par, res, man, lower, TRUE)$restrictions
(FC <- coef(res2) * (coef(res2) %*% cormat(res2)))
FC

# cohyperplanarity_1st
res <- make_restrictions(man, beta = beta, Phi = Phi,
			criteria = list("cohyperplanarity_1st"),
			methodArgs = list())
fits <- fitS4(par, res, man, lower, TRUE)
stopifnot(which(fits != -1)[1] == length(fits))
res2 <- restrictions2model(par, res, man, lower, TRUE)$restrictions
PP <- coef(res2)
C <- fitted(res2, reduced = TRUE, standardized = TRUE)
diag(C) <- 1
det(C)^(1/nrow(C))
mark <- which(!PP[,1])
det(C[mark,mark])^(1/length(mark))
mark <- which(!PP[,2])
det(C[mark,mark])^(1/length(mark))
mark <- which(!PP[,3])
det(C[mark,mark])^(1/length(mark))

# dist_cols_1st
res <- make_restrictions(man, beta = beta, Phi = Phi,
			criteria = list("dist_cols_1st"),
			methodArgs = list(cutpoint = 0.5))
fits <- fitS4(par, res, man, lower, TRUE)
stopifnot(which(fits != -1)[1] == length(fits))
res2 <- restrictions2model(par, res, man, lower, TRUE)$restrictions
PP <- coef(res2)
dist(t(!PP), method = "binary")

# volume_1st
res <- make_restrictions(man, beta = beta, Phi = Phi,
			criteria = list("volume_1st"),
			methodArgs = list())
fits <- fitS4(par, res, man, lower, TRUE)
stopifnot(which(fits != -1)[1] == length(fits))
res2 <- restrictions2model(par, res, man, lower, TRUE)$restrictions
C <- fitted(res2, reduced = FALSE, standardized = TRUE)
S <- cormat(man)
det(S)
det(C)

# block_1st
blockers <- matrix(NA, nrow = 24, ncol = 3, dimnames = dimnames(fixed))
blockers[1,1] <- TRUE
res <- make_restrictions(man, beta = beta, Phi = Phi,
			criteria = list("block_1st"),
			methodArgs = list(blockers = blockers))
fits <- fitS4(par, res, man, lower, TRUE)
stopifnot(which(fits != -1)[1] == length(fits))
res2 <- restrictions2model(par, res, man, lower, TRUE)$restrictions
coef(res2)[1,]

## Exact restrictions
# one fixed to zero
fixed[1,3] <- 0
free <- is.na(fixed)
beta <- new("parameter.coef.SEFA", x = fixed, free = free, num_free = sum(free),
		mapping_rule = mapping_rule)
res <- make_restrictions(manifest = man, beta = beta, Phi = Phi, 
			discrepancy = "MLE", criteria = list())

fits <- fitS4(par[-75], res, man, lower, TRUE)
stopifnot(which(fits != -1)[1] == length(fits))
res2 <- restrictions2model(par[-75], res, man, lower, TRUE)$restrictions
coef(res2)[1,]

# one "nonlinear" restriction (really just an equality restriction)
fixed[1,3] <- NA_real_
fixed[2,1] <- Inf
free <- is.na(fixed)
foo <- function(x) {
	x[2,1] <- x[1,1]
	return(x)
}
beta <- new("parameter.coef.SEFA.nl", x = fixed, free = free, num_free = sum(free),
		mapping_rule = mapping_rule, nonlinearities = foo)
res <- make_restrictions(manifest = man, beta = beta, Phi = Phi, 
			discrepancy = "MLE", criteria = list())
fits <- fitS4(par[-2], res, man, lower, TRUE)
stopifnot(which(fits != -1)[1] == length(fits))
res2 <- restrictions2model(par[-2], res, man, lower, TRUE)$restrictions
coef(res2)[1:2,]

# one equality the conventional way
er <- new("equality_restriction", free = 1L, fixed = 2L, dims = dim(fixed),
		rownames = rownames(cormat(man)), level = 1L)
beta <- new("parameter.coef.SEFA", x = fixed, free = free, num_free = sum(free),
		mapping_rule = mapping_rule, equalities = list(er))
fits <- fitS4(par[-2], res, man, lower, TRUE)
stopifnot(which(fits != -1)[1] == length(fits))
res2 <- restrictions2model(par[-2], res, man, lower, TRUE)$restrictions
coef(res2)[1:2,]

## redo with a matrix of data
man <- make_manifest(x = USJudgeRatings[,-c(1:5)], how = "unbiased")
man <- make_manifest(data = USJudgeRatings[,-c(1:5)], how = "unbiased")
man <- make_manifest(x = as.matrix(USJudgeRatings[,-c(1:5)]), how = "unbiased")
man <- make_manifest(data = as.matrix(USJudgeRatings[,-c(1:5)]), how = "unbiased")
set.seed(2)
par <- runif(31, max = 0.5)
par[25:31] <- log(par[25:31])

Phi  <- diag(3)
free <- lower.tri(Phi)
Phi  <- new("parameter.cormat", x = Phi, free = free, num_free = sum(free))

fixed <- matrix(NA_real_, nrow = nrow(cormat(man)), ncol = 3)
rownames(fixed) <- rownames(cormat(man))
free <- is.na(fixed)

beta <- new("parameter.coef.SEFA", x = fixed, free = free, num_free = sum(free),
		mapping_rule = mapping_rule)

# ADF
res <- make_restrictions(manifest = man, beta = beta, Phi = Phi, 
			discrepancy = "ADF", criteria = list())
fits <- fitS4(par, res, man, lower, TRUE)
stopifnot(which(fits != -1)[1] == length(fits))

# SHK
res <- make_restrictions(manifest = man, beta = beta, Phi = Phi, 
			discrepancy = "SHK", criteria = list())
fits <- fitS4(par, res, man, lower, TRUE)
stopifnot(which(fits != -1)[1] == length(fits))

# HK (not posdef)
# res <- make_restrictions(manifest = man, beta = beta, Phi = Phi, 
# 			discrepancy = "HK", criteria = list())
# fits <- fitS4(par, res, man, lower, TRUE)
# stopifnot(which(fits != -1)[1] == length(fits))

# ELLIPTICAL
res <- make_restrictions(manifest = man, beta = beta, Phi = Phi, 
			discrepancy = "ELLIPTICAL", criteria = list())
fits <- fitS4(par, res, man, lower, TRUE)
stopifnot(which(fits != -1)[1] == length(fits))

# MLE
res <- make_restrictions(manifest = man, beta = beta, Phi = Phi, 
			discrepancy = "MLE", criteria = list())
fits <- fitS4(par, res, man, lower, TRUE)
stopifnot(which(fits != -1)[1] == length(fits))

# ranks_rows_1st
ranks <- matrix(3, nrow = 7, ncol = 3, dimnames = dimnames(fixed))
ranks[1,1] <- 2
res <- make_restrictions(manifest = man, beta = beta, Phi = Phi,
			criteria = list("ranks_rows_1st"),
			methodArgs = list(row_ranks = ranks))
fits <- fitS4(par, res, man, lower, TRUE)
stopifnot(which(fits != -1)[1] == length(fits))
res2 <- restrictions2model(par, res, man, lower, TRUE)$restrictions
(FC <- coef(res2) * (coef(res2) %*% cormat(res2)))[1,]

# ranks_cols_1st
ranks <- matrix(7, nrow = 7, ncol = 3, dimnames = dimnames(fixed))
ranks[1,1] <- 6
res <- make_restrictions(manifest = man, beta = beta, Phi = Phi,
			criteria = list("ranks_cols_1st"),
			methodArgs = list(col_ranks = ranks))
fits <- fitS4(par, res, man, lower, TRUE)
stopifnot(which(fits != -1)[1] == length(fits))
res2 <- restrictions2model(par, res, man, lower, TRUE)$restrictions
(FC <- coef(res2) * (coef(res2) %*% cormat(res2)))[,1]

# indicators_1st
indicators <- apply(FC, 2, which.max)
res <- make_restrictions(manifest = man, beta = beta, Phi = Phi,
			criteria = list("indicators_1st"),
			methodArgs = list(indicators = indicators))
fits <- fitS4(par, res, man, lower, TRUE)
stopifnot(which(fits != -1)[1] == length(fits))

# evRF_1st
res <- make_restrictions(manifest = man, beta = beta, Phi = Phi,
			criteria = list("evRF_1st"),
			methodArgs = list())
fits <- fitS4(par, res, man, lower, TRUE)
stopifnot(which(fits != -1)[1] == length(fits))
res2 <- restrictions2model(par, res, man, lower, TRUE)$restrictions
PF <- cormat(res2)
RF <- cov2cor(chol2inv(chol(PF)))
C <- fitted(res2, reduced = FALSE, standardized = TRUE)
det(C)^(1/ncol(C))
det(RF)^(1/ncol(RF))

# evPF_1st
res <- make_restrictions(manifest = man, beta = beta, Phi = Phi,
			criteria = list("evPF_1st"),
			methodArgs = list())
fits <- fitS4(par, res, man, lower, TRUE)
stopifnot(which(fits != -1)[1] == length(fits))
res2 <- restrictions2model(par, res, man, lower, TRUE)$restrictions
PF <- cormat(res2)
C  <- fitted(res2, reduced = FALSE, standardized = TRUE)
det(C)^(1/ncol(C))
det(PF)^(1/ncol(PF))

# h2_over_FC_1st
res <- make_restrictions(manifest = man, beta = beta, Phi = Phi,
			criteria = list("h2_over_FC_1st"),
			methodArgs = list())
fits <- fitS4(par, res, man, lower, TRUE)
stopifnot(which(fits != -1)[1] == length(fits))
res2 <- restrictions2model(par, res, man, lower, TRUE)$restrictions
(FC <- coef(res2) * (coef(res2) %*% cormat(res2)))

# no_neg_suppressors_1st
res <- make_restrictions(manifest = man, beta = beta, Phi = Phi,
			criteria = list("no_neg_suppressors_1st"),
			methodArgs = list(FC_threshold = 0))
fits <- fitS4(par, res, man, lower, TRUE)
stopifnot(which(fits != -1)[1] == length(fits))
res2 <- restrictions2model(par, res, man, lower, TRUE)$restrictions
(FC <- coef(res2) * (coef(res2) %*% cormat(res2)))

# gv_1st
res <- make_restrictions(manifest = man, beta = beta, Phi = Phi,
			criteria = list("gv_1st"),
			methodArgs = list())
fits <- fitS4(par, res, man, lower, TRUE)
stopifnot(which(fits != -1)[1] == length(fits))
res2 <- restrictions2model(par, res, man, lower, TRUE)$restrictions
PF <- cormat(res2)
RF <- cov2cor(chol2inv(chol(PF)))
det(PF)
det(RF)

# distinguishability_1st
res <- make_restrictions(manifest = man, beta = beta, Phi = Phi,
			criteria = list("distinguishability_1st"),
			methodArgs = list())
fits <- fitS4(par, res, man, lower, TRUE)
stopifnot(which(fits != -1)[1] == length(fits))
res2 <- restrictions2model(par, res, man, lower, TRUE)$restrictions
(FC <- coef(res2) * (coef(res2) %*% cormat(res2)))

# cohyperplanarity_1st
res <- make_restrictions(manifest = man, beta = beta, Phi = Phi,
			criteria = list("cohyperplanarity_1st"),
			methodArgs = list())
fits <- fitS4(par, res, man, lower, TRUE)
stopifnot(which(fits != -1)[1] == length(fits))
res2 <- restrictions2model(par, res, man, lower, TRUE)$restrictions
PP <- loadings(res2)
C <- fitted(res2, reduced = FALSE, standardized = TRUE)
det(C)^(1/ncol(C))
mark <- which(!PP[,1])
det(C[mark,mark])^(1/length(mark))
mark <- which(!PP[,2])
det(C[mark,mark])^(1/length(mark))
mark <- which(!PP[,3])
det(C[mark,mark])^(1/length(mark))

# dist_cols_1st
res <- make_restrictions(manifest = man, beta = beta, Phi = Phi,
			criteria = list("dist_cols_1st"),
			methodArgs = list(cutpoint = 0.5))
fits <- fitS4(par, res, man, lower, TRUE)
stopifnot(which(fits != -1)[1] == length(fits))
res2 <- restrictions2model(par, res, man, lower, TRUE)$restrictions
PP <- coef(res2)
dist(t(!PP), method = "binary")

# volume_1st
res <- make_restrictions(manifest = man, beta = beta, Phi = Phi,
			criteria = list("volume_1st"),
			methodArgs = list())
fits <- fitS4(par, res, man, lower, TRUE)
stopifnot(which(fits != -1)[1] == length(fits))
res2 <- restrictions2model(par, res, man, lower, TRUE)$restrictions
C <- fitted(res2, reduced = FALSE, standardized = TRUE)
S <- cormat(man) 
det(S)
det(C)

# block_1st
blockers <- matrix(NA, nrow = 7, ncol = 3, dimnames = dimnames(fixed))
blockers[1,1] <- TRUE
res <- make_restrictions(manifest = man, beta = beta, Phi = Phi,
			criteria = list("block_1st"),
			methodArgs = list(blockers = blockers))
fits <- fitS4(par, res, man, lower, TRUE)
stopifnot(which(fits != -1)[1] == 5) # sensibly fails
res2 <- restrictions2model(par, res, man, lower, TRUE)$restrictions
coef(res2)[1,]

## general second-order factor
fixed2 <- matrix(NA_real_, nrow = 3, ncol = 1)
rownames(fixed2) <- paste("F", 1:3, sep = "")
free  <- is.na(fixed2)
Delta <- new("parameter.coef", x = fixed2, free = free, num_free = sum(free))

ranks <- matrix(c(2,3,3), ncol = 1)
res <- make_restrictions(manifest = man, beta = beta, Delta = Delta,
			criteria = list("ranks_cols_2nd"),
			methodArgs = list(col_ranks = ranks))
set.seed(2)
par <- runif(31, max = 0.5)
par[26:31] <- log(par[26:31])
fits <- fitS4(par, res, man, lower, TRUE)
stopifnot(which(fits != -1)[1] == 6) # sensibly fails
res2 <- restrictions2model(par, res, man, lower, TRUE)$restrictions
rowSums(res2@Delta@x^2)

## exact restrictions at level 2
# one fixed to zero
fixed2[1,1] <- 0
free <- is.na(fixed2)
Delta <- new("parameter.coef", x = fixed2, free = free, num_free = sum(free))

res <- make_restrictions(manifest = man, beta = beta, Delta = Delta)
fits <- fitS4(par[-1], res, man, lower, TRUE)
stopifnot(which(fits != -1)[1] == length(fits))
res2 <- restrictions2model(par[-1], res, man, lower, TRUE)$restrictions
res2@Delta@x

# one "nonlinear" restriction (really just an equality restriction)
foo <- function(x) {
	x[2,1] <- x[1,1]
	return(x)
}
fixed2[1,1] <- NA_real_
fixed2[2,1] <- Inf
free <- is.na(fixed2)
Delta <- new("parameter.coef.nl", x = fixed2, free = free, num_free = sum(free),
		nonlinearities = foo)
res <- make_restrictions(manifest = man, beta = beta, Delta = Delta)
fits <- fitS4(par[-1], res, man, lower, TRUE)
stopifnot(which(fits != -1)[1] == length(fits))
res2 <- restrictions2model(par[-2], res, man, lower, TRUE)$restrictions
res2@Delta@x

# one equality the conventional way
er <- new("equality_restriction", free = 1L, fixed = 2L, dims = dim(fixed2),
		rownames = rownames(fixed2), level = 2L)
Delta <- new("parameter.coef", x = fixed2, free = free, num_free = sum(free),
		equalities = list(er))
res <- make_restrictions(manifest = man, beta = beta, Delta = Delta)
fits <- fitS4(par[-1], res, man, lower, TRUE)
stopifnot(which(fits != -1)[1] == length(fits))
res2 <- restrictions2model(par[-2], res, man, lower, TRUE)$restrictions
res2@Delta@x

## two second-order factors
man <- make_manifest(covmat = Harman74.cor)
fixed <- matrix(NA_real_, nrow = nrow(cormat(man)), ncol = 5)
rownames(fixed) <- rownames(cormat(man))
free <- is.na(fixed)
beta <- new("parameter.coef.SEFA", x = fixed, free = free, num_free = sum(free))

fixed2 <- matrix(NA_real_, 5, 2)
rownames(fixed2) <- LETTERS[1:5]
set.seed(7)
par <- runif(155, max = 0.5)
par[132:155] <- log(par[132:155])

free <- is.na(fixed2)
Delta <- new("parameter.coef.SEFA", x = fixed2, free = free, num_free = sum(free),
		rankcheck = "howe")

Xi <- diag(2)
free <- lower.tri(Xi)
Xi <- new("parameter.cormat", x = Xi, free = free, num_free = sum(free))

# ranks_rows_2nd
ranks <- matrix(2, nrow = 5, ncol = 2, dimnames = dimnames(fixed2))
ranks[1,1] <- 1
res <- make_restrictions(manifest = man, beta = beta, Delta = Delta, Xi = Xi,
			criteria = list("ranks_rows_2nd"),
			methodArgs = list(row_ranks_2nd = ranks))
fits <- fitS4(par, res, man, lower, TRUE)
stopifnot(which(fits != -1)[1] == length(fits))
res2 <- restrictions2model(par, res, man, lower, TRUE)$restrictions
FC <- loadings(res2, level = 2) * (loadings(res2, level = 2) %*% cormat(res2, level = 2))
FC[1,]

# ranks_cols_2nd
ranks <- matrix(5, nrow = 5, ncol = 2, dimnames = dimnames(fixed2))
ranks[1,1] <- 4
res <- make_restrictions(manifest = man, beta = beta, Delta = Delta, Xi = Xi,
			criteria = list("ranks_cols_2nd"),
			methodArgs = list(col_ranks_2nd = ranks))
fits <- fitS4(par, res, man, lower, TRUE)
stopifnot(which(fits != -1)[1] == length(fits))
res2 <- restrictions2model(par, res, man, lower, TRUE)$restrictions
FC <- loadings(res2, level = 2) * (loadings(res2, level = 2) %*% cormat(res2, level = 2))
FC[,1]

# indicators_2nd
indicators <- apply(FC, 2, which.max)
res <- make_restrictions(manifest = man, beta = beta, Delta = Delta, Xi = Xi,
			criteria = list("indicators_2nd"),
			methodArgs = list(indicators_2nd = indicators))
fits <- fitS4(par, res, man, lower, TRUE)
stopifnot(which(fits != -1)[1] == length(fits))

# evRF_2nd
res <- make_restrictions(manifest = man, beta = beta, Delta = Delta, Xi = Xi,
			criteria = list("evRF_2nd"), methodArgs = list())
fits <- fitS4(par, res, man, lower, TRUE)
stopifnot(which(fits != -1)[1] == 8) # sensibly fails
res2 <- restrictions2model(par, res, man, lower, TRUE)$restrictions
PF  <- cormat(res2, level = 2)
RF  <- cov2cor(chol2inv(chol(PF)))
PF1 <- cormat(res2)
det(PF1)^(1/ncol(PF1))
det(RF)^(1/ncol(RF))

# evPF_2nd
res <- make_restrictions(manifest = man, beta = beta, Delta = Delta, Xi = Xi,
			criteria = list("evPF_2nd"), methodArgs = list())
fits <- fitS4(par, res, man, lower, TRUE)
stopifnot(which(fits != -1)[1] == 8) # sensibly fails
res2 <- restrictions2model(par, res, man, lower, TRUE)$restrictions
PF1  <- cormat(res2, level = 1)
PF   <- cormat(res2, level = 2)
det(PF1)^(1/ncol(PF1))
det(PF)^(1/ncol(PF))

# h2_over_FC_2nd
res <- make_restrictions(manifest = man, beta = beta, Delta = Delta, Xi = Xi,
			criteria = list("h2_over_FC_2nd"))
fits <- fitS4(par, res, man, lower, TRUE)
stopifnot(which(fits != -1)[1] == length(fits))
res2 <- restrictions2model(par, res, man, lower, TRUE)$restrictions
FC <- loadings(res2, level = 2) * (loadings(res2, level = 2) %*% cormat(res2, level = 2))
FC

# no_neg_suppressors_2nd
res <- make_restrictions(manifest = man, beta = beta, Delta = Delta, Xi = Xi,
			criteria = list("no_neg_suppressors_2nd"),
			methodArgs = list(FC_threshold_2nd = 0))
fits <- fitS4(par, res, man, lower, TRUE)
stopifnot(which(fits != -1)[1] == length(fits))
res2 <- restrictions2model(par, res, man, lower, TRUE)$restrictions
FC <- loadings(res2, level = 2) * (loadings(res2, level = 2) %*% cormat(res2, level = 2))
FC

# gv_2nd
res <- make_restrictions(manifest = man, beta = beta, Delta = Delta, Xi = Xi,
			criteria = list("gv_2nd"), methodArgs = list())
fits <- fitS4(par, res, man, lower, TRUE)
stopifnot(which(fits != -1)[1] == length(fits))
res2 <- restrictions2model(par, res, man, lower, TRUE)$restrictions
PF  <- cormat(res2, level = 2)
RF  <- cov2cor(chol2inv(chol(PF)))
det(PF)
det(RF)

# distinguishability_2nd
res <- make_restrictions(manifest = man, beta = beta, Delta = Delta, Xi = Xi,
			criteria = list("distinguishability_2nd"),
			methodArgs = list())
fits <- fitS4(par, res, man, lower, TRUE)
stopifnot(which(fits != -1)[1] == length(fits))
res2 <- restrictions2model(par, res, man, lower, TRUE)$restrictions
FC <- loadings(res2, level = 2) * (loadings(res2, level = 2) %*% cormat(res2, level = 2))
FC

# dist_cols_2nd
res <- make_restrictions(manifest = man, beta = beta, Delta = Delta, Xi = Xi,
			criteria = list("dist_cols_2nd"),
			methodArgs = list(cutpoint_2nd = 0.5))
fits <- fitS4(par, res, man, lower, TRUE)
stopifnot(which(fits != -1)[1] == length(fits))
res2 <- restrictions2model(par, res, man, lower, TRUE)$restrictions
dist(t(!res2@Delta@x), method = "binary")

# cohyperplanarity_2nd
res <- make_restrictions(manifest = man, beta = beta, Delta = Delta, Xi = Xi,
			criteria = list("cohyperplanarity_2nd"),
			methodArgs = list(cutpoint_2nd = 0.5))
fits <- fitS4(par, res, man, lower, TRUE)
stopifnot(which(fits != -1)[1] == length(fits))
res2 <- restrictions2model(par, res, man, lower, TRUE)$restrictions
PP <- loadings(res2, level = 2)
PF <- cormat(res2)
det(PF)^(1/ncol(PF))
mark <- which(!PP[,1])
det(PF[mark,mark])^(1/length(mark))
mark <- which(!PP[,2])
det(PF[mark,mark])^(1/length(mark))
