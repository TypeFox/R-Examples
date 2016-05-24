notests <- TRUE
if(notests) q(save = "no")
stopifnot(require(FAiR))

## Globals
tol <- 0.006

## Compare CFA with example from library(sem)
semC <- rbind(
c(1.0000000, 0.8267444, 0.7745019, 0.4859589, 0.4635265, 0.4085056, 0.4732969, 0.4365315, 0.4264270),
c(0.8267444, 1.0000000, 0.7823011, 0.4908525, 0.4681942, 0.4126192, 0.4780630, 0.4409274, 0.4307211),
c(0.7745019, 0.7823011, 1.0000000, 0.4598352, 0.4386087, 0.3865455, 0.4478539, 0.4130649, 0.4035036),
c(0.4859589, 0.4908525, 0.4598352, 1.0000000, 0.6662540, 0.5871692, 0.4158054, 0.3835059, 0.3746288),
c(0.4635265, 0.4681942, 0.4386087, 0.6662540, 1.0000000, 0.5600648, 0.3966113, 0.3658028, 0.3573354),
c(0.4085056, 0.4126192, 0.3865455, 0.5871692, 0.5600648, 1.0000000, 0.3495332, 0.3223817, 0.3149195),
c(0.4732969, 0.4780630, 0.4478539, 0.4158054, 0.3966113, 0.3495332, 1.0000000, 0.5623102, 0.5492943),
c(0.4365315, 0.4409274, 0.4130649, 0.3835059, 0.3658028, 0.3223817, 0.5623102, 1.0000000, 0.5066255),
c(0.4264270, 0.4307211, 0.4035036, 0.3746288, 0.3573354, 0.3149195, 0.5492943, 0.5066255, 1.0000001))

R <- diag(9)
count <- 2
R[count, 1:(count-1)] <-    .828; count <- count + 1
R[count, 1:(count-1)] <-  c(.776, .779); count <- count + 1
R[count, 1:(count-1)] <-  c(.439, .493, .460); count <- count + 1
R[count, 1:(count-1)] <-  c(.432, .464, .425, .674); count <- count + 1
R[count, 1:(count-1)] <-  c(.447, .489, .443, .590, .541); count <- count + 1
R[count, 1:(count-1)] <-  c(.447, .432, .401, .381, .402, .288); count <- count + 1
R[count, 1:(count-1)] <-  c(.541, .537, .534, .350, .367, .320, .555); count <- count + 1
R[count, 1:(count-1)] <-  c(.380, .358, .359, .424, .446, .325, .598, .452)
R <- R + t(R)
diag(R) <- 1

varnames <- c('Sentences','Vocabulary','Sent.Completion','First.Letters','4.Letter.Words',
		'Suffixes', 'Letter.Series','Pedigrees', 'Letter.Group')
rownames(semC) <- colnames(semC) <- rownames(R) <- colnames(R) <- varnames
fixed <- matrix(0, nrow = 9, ncol = 3)
rownames(fixed)
fixed[1:3,1] <- fixed[4:6,2] <- fixed[7:9,3] <- NA_real_

manifest_Thur <- make_manifest(covmat = R, n.obs = 213)
free <- is.na(fixed)
beta <- new("parameter.coef", x= fixed, free = free, num_free = sum(free))

fixed2 <- matrix(NA_real_, 3, 1)
rownames(fixed2) <- paste("F", 1:3, sep = "")
free <- is.na(fixed2)
Delta <- new("parameter.coef", x = fixed2, free = free, num_free = sum(free))
res.Thur.1 <- make_restrictions(manifest = manifest_Thur, Delta = Delta, beta = beta,
				discrepancy = "MLE")

CFA <- Factanal(manifest = manifest_Thur, restrictions = res.Thur.1)

library(sem)
RAM <- FA2RAM(CFA)
(sem_CFA <- sem(RAM, R, N = 213))
deviance(CFA)
df.residual(CFA)

FAiR:::FAiR_stress_test(CFA)

orders_row <- fixed + 3
orders_row[is.na(orders_row)] <- 1
orders_col <- fixed + 9
orders_col[is.na(orders_col)] <- 3

x <- matrix(NA_real_, nrow = 9, ncol = 3)
rownames(x) <- rownames(cormat(manifest_Thur))
free <- is.na(x)
beta <- new("parameter.coef.SEFA", x = x, free = free, num_free = sum(free))

x <- matrix(NA_real_, nrow = 3, ncol = 1)
rownames(x) <- paste("F", 1:3, sep = "")
free <- is.na(x)
Delta <- new("parameter.coef", x = x, free = free, num_free = sum(free))

res.Thur.2 <- make_restrictions(manifest = manifest_Thur, beta = beta, Delta = Delta,
				criteria = list("no_neg_suppressors_1st",
						"ranks_rows_1st", "ranks_cols_1st"),
				methodArgs = list(FC_threshold = 0, 
						row_ranks = orders_row,
						col_ranks = orders_col))
# set.seed(12345)
# starts <- matrix(runif(1000 * res.Thur.2@nvars), nrow = 1000)
# starts <- c(loadings(CFA, level = 2), loadings(CFA, level = 1))
SEFA <- Factanal(manifest = manifest_Thur, restrictions = res.Thur.2, 
		boundary.enforcement = 1, starting.values = CFA)
show(SEFA)
summary(SEFA)

RAM <- FA2RAM(SEFA)
(sem_SEFA <- sem(RAM, R, N = 213))
deviance(SEFA)
df.residual(SEFA)

FAiR:::FAiR_stress_test(SEFA)

## Holzinger example
R.Holzinger <- matrix(c(
       1,     0,      0,      0,      0,      0,      0,      0,      0,
     .75,     1,      0,      0,      0,      0,      0,      0,      0,
     .78,   .72,      1,      0,      0,      0,      0,      0,      0,
     .44,   .52,    .47,      1,      0,      0,      0,      0,      0,
     .45,   .53,    .48,    .82,      1,      0,      0,      0,      0,
     .51,   .58,    .54,    .82,    .74,      1,      0,      0,      0,
     .21,   .23,    .28,    .33,    .37,    .35,      1,      0,      0,
     .30,   .32,    .37,    .33,    .36,    .38,    .45,      1,      0,
     .31,   .30,    .37,    .31,    .36,    .38,    .52,    .67,      1),
     9, 9, byrow=TRUE)
rownames(R.Holzinger) <- colnames(R.Holzinger) <-
    c("Word.meaning", "Sentence.completion", "Odd.words", "Mixed.arithmetic",
      "Remainders", "Missing.numbers", "Gloves", "Boots", "Hatchets")

R.Holzinger <- R.Holzinger + t(R.Holzinger)
diag(R.Holzinger) <- 1

fixed <- matrix(0, nrow = 9, ncol = 3)
fixed[1:3,1] <- fixed[4:6,2] <- fixed[7:9,3] <- NA_real_
free <- is.na(fixed)
beta <- new("parameter.coef", x = fixed, free = free, num_free = sum(free))

x <- diag(3)
free <- lower.tri(x)
Phi <- new("parameter.cormat", x = x, free = free, num_free = sum(free))

manifest.Holzinger <- make_manifest(covmat = R.Holzinger, n.obs = 696)
res_Holzinger.1 <- make_restrictions(manifest = manifest.Holzinger, beta = beta,
					Phi = Phi)

CFA <- Factanal(manifest = manifest.Holzinger, restrictions = res_Holzinger.1)
summary(CFA)
RAM <- FA2RAM(CFA)
library(sem)
(sem_CFA <- sem(RAM, R.Holzinger, N = 696))
deviance(CFA)
df.residual(CFA)

FAiR:::FAiR_stress_test(CFA)

orders_row <- fixed + 3
orders_row[is.na(orders_row)] <- 1
orders_col <- fixed + 9
orders_col[is.na(orders_col)] <- 3

x <- matrix(NA_real_, 9, 3)
free <- is.na(x)
beta <- new("parameter.coef.SEFA", x = x, free = free, num_free = sum(free))

x <- diag(3)
free <- lower.tri(x)
lower <- 0 * x
upper <- 1 * free
Domains <- array(cbind(lower, upper), c(dim(x), 2))
Phi <- new("parameter.cormat", x = x, free = free, num_free = sum(free),
		Domains = Domains)
res_Holzinger.2 <- make_restrictions(manifest = manifest.Holzinger,
				beta = beta, Phi = Phi,
				criteria = list("no_neg_suppressors_1st",
						"ranks_rows_1st", "ranks_cols_1st"),
				methodArgs = list(FC_threshold = 0, 
						row_ranks = orders_row,
						col_ranks = orders_col))

Phi <- cormat(CFA)
starts <- c(Phi[lower.tri(Phi)], loadings(CFA@restrictions), 
		log(sqrt(diag(manifest.Holzinger@cov))))
SEFA <- Factanal(manifest = manifest.Holzinger, restrictions = res_Holzinger.2,
		starting.values = starts)

RAM <- FA2RAM(SEFA)
(sem_SEFA <- sem(RAM, R.Holzinger, N = 696))
deviance(SEFA)
df.residual(SEFA)

FAiR:::FAiR_stress_test(SEFA)

## Hu et. all
library(MASS)

Lambda <- matrix(0, nrow = 15, ncol = 3)
Lambda[1:5,1]  <-  c(.7, .7, .75, .8, .8)
Lambda[6:10,2] <-  c(.7, .7, .75, .8, .8)
Lambda[11:15,3] <- c(.7, .7, .75, .8, .8)

Phi <- diag(3)
Phi[1,2] <- Phi[2,1] <- 0.5
Phi[1,3] <- Phi[3,1] <- 0.4
Phi[2,3] <- Phi[3,2] <- 0.3

Psi <- diag(c(.51, .51, .4375, .36, .36, .51, .51, .4375, .36, 
	      .36, .51, .51, .4375, .36, .36))

Sigma <- Lambda %*% Phi %*% t(Lambda) + Psi
N <- 250
fixed <- Lambda
fixed[fixed > 0] <- NA_real_

X <- mvrnorm(N, rep(0, ncol(Sigma)), Sigma)
X <- sweep(X, 2, colMeans(X), "-", FALSE)
colnames(X) <- paste("X", 1:ncol(X), sep = "")
man <- make_manifest(X, how = "mle", seed = NULL)

free <- is.na(fixed)
beta <- new("parameter.coef", x = fixed, free = free, num_free = sum(free))

free <- lower.tri(Phi)
Phi <- new("parameter.cormat", x = Phi, free = free, num_free = sum(free))

res.mle <- make_restrictions(man, beta = beta, Phi = Phi, discrepancy = "MLE")

starts <- apply(res.mle@Domains, 1, FUN = function(x) runif(1000, x[1], x[2]))
starts[1,] <- c(.5, .4, .3, Lambda[Lambda > 0], log(sqrt(diag(Sigma))))

opt.mle <- Factanal(man, res.mle, starting.values = starts, 
			max.generations = 200, seeds = NULL)
model_comparison(opt.mle)
if(TRUE) {
	library(sem)
	RAM <- FA2RAM(opt.mle)
	sem_CFA <- sem(RAM, model.matrix(man, standardized = FALSE), N)
	print(summary(sem_CFA))
	deviance(opt.mle)
	FAiR:::FAiR_stress_test(opt.mle)
}

starts[2,] <- opt.mle@optimization$extraction$par
res.adf <-  make_restrictions(man, beta = beta, Phi = Phi, discrepancy = "ADF")

opt.adf <- Factanal(man, res.adf, starting.values = starts, seeds = NULL)

summary(opt.adf)
model_comparison(opt.adf)
FAiR:::FAiR_stress_test(opt.adf)
