## ----setup, echo = FALSE, message = FALSE--------------------------------
library(knitr);
opts_chunk$set(cache = TRUE, autodep = TRUE);

## ----install, echo = TRUE, eval = FALSE----------------------------------
#  install.packages("NNLM");
#  
#  # or get a dev-version from Github
#  library(devtools);
#  install_github('linxihui/NNLM')

## ----interface, eval = FALSE---------------------------------------------
#  nnlm(x, y, alpha = rep(0,3), method = c('scd', 'lee'), loss = c('mse', 'mkl'), init = NULL,
#  	mask = NULL, check.x = TRUE, max.iter = 10000L, rel.tol = 1e-12, n.threads = 1L)
#  
#  nnmf(A, k = 1L, alpha = rep(0,3), beta = rep(0,3), method = c("scd", "lee"), loss = c("mse", "mkl"),
#  	init = NULL, mask = NULL, W.norm = -1L, check.k = TRUE, max.iter = 500L, rel.tol = 1e-4,
#  	n.threads = 1L, trace = 10L, verbose = 1L, show.warning = TRUE,
#  	inner.max.iter = ifelse('mse' == loss, 50L, 1L), inner.rel.tol = 1e-09)

## ----alg-compare,  warning = FALSE---------------------------------------
library(NNLM);
set.seed(123);

k <- 15;
init <- list(W = matrix(runif(nrow(nsclc)*k), ncol = k),
	H = matrix(runif(ncol(nsclc)*k), nrow = k));
scd.mse  <- nnmf(nsclc, k, init = init, max.iter = 100, rel.tol = -1);
lee.mse  <- nnmf(nsclc, k, init = init, max.iter = 100, rel.tol = -1, method = 'lee');
scd.mkl  <- nnmf(nsclc, k, init = init, max.iter = 5000, rel.tol = -1, loss = 'mkl');
lee.mkl  <- nnmf(nsclc, k, init = init, max.iter = 5000, rel.tol = -1, loss = 'mkl',
	method = 'lee');
lee.mse1 <- nnmf(nsclc, k, init = init, max.iter = 5000, rel.tol = -1, method = 'lee',
	inner.max.iter = 1);

## ----alg-comp-plot, fig.show = 'hold'------------------------------------
plot(NULL, xlim = c(1, 3000), ylim = c(0.15, 0.45), xlab = 'Epochs', ylab = 'MSE');
lines(cumsum(scd.mse$average.epochs), scd.mse$mse, col = 'firebrick1');
lines(cumsum(lee.mse$average.epochs), lee.mse$mse, col = 'orange');
lines(cumsum(scd.mkl$average.epochs), scd.mkl$mse, col = 'chartreuse3');
lines(cumsum(lee.mkl$average.epochs), lee.mkl$mse, col = 'deepskyblue4');
lines(cumsum(lee.mse1$average.epochs), lee.mse1$mse, col = 'orange', lty = 2);
legend('topright', bty = 'n', lwd = 1, lty = c(1,1,2,1,1),
	legend = c('SCD-MSE', 'LEE-MSE', 'LEE-MSE-1', 'SCD-MKL', 'LEE-MKL'),
	col = c('firebrick1', 'orange', 'orange', 'chartreuse3', 'deepskyblue4'));

plot(NULL, xlim = c(1, 3000), ylim = c(0.01, 0.034), xlab = 'Epochs', ylab = 'MKL');
lines(cumsum(scd.mse$average.epochs), scd.mse$mkl, col = 'firebrick1');
lines(cumsum(lee.mse$average.epochs), lee.mse$mkl, col = 'orange');
lines(cumsum(scd.mkl$average.epochs), scd.mkl$mkl, col = 'chartreuse3');
lines(cumsum(lee.mkl$average.epochs), lee.mkl$mkl, col = 'deepskyblue4');
lines(cumsum(lee.mse1$average.epochs), lee.mse1$mkl, col = 'orange', lty = 2);
legend('topright', bty = 'n', lwd = 1, lty = c(1,1,2,1,1),
	legend = c('SCD-MSE', 'LEE-MSE', 'LEE-MSE-1', 'SCD-MKL', 'LEE-MKL'),
	col = c('firebrick1', 'orange', 'orange', 'chartreuse3', 'deepskyblue4'));

summary.nnmf <- function(x) {
	if (x$n.iteration < 2) {
		rel.tol <- NA_real_;
	} else {
		err <- tail(x$target.loss, 2);
		rel.tol <- diff(err)/mean(err); 
		}
	return(c(
		'MSE' = tail(x$mse, 1), 'MKL' = tail(x$mkl, 1), 'Target' = tail(x$target.loss, 1),
		'Rel. tol.' = abs(rel.tol), 'Total epochs' = sum(x$average.epochs),
		'# Interation' = x$n.iteration, x$run.time[1:3]));
	};

library(knitr);
kable(
	sapply(
		X = list(
			'SCD-MSE' = scd.mse,
			'LEE-MSE' = lee.mse,
			'LEE-MSE-1' = lee.mse1,
			'SCD-MKL' = scd.mkl,
			'LEE-MKL' = lee.mkl
			),
		FUN = function(x) {z <- summary(x); sapply(z, sprintf, fmt = '%.4g')}
		),
	align = rep('r', 5),
	caption = 'Compare performance of different algorithms'
	);

## ----nsclc, message = TRUE-----------------------------------------------
library(NNLM);
set.seed(123);

data(nsclc, package = 'NNLM')
str(nsclc)

decomp <- nnmf(nsclc[, 1:80], 3, rel.tol = 1e-5);
decomp

## ----nsclc1, fig.show = 'hold', message = TRUE---------------------------
heatmap(decomp$W, Colv = NA, xlab = 'Meta-gene', ylab = 'Gene', margins = c(2,2),
	labRow = '', labCol = '', scale = 'column', col = cm.colors(100));
heatmap(decomp$H, Rowv = NA, ylab = 'Meta-gene', xlab = 'Patient', margins = c(2,2),
	labRow = '', labCol = '', scale = 'row', col = cm.colors(100));

## ----nsclc2, fig.show = 'hold', fig.width = 5,fig.height = 5, fig.align = 'center'----
par(mar = c(5,5,4,5))
plot(decomp$mse, type = 'l', lwd = 2, xlab = 'Iteration', ylab = '',
	col = 'deepskyblue4', axes = FALSE, ylim = c(0.35, 0.45));
axis(1, ylim = seq_along(decomp$mse));
axis(2, ylim = c(0.35, 0.45), col = 'deepskyblue4', col.axis = 'deepskyblue4');
mtext("Mean Square Error", side = 2, col = 'deepskyblue4', line = 3)
par(new = TRUE);
plot(decomp$mkl, type = 'l', lwd = 2, xlab = '', ylab = '', col = 'firebrick4',
	axes = FALSE, ylim = c(0.025, 0.035));
axis(4, ylim = c(0.025, 0.035), col = 'firebrick4', col.axis = 'firebrick4');
mtext("Mean KL-divergence", side = 4, col = 'firebrick4', line = 3)

## ----nsclc3, fig.show = 'hold', message = TRUE---------------------------
# find the expressions of meta-genes for patient 81-100
newH <- predict(decomp, nsclc[, 81:100], which = 'H');
str(newH)

## ----isopure-data, fig.align = 'center'----------------------------------
set.seed(123);
if (!require(ISOpureR)) {
	install.packages('ISOpureR', repos = "http://cran.stat.sfu.ca/");
	library(ISOpureR);
	}

path.to.data <- file.path(system.file(package = 'ISOpureR'), 'extdata/Beer');
# normal profile
load(file.path(path.to.data, 'beer.normaldata.1000.transcripts.RData'));
# transcriptome of 30 patients (part of the Beer dataset)
load(file.path(path.to.data, 'beer.tumordata.1000.transcripts.30.patients.RData'));

# assume k = 3, beer.normaldata is the known healthy profile
beer.nmf <- nnmf(beer.tumordata, k = 3, init = list(W0 = beer.normaldata));

# compute proportion of tumour content
tm.frac.nmf <- with(beer.nmf, colSums(W[,1:3] %*% H[1:3,])/colSums(W %*% H));

# tumour content from ISOpureR using the full dataset
tm.frac.full <- read.delim(file.path(path.to.data, "alphapurities_full_dataset.txt"));

## ----isopure-plot, fig.align = 'center', fig.width = 5, fig.height = 5----
plot(tm.frac.full$alphapurities*100, tm.frac.nmf*100,
	xlim = c(20, 100), ylim = c(20, 100), col = 'firebrick1',
	xlab = "% tumour from ISOpureR full data",
	ylab = "% tumour from NMF partial data");
abline(a = 0, b = 1, lty = 2)

## ----missing-data, fig.show = 'hold', message = FALSE--------------------
set.seed(123);
nsclc2 <- nsclc;
index <- sample(length(nsclc2), length(nsclc2)*0.3);
nsclc2[index] <- NA;

## ----missing-nmf, results = 'hold', fig.show = 'hold', message = FALSE----
# NMF imputation
system.time(nsclc2.nmf <- nnmf(nsclc2, 2));
nsclc2.hat.nmf <- with(nsclc2.nmf, W %*% H);

## ----missing-mic, results = 'hold', fig.show = 'hold', message = FALSE----
# multivariate imputation by chained equations (MICE)
if(!require(mice)) {
	install.packages('mice', repos = "http://cran.stat.sfu.ca/");
	library(mice);
	}

# logarithm for positivity and the normality assumption
system.time(nsclc2.mice <- mice(log(nsclc2), m = 1, printFlag = FALSE)); # one imputation
nsclc2.hat.mice <- exp(as.matrix(complete(nsclc2.mice)));

## ----missing-missforest, results = 'hold', fig.show = 'hold', message = FALSE----
# imputation using random forest
if(!require(missForest)) {
	install.packages('missForest', repos = "http://cran.stat.sfu.ca/");
	library(missForest);
	}

system.time(capture.output(nsclc2.missForest <- missForest(log(nsclc2))));
nsclc2.hat.missForest <- exp(nsclc2.missForest$ximp);

## ----missing-plot, fig.show = 'hold'-------------------------------------
# simple imputation, fill-in median expression values of genes
nsclc2.hat.median <- matrix(
	apply(nsclc2, 1, median, na.rm = TRUE),
	nrow = nrow(nsclc2), ncol = ncol(nsclc2)
	);

# compare different imputations
library(knitr);
kable(
	sapply(
		X = list(
			Baseline = mean(nsclc2, na.rm = TRUE),
			Medians = nsclc2.hat.median[index],
			MICE = nsclc2.hat.mice[index],
			MissForest = nsclc2.hat.missForest[index],
			NMF = nsclc2.hat.nmf[index]
			),
		FUN = mse.mkl, # mean square error, mean KL-divergence
		obs = nsclc[index]
		),
	caption = "A comparison of different imputation methods",
	digits = 4
	);

plot(nsclc[index], nsclc2.hat.median[index], col = 'chartreuse3', cex = 0.3, pch = 20,
	ylim = c(2, 14), xlab = "Observed", ylab = "Medians")
abline(a = 0, b = 1, lwd = 2, lty = 2)

plot(nsclc[index], nsclc2.hat.mice[index], col = 'orange', cex = 0.3, pch = 20,
	ylim = c(2, 14), xlab = "Observed", ylab = "MICE imputed")
abline(a = 0, b = 1, lwd = 2, lty = 2)

plot(nsclc[index], nsclc2.hat.missForest[index], col = 'darkorchid', cex = 0.3, pch = 20,
	ylim = c(2, 14), xlab = "Observed", ylab = "MissForest imputed")
abline(a = 0, b = 1, lwd = 2, lty = 2)

plot(nsclc[index], nsclc2.hat.nmf[index], col = 'deepskyblue', cex = 0.3, pch = 20,
	ylim = c(2, 14), xlab = "Observed", ylab = "NMF imputed")
abline(a = 0, b = 1, lwd = 2, lty = 2)

## ----rank-sim, fig.width = 5, fig.height = 5, fig.align = 'center', warning = FALSE----
set.seed(678);

n <- 400;
m <- 50;
k <- 3;
W  <- matrix(runif(n*k),  n, k);
H  <- matrix(10*runif(k*m),  k, m);
noise <- matrix(rnorm(n*m), n, m);
A <- W %*% H + noise;
A[A < 0] <- 0;

plot(-1, xlim = c(1,6), ylim = c(0.5, 2.5), xlab = "Rank", ylab = "MSE")
cols <- c('deepskyblue', 'orange', 'firebrick1', 'chartreuse3');
for (col in cols) {
	ind <- sample(length(A), 0.3*length(A));
	A2 <- A;
	A2[ind] <- NA;
	err <- sapply(X = 1:6,
		FUN = function(k) {
			z <- nnmf(A2, k);
			c(mean((with(z, W %*% H)[ind] - A[ind])^2), tail(z$mse, 1));
			}
		);
	invisible(lines(err[1,], col = col, type = 'b'));
	invisible(lines(err[2,], col = col, type = 'b', lty = 2));
	}

## ----rank-nsclc, fig.width = 5, fig.height = 5, fig.align = 'center', warning = FALSE----
set.seed(567);

plot(0, xlim = c(1,10), ylim = c(0.4, 1.4), xlab = "Rank", ylab = "MSE")
cols <- c('deepskyblue', 'orange', 'firebrick1', 'chartreuse3');
for (col in cols) {
	index2 <- sample(which(!is.na(nsclc2)), 2000);
	nsclc3 <- nsclc2;
	nsclc3[index2] <- NA;
	err <- sapply(X = 1:10,
		FUN = function(k, A) {
			z <- nnmf(A, k, verbose = FALSE);
			mean((with(z, W%*%H)[index2] - nsclc2[index2])^2)
			},
		A = nsclc3
		);
	invisible(lines(err, col = col, type='b', lwd = 2, cex = 1));
	}

## ----deconvol, message = TRUE--------------------------------------------
n <- 200; m <- 50;
k <- 5; k1 <- 2; k2 <- 1;

set.seed(123);
W  <- matrix(runif(n*k),  n, k); # unknown heterogeneous cancer profile
H  <- matrix(runif(k*m),  k, m);
W0 <- matrix(runif(n*k1), n, k1); # known healthy profile
H1 <- matrix(runif(k1*m), k1, m);
W1 <- matrix(runif(n*k2), n, k2); # unknown common cancer profile
H0 <- matrix(1, k2, m);
noise <- 0.1*matrix(runif(n*m), n, m);

# A is the observed profile to be de-convoluted
A <- W %*% H + W0 %*% H1 + W1 %*% H0 + noise;

deconvol <- nnmf(A, k = 5, init = list(W0 = W0, H0 = H0), method = 'lee');

## ----deconvol2, message = TRUE-------------------------------------------
round(cor(W, deconvol$W[, seq_len(k)]), 2);
# round(cor(t(H), t(deconvol$H[seq_len(k), ])), 2);

## ----deconvol3, eval = c(T, F, T), echo = c(TRUE, FALSE, TRUE),  message = TRUE----
permutation <- c(3, 1, 3, 5, 4);
round(cor(W, deconvol$W[, permutation]), 2);
#round(cor(t(H), t(deconvol$H[permutation, ])), 2);

## ----deconvol4, message = TRUE-------------------------------------------
round(cor(t(H1)), 2);
round(cor(t(H1), t(deconvol$H[6:7, ])), 2);
round(cor(W1, deconvol$W[, 8]), 2);

