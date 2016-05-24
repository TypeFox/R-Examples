## library(growcurves)

##
## perform tests of XZcov function
##


context("XZcov returns correct X and Z with accurate labels")

## Two treatment levels, {0,1}, and no nuisance covariates
data(datsim)

## extract labels
lab  <- relabel(label.input = datsim$trt, start = 0)

n.fix.nn <- function(x){
		out 	<- XZcov(time = datsim$time , trt = datsim$trt, trt.lab = lab$labeli.u, subject = datsim$subject, n.random = datsim$n.random, 
		n.fix_degree = x, formula = NULL, 
		random.only = FALSE, data = NULL)
		return(out)
		}

test_that("XZcov removes collinear covariates", {

	out <- suppressWarnings(n.fix.nn(3))

	## expect_that(n.fix.nn(3), gives_warning()) ## warning indicating nuisance covariates dropped - not used in testing bc will fail if warnings are suppressed.
	expect_that(ncol(out$X), equals(5)) ## drops nuisance 3rd order covariates

})

test_that("Correct dimensions and labels returned for X and Z", {

	out <- n.fix.nn(2)

	expect_that(colnames(out$X)[3], matches("trt_1")) ## Smallest numeric trt value used as baseline for X
	expect_that(is.null(out$X.n), is_true()) ## nuisance covariate matrix is null
	expect_that(ncol(out$X), equals(5))
	expect_that(ncol(out$Z), equals(3))
	expect_that(is.null(out$Z.n), is_true())
	expect_that(colnames(out$Z)[2], matches("time"))
	expect_that(is.null(out$y), is_true()) ## no response input to XZcov

})

test_that("Correct dimensions for X and Z with nuisance covariates",{

	data(datsimcov)
	out 	<- XZcov(time = datsimcov$time , trt = datsimcov$trt, trt.lab = lab$labeli.u, subject = datsimcov$subject, n.random = datsimcov$n.random, 
		n.fix_degree = 2, formula = datsimcov$formula, 
		random.only = FALSE, data = datsimcov$data)

	expect_that(!is.null(out$X.n),is_true())
	expect_that(!is.null(out$y),is_true())
	expect_that(ncol(out$X), equals(7)) ## account for 2 nuisance covariates
	expect_that(nrow(out$X), equals(length(out$y)))

})

		
