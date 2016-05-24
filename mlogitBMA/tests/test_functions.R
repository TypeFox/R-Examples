library(mlogitBMA)

test.specification <- function() {
	spec <- mnl.spec(choice ~ fuel + price + cost | hsg2, Car,
					varying=5:ncol(Car), sep='')
	stopifnot(all(dim(spec$variable.used) == c(6,4)))
	stopifnot(all(spec$varying.names == c('fuel', 'price', 'cost')))
	stopifnot(sum(spec$same.coefs) == 3 & !spec$same.coefs['hsg2'])
	stopifnot(sum(spec$intercepts)==5 &  !spec$intercepts[1])
	cat('\nSpecification test OK.\n')
}

test.estimation.car <- function() {
	cat('\nRunning test for MNL estimation on Car dataset ... ')
	est <- estimate.mlogit(choice ~ price + cost | coml5, Car, 
						varying=5:ncol(Car), sep='')
	sest <- summary(est)
	stopifnot(all(dim(sest$coefs) == c(12, 4)))
	stopifnot(sest$lratio > 0.1)
	cat(' OK.\n')
}

test.estimation.catsup <- function() {
	cat('\nRunning test for MNL estimation on Catsup dataset ... ')
	est <- estimate.mlogit(choice ~ disp + feat + price, Catsup, 
						varying=2:(ncol(Catsup)-1), sep='.')
	sest <- summary(est)
	stopifnot(all(dim(sest$coefs) == c(6, 4)))
	stopifnot(sest$lratio > 0.3)
	stopifnot(sest$bic > 5083)
	print(sest)
	cat(' OK.\n')
}

test.bic.mlogit.car <- function() {
	cat('\nRunning test for BMA on Car dataset ... ')
	res <- bic.mlogit(choice ~ price + cost  + speed + acc + size | hsg2, Car, 
						varying=5:ncol(Car), sep='', include.intercepts = FALSE, 
						verbose = TRUE)
	stopifnot(all(dim(res$bic.glm$which)==c(1,14))) # 1 model selected, 14 variables in total
	cat('... BMA test OK.\n')
}

test.bic.mlogit.catsup <- function() {
	cat('\nRunning test for BMA on Catsup dataset ... ')
	res <- bic.mlogit(choice ~ 1 | disp + feat + price, Catsup, 
						varying=2:(ncol(Catsup)-1), sep='.', 
						base.choice = 4, 
						include.intercepts = FALSE, 
						verbose = TRUE)
	summary(res)
	stopifnot(all(dim(res$bic.glm$which)==c(2,11))) # 2 models selected, 11 variables in total
	cat('... BMA test OK.\n')
	est.res <- estimate.mlogit(res, Catsup)
	stopifnot(length(est.res)==2)
	stopifnot(all(dim(est.res[[1]]$coefs) == c(12, 4)))
	stopifnot(all(dim(est.res[[2]]$coefs) == c(11, 4)))
	cat('Estimation on the BMA object OK.\n')
}

# load data
data('Car', package='mlogit')
# convert the choice column into a numerical code,
# since it is the way the alternative-specific variables are constructed
Car[,'choice'] <- as.integer(gsub('^choice', '', Car[,'choice']))
data('Catsup', package='mlogit')

test.specification()
test.estimation.car()
test.estimation.catsup()
test.bic.mlogit.car()
# turn off to speed the tests up
#test.bic.mlogit.catsup()