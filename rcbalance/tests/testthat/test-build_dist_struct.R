library(rcbalance)
library(optmatch)
context('Building distance objects')

data(nuclearplants)

#reorder nuclearplants dataframe so treated units come first
nuclearplants <- nuclearplants[order(nuclearplants$pr, decreasing = TRUE),]

my.dist.struct <- build.dist.struct(z = nuclearplants$pr, 
	X = subset(nuclearplants[c('date', 't1', 't2', 'cap', 'bw', 'cum.n')]), 
	exact = nuclearplants$ne)	


test_that('Distance has expected number of pairings', {
	#sparsity via blocks
	dist.nocalip <- build.dist.struct(z = nuclearplants$pr, 
	X = subset(nuclearplants[c('date', 't1', 't2', 'cap', 'bw', 'cum.n')]), 
	exact = nuclearplants$ne, calip.option = 'none')
	expect_equal(length(unlist(dist.nocalip)), 
		count.pairings(z = nuclearplants$pr, exact = nuclearplants$ne))
		
	#sparsity via calipers
	dist.calip.only <- build.dist.struct(z = nuclearplants$pr, 
		X = subset(nuclearplants[c('date', 't1', 't2', 'cap', 'bw', 'cum.n')]))
	prop.score <- predict(glm(pr ~ date + t1 + t2 + cap + bw + cum.n, 
		data = nuclearplants, family = binomial()))
	calip.mat <- pmax(abs(outer(prop.score[nuclearplants$pr ==1], 
		prop.score[nuclearplants$pr == 0], FUN = "-")) - 0.2*sd(prop.score),0)
	calip.mat[calip.mat > 0] <- Inf
	expect_equal(length(unlist(dist.calip.only)), sum(is.finite(calip.mat)))
})		


test_that('Constant variables are removed', {
	nuclearplants$constant <- rep(1, nrow(nuclearplants))
	expect_identical(my.dist.struct, build.dist.struct(z = nuclearplants$pr, 
	X = subset(nuclearplants[c('date', 't1', 't2', 'cap', 'bw', 'cum.n', 'constant')]), 
	exact = nuclearplants$ne)) 
})

nuclearplants$not_bw <- 1- nuclearplants$bw
nuclearplants$bw_factor <- factor(nuclearplants$bw + 2*nuclearplants$not_bw)

test_that('Factors are converted properly to dummies', {
	expect_equal(
	  build.dist.struct(z = nuclearplants$pr, X = subset(nuclearplants[c('date', 't1', 't2', 'cap', 'bw', 'not_bw', 'cum.n')]), exact = nuclearplants$ne), 
	  build.dist.struct(z = nuclearplants$pr,X = subset(nuclearplants[c('date', 't1', 't2', 'cap', 'bw_factor', 'cum.n')]), exact = nuclearplants$ne)
	)
})

test_that('Missing values are handled properly', {
	nuclearplants$t2[1:4] <- NA
	nuclearplants$t2_missing <- is.na(nuclearplants$t2)
	nuclearplants$t2_impute <- nuclearplants$t2
	nuclearplants$t2_impute[nuclearplants$t2_missing] <- mean(nuclearplants$t2, na.rm = TRUE)
	nuclearplants$bw_factor[1:4] <- NA
	nuclearplants$bw_factor_na <- addNA(nuclearplants$bw_factor)
	
	expect_equal(
	  build.dist.struct(z = nuclearplants$pr, X = subset(nuclearplants[c('date', 't1', 't2', 'cap', 'bw_factor', 'cum.n')]), exact = nuclearplants$ne), 
	  build.dist.struct(z = nuclearplants$pr,X = subset(nuclearplants[c('date', 't1', 't2_impute','t2_missing', 'cap', 'bw_factor_na', 'cum.n')]), exact = nuclearplants$ne)	
	)
})


