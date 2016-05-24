library(rcbalance)
library(rcbsubset)
library(optmatch)
context('Precision for distances and penalties')

#skip_on_cran()

data(nuclearplants)
#reorder nuclearplants dataframe so treated units come first
nuclearplants <- nuclearplants[order(nuclearplants$pr, decreasing = TRUE),]

my.dist.struct <- build.dist.struct(z = nuclearplants$pr, 
	X = subset(nuclearplants[c('date', 't1', 't2', 'cap','cum.n')]), 
	exact = nuclearplants$ne, calip.option = 'none')
	
match.out <- rcbsubset(my.dist.struct, treated.info = nuclearplants[1:10,], control.info = nuclearplants[11:32,], fb.list = list('ct',c('ct','bw')))	

test_that('Fine balance is achieved', {
	expect_equal(match.out$fb.tables[[1]][,1], match.out$fb.tables[[1]][,2])
	expect_equal(match.out$fb.tables[[2]][,1], match.out$fb.tables[[2]][,2])
})

test_that('Overly large distances are recognized and caught', {
	new.dist <- matrix(c(1:4)*1e-5, nrow=2, ncol = 2) + .Machine$integer.max - 2
	expect_warning(expect_error(rcbalance(new.dist), 'Integer overflow in penalties!  Run with a higher tolerance, a lower penalty value, or fewer levels of fine balance.'), 'NAs introduced by coercion to integer range')
})






