library(rcbalance)
library(rcbsubset)
library(optmatch)
context('Excluding treated units')

data(nuclearplants)

#reorder nuclearplants dataframe so treated units come first
nuclearplants <- nuclearplants[order(nuclearplants$pr, decreasing = TRUE),]

my.dist.struct <- build.dist.struct(z = nuclearplants$pr, 
	X = subset(nuclearplants[c('date','t1','t2','cap','bw','ne')]), 
	exact = nuclearplants$cum.n, calip.option = 'none')

match.out <- rcbsubset(my.dist.struct)$matches


test_that('Distance object has empty rows', {
	expect_equal(length(my.dist.struct[[6]]),0)
})


test_that('Correct number excluded', {
		match.tab <- table(nuclearplants$cum.n, nuclearplants$pr)
		excl.count <- sum(pmax(match.tab[,2] - match.tab[,1],0))
		expect_equal(sum(nuclearplants$pr) - nrow(match.out), excl.count)
})

#match with only a single treated unit
dist.struct2 <- build.dist.struct(z = nuclearplants$pr[10:32], X = nuclearplants[10:32,c(3:5,8:9)], calip.option = 'none')

test_that('Matches with one pair are returned correctly', {
	match.out <- rcbsubset(dist.struct2, fb.list = list('ct'), treated.info = nuclearplants[10,], control.info = nuclearplants[nuclearplants$pr == 0,])
	expect_equal(rownames(match.out$matches), '1')
	expect_equal(dim(match.out$fb.tables[[1]]), c(1,2))
})