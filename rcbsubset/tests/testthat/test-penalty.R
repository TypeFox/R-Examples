library(rcbalance)
library(rcbsubset)
library(optmatch)
context('Setting network penalties')

data(nuclearplants)

#reorder nuclearplants dataframe so treated units come first
nuclearplants <- nuclearplants[order(nuclearplants$pr, decreasing = TRUE),]

my.dist.struct <- build.dist.struct(z = nuclearplants$pr, 
	X = subset(nuclearplants[c('date','t1','t2','cap','bw','cum.n')]), 
	exact = nuclearplants$ne, calip.option = 'none')

#create network with fine balance constraints	
match.net <- dist2net(my.dist.struct, k =1, ncontrol = sum(nuclearplants$pr == 0))	
match.net <- add.layer(match.net, nuclearplants$ct)
match.net <- add.layer(match.net, paste(nuclearplants$ct, nuclearplants$bw, sep = '.'))


test_that('penalty.update changes appropriate part of matching network', {
	p.orig <- match.net$p
	theta.orig <- match.net$theta
	match.net2 <- penalty.update(match.net, newtheta = 2*theta.orig, newp = 2*p.orig)
	cost.objs <- which(names(match.net) %in% c('cost', 'penalties', 'theta', 'p'))
	expect_identical(match.net[-cost.objs], match.net2[-cost.objs])
	#treated-control edges should be unaffected
	expect_equal(match.net$cost[1:match.net$tcarcs], 
	  match.net2$cost[1:match.net2$tcarcs])
})

test_that('penalize.near.exact changes appropriate part of matching network', {
	match.net2 <- penalize.near.exact(match.net, nuclearplants$pt)
	cost.objs <- which(names(match.net) %in% c('cost', 'penalties', 'theta', 'p'))
	expect_identical(match.net[-cost.objs], match.net2[-cost.objs])
	#only treated-control edges should be unaffected
	expect_equal(match.net$cost[-c(1:match.net$tcarcs)], 
	  match.net2$cost[-c(1:match.net2$tcarcs)])		
})

