context("ccaR with coordinate list from RasterLayer")

test_that("ccaR works with cordinate list for landcover", {
	data(landcover)
	test <- cca(landcover, cell.class=1,s=0.01)
	expect_equivalent(length(test$size), 20)
	expect_equivalent(table(test$cluster[,3])[1], 2)
	test2 <- cca(landcover, s=2*res(landcover)[1]-0.00001, cell.class=c(2), compare="s")
	expect_equivalent(length(test2$size), 60)
	expect_equivalent(table(test2$cluster[,3])[1], 41)
	test3 <- cca(landcover, s=2*res(landcover)[1]-0.00001, cell.class=c(0,1), compare="")
	expect_equivalent(length(test3$size), 60)
	expect_equivalent(table(test3$cluster[,3])[1], 41)
})
