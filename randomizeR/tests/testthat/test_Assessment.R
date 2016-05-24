###################################################################
# ----------------------------------------------------------------#
# Tests for the assess function                                   #
# ----------------------------------------------------------------#
###################################################################

context("Assessment")

test_that("assess returns valid object", {
	seqs <- getAllSeq(pbrPar(bc = c(2, 2, 2)))
		
	type <- sample(c("CS", "DS"), 1)
	i1 <- corGuess(type)
	method <- sample(c("sim", "exact"), 1)
	i2 <- selBias(type, 1, method)
	# endp <- normEndp(c(0,0), c(1,1))
	expect_is(assess(seqs, i1), "assessment")
	expect_is(assess(seqs, i1, i2, endp = normEndp(c(0,0), c(1,1))), "assessment")
	
	type <- sample(c("imb", "absImb", "loss", "maxImb"),1)
	i3 <- imbal(type)
	expect_is(assess(seqs, i3), "assessment")

	expect_error(assess(seqs))
	expect_error(assess(seqs, "blubb"))
	expect_error(assess(seqs, "issue"))
	expect_error(assess(seqs, 42))
	
})

test_that("issues are computed right", {
  RP1 <- getAllSeq(pbrPar(bc = c(2, 2, 2, 2)))
  RP2 <- getAllSeq(bsdPar(10, 1))
  RP3 <- genSeq(rpbrPar(N = 8, rb = 2), r = 100)
  RP4 <- getAllSeq(mpPar(10,1))
  i1 <- corGuess("CS")
  i2 <- corGuess("DS")
  i3 <- imbal("imb")
  i4<- imbal("absImb")
  i5 <- imbal("loss")
  i6 <- imbal("maxImb")
  
  # test for all RPs if issues are calculated exact
  for (RP in c("RP1", "RP2", "RP3", "RP4")) {
    RP <- eval(parse(text = RP[1]))
    expect_true(all(assess(RP, i1, i2, i3, i4, i5, i6)@D[,3] == 0.75))
    expect_true(all(assess(RP, i1, i2, i3, i4, i5, i6)@D[,4] == 0.25))
    expect_true(all(assess(RP, i1, i2, i3, i4, i5, i6)@D[,5] == 0))
    expect_true(all(assess(RP, i1, i2, i3, i4, i5, i6)@D[,6] == 0))
    expect_true(all(assess(RP, i1, i2, i3, i4, i5, i6)@D[,7] == 0))
    expect_true(all(assess(RP, i1, i2, i3, i4, i5, i6)@D[,8] == 1))
  }  
})