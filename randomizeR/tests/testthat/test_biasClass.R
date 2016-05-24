###################################################################
# ----------------------------------------------------------------#
# Tests for objects from the bias class and associated functions  #
# ----------------------------------------------------------------#
###################################################################
context("Biases")


test_that("corGuess returns valid object", {
    expect_is(corGuess("CS"), "issue")
    expect_is(corGuess("CS"), "corGuess")
    csError <- corGuess("CS")
    csError@type <- "david"
    expect_error(print(csError))
  }
)

test_that("tests for correctness of corGuess", {
  randSeq <- genSeq(rpbrPar(N = 12, rb = 2))
  cs <- corGuess("CS")
  ds <- corGuess("DS")
  expect_true(sum(getCorGuesses(randSeq, cs) == "nG") == 6)
  expect_true(sum(getCorGuesses(randSeq, ds) == "nG") == 6)
  expect_true(sum(getCorGuesses(randSeq, cs) == getRandList(randSeq)) == 6)
  expect_true(sum(getCorGuesses(randSeq, ds) == getRandList(randSeq)) == 0)
  }
)


test_that("selBias returns valid object", {
    expect_is(selBias("CS", 2, "exact"), "selBias")
    csError <- corGuess("CS")
    csError@type <- "david"
    expect_error(print(csError))
    expect_error(selBias("CS", "a", "exact"))
  }
)


test_that("chronBias returns valid object", {
  expect_is(chronBias("linT", 2, "exact"), "chronBias")
  expect_is(chronBias("logT", 2, "exact"), "chronBias")
  expect_error(chronBias("stepT", 2, "exact"))
  expect_is(chronBias("stepT", 2, "exact", saltus = 2), "chronBias")
  }
)
