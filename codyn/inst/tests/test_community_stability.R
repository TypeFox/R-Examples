context("community_stability")

test_that("community_stability loads and returns correct result", {

  # Load our example data set
  data(knz_001d)

  expect_that(names(knz_001d)[4], equals("abundance"))

  #give new column names
  knz_001d2 <- knz_001d
  names(knz_001d2) <- c("sp", "yr", "sub", "abund")

  #add a random character and factor column
  knz_001d2$randcharacter <- "rchar"
  knz_001d2$randfactor <- as.factor(knz_001d2$randcharacter)

  #take a subset
  dat1 <- subset(knz_001d, knz_001d$subplot == "A_1")

  #rename the subset
  dat2 <- dat1
  names(dat2) <- c("sp", "yr", "sub", "abund")

  #make subplot a character
  dat3 <- dat1
  dat3$subplot <- as.character(dat3$subplot)

  #test the stability_onerep function
  dat1agg <- aggregate(abundance~year + subplot, data = dat1, sum)
  myresults <- codyn:::stability_onerep(dat1agg,  "abundance")
  expect_is(myresults, "numeric")
  expect_equal(length(myresults), 1)
  expect_equal(myresults, expected = 4.1233, tolerance = 0.00001)


  #test that works with different column names
  dat2agg <- aggregate(abund~yr + sub, data = dat2, sum)
  myresults2 <- codyn:::stability_onerep(dat2agg, "abund")
  expect_equivalent(myresults, myresults2)

  ## check that x cannnto be character
  expect_error(codyn:::stability_onerep(dat2agg, "sub"),
               "is not a numeric or integer")

  #test that errors if subplot does not exist
  expect_error(stability_onerep(dat2agg, "subplot"),
               "df does not have name subplot")


  #test the community_stability function
  #test that works on a single replicate
  myresults3 <- community_stability(dat1, replicate.var = NA,
                                    time.var = "year",
                                    abundance.var = "abundance")
  expect_equivalent(myresults3, myresults2)

  #test that will still run if there are missing levels in a factor "replicate"; deleting levels that are NaN
  myresults4 <- community_stability(dat1, replicate.var = "subplot",
                                    time.var = "year", abundance.var = "abundance")

  expect_equal(myresults3, myresults4$stability)

  #test that works whether replicate is a character or factor
  myresults6 <- community_stability(dat3,
                                    replicate.var = "subplot",
                                    time.var = "year",
                                    abundance.var = "abundance")

  expect_is(myresults4$subplot, "factor")
  expect_is(myresults6$subplot, "character")
  expect_equivalent(myresults6$stability, myresults4$stability)

  #test that works with multiple replicates
  myresults7 <- community_stability(knz_001d, replicate.var = "subplot",
                                    time.var = "year", abundance.var = "abundance")


  expect_equivalent(myresults4, myresults7[1,])

  #test that works with different column names
  myresults8 <- community_stability(knz_001d2,
                                    replicate.var = "sub",
                                    time.var = "yr", abundance.var = "abund")

  expect_equivalent(myresults7, myresults8)

  #test that gives error when abundance column is a character or factor
  expect_error(community_stability(knz_001d2, replicate.var = "sub",
                                   time.var = "yr", abundance.var = "randcharacter"))
  expect_error(community_stability(knz_001d2, replicate.var = "sub",
                                   time.var = "yr", abundance.var = "randfactor"))


  # test that works regardless of order of the input replicates
  knz_001dreorder <- knz_001d[order(knz_001d$abundance, knz_001d$year, knz_001d$species),]
  myresults_original <- community_stability(knz_001d, time.var  =  "year", abundance.var = "abundance",
                             replicate.var = "subplot")

  myresults_reorder <- community_stability(knz_001dreorder, time.var = "year",
                                         abundance.var = "abundance",
                                         replicate.var = "subplot")
  expect_equal(myresults_original, myresults_reorder)

  })