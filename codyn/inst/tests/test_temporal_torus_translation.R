context("temporal_torus_translation")

# set up testing data ----

# Load our example data set
data(knz_001d)

## check dataset
expect_that(names(knz_001d)[4], equals("abundance"))

#give new column names
knz_001d2 <- knz_001d
names(knz_001d2) <- c("sp", "yr", "sub", "abund")

#add a random character and factor column
knz_001d2$randcharacter <- sample(letters,
                                  size = nrow(knz_001d2),
                                  replace = TRUE)
knz_001d2$randfactor <- as.factor(knz_001d2$randcharacter)

#take a subset
dat1 <- subset(knz_001d, knz_001d$subplot == "A_1")

#rename the subset
dat2<-dat1
names(dat2)=c("sp", "yr", "sub", "abund")

#make subplot a character
dat3<-dat1
dat3$subplot<-as.character(dat3$subplot)

#take a two replicate subset
dat4<-subset(knz_001d, subplot=="A_1"|subplot=="A_2")


# tests -----------------------------------------------




test_that("utility functions work properly", {
  #make a species matrix

  datmat <- codyn:::transpose_community(dat1, "year", "species", "abundance")

  #Test the shuffle_community function
  myresults <- codyn:::shuffle_community(datmat)
  expect_is(myresults, "data.frame")

  #test that does not generate the same matrix each time
  myresults2 <- codyn:::shuffle_community(datmat)
  myresults3 <- sort(as.vector(myresults == myresults2))

  expect_true(any(as.vector(myresults == myresults2) == FALSE))

  #test that random matrix has same dimensions as datmat
  expect_equal(dim(myresults), dim(datmat))

  #test that myresults has the same values as datmat, if both are sorted
  rand.col.1 <- sort(as.matrix(myresults[,1]))
  dat.col.1 <- sort(as.matrix(datmat[,1]))
  expect_equal(rand.col.1, dat.col.1)
})

test_that("cyclic_shift returns correct result", {
  expect_error(cyclic_shift(dat1, time.var = "year", species.var = "species",
                            abundance.var = "subplot",
                            codyn:::variance_ratio_matrixdata,
                            bootnumber = 1),
               "Abundance variable is not numeric")

  # Test the cyclic_shift function
  myresults <- cyclic_shift(dat1, time.var = "year", species.var = "species",
                            abundance.var = "abundance",
                            FUN = codyn:::variance_ratio_matrixdata,
                            bootnumber = 1)

  expect_is(myresults, "cyclic_shift")

  #test structure of this class
  expect_equal(length(myresults), 1) ## will change soon

  expect_is(myresults$out, "numeric")
  #test that does not generate the same value every time
  myresults2 <- cyclic_shift(dat1, time.var = "year",
                             species.var = "species",
                             abundance.var = "abundance",
                             FUN = codyn:::variance_ratio_matrixdata,
                             bootnumber = 1)

  expect_false(myresults$out == myresults2$out)

  #test that is not sensitive to different column names
  myresults3 <- cyclic_shift(dat2, time.var = "yr",
                             species.var = "sp", abundance.var = "abund",
                             FUN = codyn:::variance_ratio_matrixdata,
                             bootnumber = 1)
  expect_equal(length(myresults3), 1)
  expect_is(myresults3$out, "numeric")
})

test_that("cyclic_shift works with replicates", {

  shift_with_rep <- cyclic_shift(knz_001d2,
                                 replicate.var = "sub", species.var = "sp",
                                 time.var = "yr",
                                 abundance.var = "abund",
                                 FUN = variance_ratio_matrixdata,
                                 bootnumber = 10)

  ## supposed to be a test for "it ran, and the output is the right length"
  expect_that(length(shift_with_rep$out), equals(10))

})

test_that("confint.cyclic_shift works as expected",{

  # Test the cyclic_shift function
  myresults_ <- cyclic_shift(dat1, time.var = "year", species.var = "species",
                            abundance.var = "abundance",
                            FUN = codyn:::variance_ratio_matrixdata,
                            bootnumber = 10)

  myresults_conf <- confint(myresults_)
  expect_is(myresults_conf, "data.frame")
  expect_named(myresults_conf, c("lowerCI", "upperCI", "nullmean"))
  expect_that(nrow(myresults_conf), equals(1))
})


  #   expect_error(confint.cyclic_shift(knz_001d2, replicate.var="sub", species.var="sp", time.var="yr", abundance.var="randfactor", method=variance_ratio_matrixdata, bootnumber=2, li=0.025, ui=0.975, average.replicates=FALSE))



  ##Test that missing levels in a factor replicate are dropped (so that the same number of subplots if run as a character
  #or as a factor with additional, missing levels)

  # #For example :
  # myresults2<-confint.cyclic_shift(dat1, replicate.var="subplot", species.var="species", time.var="year", abundance.var="abundance", method=variance_ratio_matrixdata, bootnumber=2, li=0.025, ui=0.975, average.replicates=FALSE)
  # expect_true(is.factor(myresults2$subplot))
  #
  # myresults2b<-confint.cyclic_shift(dat3, replicate.var="subplot", species.var="species", time.var="year", abundance.var="abundance", method=variance_ratio_matrixdata, bootnumber=2, li=0.025, ui=0.975, average.replicates=FALSE)
  # expect_true(is.character(myresults2b$subplot))
  #
  # expect_equal(as.character(myresults2$subplot), myresults2b$subplot)
  #
  # myresults2b<-confint.cyclic_shift(knz_001d, replicate.var="subplot", species.var="species", time.var="year", abundance.var="abundance", method=variance_ratio_matrixdata, bootnumber=2, li=0.025, ui=0.975, average.replicates=FALSE)
  #
  # #Test that is correct for whether "averagereps" is true or false
  # dat5<-dat4
  # dat5$subplot<-as.character(dat5$subplot)
  # myresults3<-confint.cyclic_shift(dat5, replicate.var="subplot", species.var="species", time.var="year", abundance.var="abundance", method=variance_ratio_matrixdata, bootnumber=2, li=0.025, ui=0.975, average.replicates=TRUE)
  # myresults4<-confint.cyclic_shift(dat5, replicate.var="subplot", species.var="species", time.var="year", abundance.var="abundance", method=variance_ratio_matrixdata, bootnumber=2, li=0.025, ui=0.975, average.replicates=FALSE)
  #
  # #expect that subplot names are the same between input and output
  # expect_equal(unique(dat5$subplot), myresults4$subplot)
  #
  # #Test both options for "averagereps"
  # myresults4<-confint.cyclic_shift(knz_001d, replicate.var="subplot", species.var="species", time.var="year", abundance.var="abundance", method=variance_ratio_matrixdata, bootnumber=2, li=0.025, ui=0.975, average.replicates=TRUE)
  # expect_that(nrow(myresults4), equals(1))
  # myresults5<-confint.cyclic_shift(knz_001d, replicate.var="subplot", species.var="species", time.var="year", abundance.var="abundance", method=variance_ratio_matrixdata, bootnumber=2, li=0.025, ui=0.975, average.replicates=FALSE)
  # expect_that(nrow(myresults5), equals(length(unique(knz_001d$subplot))))


