context("rank_shift")

test_that("rank_shift loads and returns correct result", {
    # Load our example data set
    data(knz_001d)
    #knz_001d <- read.csv(system.file("extdata", "knz_001d.csv", package="codyn"), sep=",", header=TRUE)
    expect_that(names(knz_001d)[4], equals("abundance"))


    
    # Basic test if mean rank produces data frame with right structure and values
    result <- rank_onerep(knz_001d, "year", "species", "abundance")
    expect_that(length(names(result)), equals(2))
    expect_that(names(result)[[1]], matches("year_pair"))
    expect_that(class(result[,1]), matches("factor"))
    expect_that(names(result)[[2]], matches("MRS"))
    expect_that(class(result[,2]), matches("numeric"))
    expect_that(nrow(result), equals(23))
    expect_that(result[[1,2]], equals(332.9921, tolerance=.00001))
    expect_that(result[[23,2]], equals(670.6388, tolerance=.00001))

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
  	names(dat2) <-  c("sp", "yr", "sub", "abund")
  	#make subplot a character
  	dat3 <- dat1
  	dat3$subplot <- as.character(dat3$subplot)

  	## rank_onerep() recognizes correct names

  	expect_error(rank_onerep(dat1, time.var = "pippin",
  	                      species.var = "species",
  	                      abundance.var = "abundance"))
  	expect_error(rank_onerep(dat1, time.var = "year",
  	                      species.var = "meriadoc",
  	                      abundance.var = "abundance"))
  	expect_error(rank_onerep(dat1, time.var = "year",
  	                      species.var = "species",
  	                      abundance.var = "samwise"))

	#test the get_slope function
  	myresults <- rank_onerep(dat1, time.var = "year",
  	                      species.var = "species",
  	                      abundance.var = "abundance")
  	expect_that(class(myresults[,2]), equals("numeric"))
  	expect_that(length(myresults), equals(2))
  	#test that rank_onerep function works with different column names
  	myresults2 <- rank_onerep(dat2,  time.var = "yr",
  	                       species.var = "sp",
  	                       abundance.var = "abund")
  	expect_that(myresults2, equals(myresults))
  	#test that gives a warning if running on factor instead of numeric
  	expect_error(rank_onerep(dat2, time.var = "yr",
  	                      species.var = "sp",
  	                      abundance.var = "subplot"))

	#test the mean_rank_shift function
  	#test that works on a single replicate
  	myresults3<-mean_rank_shift(dat1, replicate.var=NA, time.var="year", species.var="species", abundance.var="abundance")
  	expect_that(myresults3, equals(myresults2))

  	#test that works whether replicate is a character or factor
  	myresults4<-mean_rank_shift(dat3, replicate.var="subplot", time.var="year", species.var="species", abundance.var="abundance")
  	expect_that((myresults4[2,2]), equals(myresults3[2,2]))

  	#test that works with multiple replicates
  	myresults5<-mean_rank_shift(knz_001d, replicate.var="subplot", time.var="year", species.var="species", abundance.var="abundance")
  	expect_that(myresults4[2,2], equals(myresults5[2,2]))

    #test that works with different column names
  	myresults6<-mean_rank_shift(knz_001d2, replicate.var="sub", time.var="yr", species.var="sp", abundance.var="abund")
  	expect_that(myresults6[2,2], equals(myresults5[2,2]))

	#test that works regardless of whether parameter is specified or just ordered
  	myresults7<-mean_rank_shift(knz_001d,  "year", "species", "abundance", "subplot")
  	expect_that(myresults7, is_identical_to(myresults5))

	#test that works with different column orders if names specified
  	myresults8<-mean_rank_shift(knz_001d, abundance.var="abundance", replicate.var="subplot", species.var="species", time.var="year")
  	expect_that(myresults8, is_identical_to(myresults5))

  	#test that works with different column names
  	myresults9<-mean_rank_shift(knz_001d2, replicate.var="sub", time.var="yr", species.var="sp", abundance.var="abund")
  	expect_that(myresults9[2,2], equals(myresults5[2,2]))

	#test that it works even if there are additional unused columns
  	knz_001d3<-knz_001d
  	knz_001d3$site<-"KNZ"
  	myresults10<-mean_rank_shift(knz_001d3, "year", "species", "abundance", "subplot")
  	expect_that(myresults10, is_identical_to(myresults5))

  	#test that gives error when abundance column is a character or factor
  	expect_error(mean_rank_shift(knz_001d2, replicate.var="sub", time.var="yr", species.var="sp", abundance.var="randcharacter"))
  	expect_error(mean_rank_shift(knz_001d2, replicate.var="sub", time.var="yr", species.var="sp", abundance.var="randfactor"))

  	# test that works regardless of order of the input replicates
  	knz_001dreorder <-knz_001d[order(knz_001d$abundance, knz_001d$year, knz_001d$species),]
  	myresults<-mean_rank_shift(knz_001d, time.var="year", species.var="species", abundance.var="abundance",
  	                            replicate.var="subplot")

  	myresults_reorder<-mean_rank_shift(knz_001dreorder, time.var="year", species.var="species", abundance.var="abundance",
  	                            replicate.var="subplot")
  	expect_equal(myresults, myresults_reorder)

})