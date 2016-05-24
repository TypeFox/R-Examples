context("variance ratio")


# prepare data ----------------------------------------


# Load our example data set
data("knz_001d", package = "codyn")  # This doesn't work for CSV files :(
expect_that(names(knz_001d)[4], equals("abundance"))


#give new column names
knz_001d2 <- knz_001d
names(knz_001d2) = c("sp", "yr", "sub", "abund")

#take a subset
dat1 <- subset(knz_001d, knz_001d$subplot == "A_1")

#rename the subset
dat2 <- dat1
names(dat2) = c("sp", "yr", "sub", "abund")

#make subplot a character
dat3 <- dat1
dat3$subplot <- as.character(dat3$subplot)

#take a two replicate subset
dat4 <- subset(knz_001d, subplot == "A_1" | subplot == "A_2")

#make a species matrix
datmat <- codyn:::transpose_community(dat1, "year",  "species", "abundance")


# run tests -------------------------------------------


test_that("variance_ratio function returns correct result", {

    #test the class returned with default settings
    myresults <- variance_ratio(knz_001d, time.var = "year",
                              species.var = "species",
                              abundance.var = "abundance",
                              bootnumber = 1,
                              replicate.var = "subplot")


    expect_is(myresults, "data.frame")
    expect_equal(nrow(myresults), 1)
    expect_equal(myresults$VR, 1.01443, tolerance = 0.00001)

    expect_warning(variance_ratio(knz_001d, time.var = "year",
                                  species.var = "species",
                                  abundance.var = "abundance",
                                  bootnumber = 1,
                                  replicate.var = "subplot",
                                  li = 0.025))


    # test that works regardless of order of the input replicates
    knz_001dreorder <- knz_001d[order(knz_001d$abundance,
                                      knz_001d$year,
                                      knz_001d$species),]

    myresults_reorder <- variance_ratio(knz_001dreorder,
                                        time.var = "year",
                                        species.var = "species",
                                      abundance.var = "abundance",
                                      bootnumber = 1,
                                      replicate = "subplot")


    expect_equal(myresults[,4], myresults_reorder[,4])

    ##test that works regardless of order of the input replicates if average.replicates=F
    #test the class returned with default settings
    myresultsnoavg <- variance_ratio(knz_001d,
                                     time.var = "year",
                                     species.var = "species",
                                     abundance.var = "abundance",
                                     bootnumber = 1, replicate = "subplot",
                                     average.replicates = FALSE)

    expect_equal(class(myresultsnoavg), "data.frame")
    expect_equal(nrow(myresultsnoavg), 20)
    expect_equal(myresultsnoavg[1,5], 0.9694810, tolerance=0.00001)


    # test that works regardless of order of the input replicates
    knz_001dreorder <- knz_001d[order(knz_001d$abundance, knz_001d$year, knz_001d$species),]
    myresults_reordernoavg <- variance_ratio(knz_001dreorder, time.var="year",
                                             species.var="species",
                                             abundance.var="abundance",
                                             bootnumber=1, replicate="subplot",
                                             average.replicates = FALSE)
    expect_equal(myresultsnoavg[,5], myresults_reordernoavg[,5])
    expect_equal(myresultsnoavg[,1], myresults_reordernoavg[,1])


    #test that it also works with alternate column names
    knz_001d2<-knz_001d
    names(knz_001d2) = c("sp", "yr", "sub", "abund")
    myresults2<-variance_ratio(knz_001d2,"yr", "sp", "abund", 1, "sub")
    expect_that(sum(myresults2$VR), equals(sum(myresults$VR)))
    myresults2.2<-variance_ratio(knz_001d2, "yr", "sp",  "abund", 1, "sub", average.replicates = FALSE)
    expect_warning(variance_ratio(knz_001d2, "yr", "sp",  "abund", 1, NA, average.replicates = FALSE))
    expect_warning(variance_ratio(knz_001d2,"yr", "sp",  "abund", 1, NA, average.replicates = TRUE))

    #test that it works even if there are additional unused columns
    knz_001d3<-knz_001d
    knz_001d3$site<-"KNZ"
    myresults3<-variance_ratio(knz_001d3, "year", "species",  "abundance",1, "subplot")
    expect_that(myresults3$VR, is_identical_to(myresults$VR))

    ##test that works with replicate = NA
    myresults4 <- variance_ratio(dat1, time.var = "year",
                                 species.var = "species",
                                 abundance.var = "abundance",
                                 bootnumber = 10,
                                 replicate.var = NA)

    #test that works with one, specified replicate
    myresults5 <- variance_ratio(dat1, time.var = "year",
                                 species.var = "species",
                                 abundance.var = "abundance",
                                 bootnumber = 10,
                                 replicate.var = "subplot")

    #test that works with replicate as a character
    myresults6 <- variance_ratio(dat3, time.var = "year",
                                species.var = "species",
                                abundance.var = "abundance",
                                bootnumber = 10,
                                replicate.var = "subplot")

    kzcharacter <- knz_001d
    kzcharacter$subplot <- as.character(kzcharacter$subplot)
    myresults6.5 <- variance_ratio(kzcharacter, time.var = "year",
                                 species.var = "species",
                                 abundance.var = "abundance",
                                 bootnumber = 10,
                                 replicate.var = "subplot",
                                 average.replicates = F)

    #test that gives same for a single rep if average reps is false
    myresults7 <- variance_ratio(dat3, time.var = "year",
                                 species.var = "species",
                                 abundance.var = "abundance",
                                 bootnumber = 1,
                                 replicate.var = "subplot",
                                 average.replicates = FALSE)

    #test that all give the same VR value
    expect_that(myresults4$VR, is_identical_to(myresults5$VR))
    expect_that(myresults4$VR, is_identical_to(myresults6$VR))
    expect_that(myresults4$VR, is_identical_to(myresults7$VR))

    ##test variance_ratio_matrixdata
    calVRresults<-variance_ratio_matrixdata(datmat)
    expect_that(calVRresults, equals(myresults7$VR))
    expect_true(calVRresults>=0)
    expect_true(is.numeric(calVRresults))

    ##test variance_ratio_longformdata
    lfresults<-variance_ratio_longformdata(dat1,  time.var="year",species.var="species", abundance.var="abundance")
    expect_that(lfresults, equals(myresults7$VR))
    expect_true(lfresults>=0)
    expect_true(is.numeric(lfresults))

    ## test that bad data values are handled with graceful error messages
    bad_data <- knz_001d
    bad_data['abundance'][1,] <- 'dung'
    expect_error(variance_ratio(bad_data, "year", "species",
                                "abundance",  bootnumber=1,
                                replicate="subplot"),
                 "not numeric")
})
