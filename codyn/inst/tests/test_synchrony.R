context("synchrony")

test_that("synchrony loads and returns correct result", {

    # Load our example data set
    data(knz_001d)
    expect_that(names(knz_001d)[4], equals("abundance"))

    # From community_stability

    # make new data frame with different column names
    knz_001d2 <- knz_001d
    names(knz_001d2)=c("sp", "yr", "sub", "abund")

    #add a random character and factor column
    knz_001d2$randcharacter<-sample(letters, size = nrow(knz_001d2), replace = T)
    knz_001d2$randfactor<-as.factor(knz_001d2$randcharacter)

    #take a subset
    dat1 <- subset(knz_001d, knz_001d$subplot=="A_1")

    #rename the subset
    dat2<-dat1
    names(dat2)=c("sp", "yr", "sub", "abund")

    #make subplot a character
    dat3<-dat1
    dat3$subplot<-as.character(dat3$subplot)

    #test the synch_onerep function
    myresults<-synch_onerep(dat1, time.var = "year", abundance.var="abundance", species.var="species")
    expect_that(class(myresults), equals("numeric"))
    expect_that(length(myresults), equals(1))
    expect_that(myresults, equals(0.114323, tolerance=.0001))

    #test that works with different column names
    myresults2<-synch_onerep(dat2, time.var = "yr", abundance.var="abund", species.var="sp")
    expect_that(myresults, equals(myresults2))

    #test that gives a warning if running on factor instead of numeric
    expect_error(synchrony(dat2, time.var = "yr", abundance.var="sub", species.var="sp"), "not numeric")

    #test that works on a single replicate
    myresults3 <- synchrony(dat1, time.var="year", species.var="species", abundance.var="abundance", replicate.var=NA)
    expect_that(myresults3, equals(myresults2))

    # test that stops if more than one record for a species within a year, in one replicate
    dat4 = rbind(dat1, dat1[nrow(dat1),])
    expect_warning(synchrony(dat4, time.var="year", species.var="species", abundance.var="abundance", replicate.var=NA))

    dat5 = rbind(knz_001d, knz_001d[nrow(knz_001d),], knz_001d[1,])

    expect_error(synchrony(dat5, replicate.var="subplot"))

    #test that will still run if there are missing levels in a factor "replicate"; deleting levels that are NaN
    myresults4<-synchrony(dat1, time.var="year", species.var = "species", abundance.var="abundance", replicate.var= "subplot")

    #this will give a warning because replicate is a factor without all values present in dat1 - the warning is a good thing
    myresults5<-as.numeric(myresults4[2])
    expect_that(myresults5, equals(myresults3))

    #test that works whether replicate is a character or factor
    myresults6<-synchrony(dat3, replicate.var="subplot", species.var = "species", time.var="year", abundance.var="abundance")
    expect_that((myresults6[1,2]), equals(myresults3))

    #test that works with multiple replicates
    myresults7<-synchrony(knz_001d, replicate.var="subplot", time.var="year", abundance.var="abundance", species.var = "species")
    expect_that(myresults6[1,2], equals(myresults7[1,2]))

    # test that works regardless of order of the input replicates
    knz_001dreorder <-knz_001d[order(knz_001d$abundance, knz_001d$year, knz_001d$species),]
    myresults_reorder<-synchrony(knz_001dreorder, replicate.var="subplot", time.var="year", abundance.var="abundance", species.var = "species")
    expect_equal(myresults7, myresults_reorder)


    #test that works with different column names
    myresults8<-synchrony(knz_001d2, replicate.var="sub", time.var="yr", abundance.var="abund", species.var = "sp")

    expect_that(myresults7[1,2], equals(myresults8[1,2]))

    # test that if species name is not specified, will cause error
    expect_error(synchrony(knz_001d2, replicate.var="sub", time.var="yr", abundance.var="abund"))

    #test that gives error when abundance column is a character or factor
    expect_error(synchrony(knz_001d2, replicate.var="sub", time.var="yr", abundance.var="randcharacter"))
    expect_error(synchrony(knz_001d2, replicate.var="sub", time.var="yr", abundance.var="randfactor"))

    # For the "Gross" metric, if a species doesn't vary at all over the time series it omits that species from focal species list and gives warning.
    # E.g, in KNZ_001d subplot D_5, salvia azurea has an abundance of 3 across the time series.
    dat5 <- knz_001d[knz_001d$subplot == "D_5",]
    expect_warning(synchrony(dat5, replicate.var = NA, species.var = "species", time.var = "year", 
                             abundance.var = "abundance", metric = "G"))

    # Checks Loreau method with the same subplot, for a single plot vs whole dataset
    myres7 <- synchrony(dat5, replicate.var = NA, species.var = "species", time.var = "year", 
                        abundance.var = "abundance",
                        metric = "L")
    myres8 <- synchrony(knz_001d, replicate.var = "subplot", 
                        species.var = "species", time.var = "year", 
                        abundance.var = "abundance", metric = "L")
    expect_equal(myres7, myres8[myres8$subplot=="D_5","synchrony"])

    #test that Gross works with multiple replicates

    expect_warning(synchrony(knz_001d, species.var = "species", replicate.var="subplot", time.var="year", abundance.var="abundance", metric="Gross"))

    myresults7gross <- suppressWarnings(synchrony(knz_001d, species.var = "species", replicate.var="subplot", time.var="year", abundance.var="abundance", metric="Gross"))

    # test that works regardless of order of the input replicates
    knz_001dreorder <- knz_001d[order(knz_001d$abundance, knz_001d$year, knz_001d$species),]
    myresults_reordergross <- suppressWarnings(synchrony(knz_001dreorder, species.var = "species", replicate.var = "subplot", time.var = "year", abundance.var = "abundance", metric="Gross"))

    expect_equal(myresults7gross, myresults_reordergross)

})