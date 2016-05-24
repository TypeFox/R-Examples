context("rate_change")

test_that("rate_change loads and returns correct result", {
    # Ensure that trivial tests work correctly
    expect_equal(length("a"), 1)
    
    # Load our example data set
    data(knz_001d)
    expect_equal(names(knz_001d)[4], "abundance")
    
    #give new column names
    knz_001d2 <- knz_001d
    names(knz_001d2)=c("sp", "yr", "sub", "abund")
    
    #add a random character and factor column
    knz_001d2$randcharacter<-"rchar"
    knz_001d2$randfactor<-as.factor(knz_001d2$randcharacter)
    
    #take a subset
    dat1 <- subset(knz_001d, knz_001d$subplot=="A_1")
    
    #rename the subset
    dat2<-dat1
    names(dat2)=c("sp", "yr", "sub", "abund")
    
    #make subplot a character
    dat3<-dat1
    dat3$subplot<-as.character(dat3$subplot)
    
    #test the lagged_slope function
    myresults <- lagged_slope(dat1, "year", "species", "abundance")
    expect_equal(class(myresults), "numeric")
    expect_equal(length(myresults), 1)
    expect_equal(myresults, 0.6706902, tolerance=0.00001)
    #test that works with different column names
    myresults2<-lagged_slope(dat2,  "yr", "sp", "abund")
    expect_equal(myresults, myresults2)
    
    #test that gives a warning if running on factor instead of numeric
    expect_error(lagged_slope(dat2, "yr", "sp", "subplot"))
    
    #test the rate_change function
    #test that works on a single replicate
    myresults3<-rate_change(dat1, replicate.var=NA, time.var="year", species.var="species", abundance.var="abundance")
    expect_equal(myresults3, myresults2)
    
    # Test getting the lagged distance data as a data frame, by replicates
    dist1 <- rate_change_interval(knz_001d, replicate.var="subplot", time.var="year", species.var="species", abundance.var="abundance")
    expect_equal(class(dist1), "data.frame")
    expect_equal(ncol(dist1), 3)
    expect_equal(nrow(dist1), 5520)
    
    #test that will still run if there are missing levels in a factor "replicate"; deleting levels that are NaN
    myresults4<-rate_change(dat1, replicate.var="subplot", time.var="year", species.var="species", abundance.var="abundance")
    #this will give a warning because replicate is a factor without all values present in dat1 - the warning is a good thing
    myresults5<-as.numeric(myresults4[2])
    expect_equal(myresults5, myresults3)
    
    #test that works whether replicate is a character or factor
    myresults6<-rate_change(dat3, replicate.var="subplot", time.var="year", species.var="species", abundance.var="abundance")
    expect_equal((myresults6[1,2]), myresults3)
    
    #test that works with multiple replicates
    myresults7<-rate_change(knz_001d, replicate.var="subplot", time.var="year", species.var="species", abundance.var="abundance")
    expect_equal(myresults6[1,2], myresults7[1,2])
    
    #test that works with different column names
    myresults8<-rate_change(knz_001d2, replicate.var="sub", time.var="yr", species.var="sp", abundance.var="abund")
    expect_equal(myresults7[1,2], myresults8[1,2])
    
    #test that works regardless of whether parameter is specified or just ordered
    myresults8<-rate_change(knz_001d, replicate.var="subplot", time.var = "year", species.var = "species", abundance.var = "abundance")
    expect_identical(myresults8, myresults7)
    
    #test that works with different column orders if names specified
    myresults9<-rate_change(knz_001d, abundance.var="abundance", replicate.var="subplot", species.var="species", time.var="year")
    expect_identical(myresults9, myresults7)
    
    #test that works with different column names
    myresults10<-rate_change(knz_001d2, replicate.var="sub", time.var="yr", species.var="sp", abundance.var="abund")
    expect_equal(myresults10[1,2], myresults7[1,2])
    
    #test that it works even if there are additional unused columns
    knz_001d3<-knz_001d
    knz_001d3$site<-"KNZ"
    myresults11<-rate_change(knz_001d3, replicate.var = "subplot", time.var = "year", species.var = "species", abundance.var = "abundance")
    expect_identical(myresults11, myresults7)
    
    #test that gives error when abundance column is a character or factor
    expect_error(rate_change(knz_001d2, replicate.var="sub", time.var="yr", species.var="sp", abundance.var="randcharacter"))
    expect_error(rate_change(knz_001d2, replicate.var="sub", time.var="yr", species.var="sp", abundance.var="randfactor"))
    
    
    #test that works with multiple replicates regardless of the order
    myresults7<-rate_change(knz_001d, replicate.var="subplot", time.var="year", species.var="species", abundance.var="abundance")
    
    knz_001dreorder <-knz_001d[order(knz_001d$abundance, knz_001d$year, knz_001d$species),]
    myresults7reorder<-rate_change(knz_001dreorder, replicate.var="subplot", time.var="year", species.var="species", abundance.var="abundance")
    
    expect_equal(myresults7[,1], myresults7reorder[,1])
    expect_equal(myresults7[,2], myresults7reorder[,2])
    
    
    
    #test that rate_change_interval works with multiple replicates regardless of the order
    myresults7interval<-rate_change_interval(knz_001d, replicate.var="subplot", time.var="year", species.var="species", abundance.var="abundance")
    
    knz_001dreorder <-knz_001d[order(knz_001d$abundance, knz_001d$year, knz_001d$species),]
    myresults7reorderinterval<-rate_change_interval(knz_001dreorder, replicate.var="subplot", time.var="year", species.var="species", abundance.var="abundance")
    
    expect_equal(myresults7interval[,1], myresults7reorderinterval[,1])
    expect_equal(myresults7interval[,2], myresults7reorderinterval[,2])
    expect_equal(myresults7interval[,3], myresults7reorderinterval[,3])
    
})