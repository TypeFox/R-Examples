context("turnover")

test_that("turnover loads and returns correct result", {
    # Ensure that trivial tests work correctly
    expect_that(length("a"), equals(1))
    library(codyn)
    # Load our example data set
    data(knz_001d)
    expect_that(names(knz_001d)[4], equals("abundance"))
    
    #give new column names
    knz_001d2 <- knz_001d
    names(knz_001d2) = c("sp", "yr", "sub", "abund")
    
    #add a random character and factor column
    knz_001d2$randcharacter <- sample(letters, size = nrow(knz_001d2), replace = T)
    knz_001d2$randfactor <- as.factor(knz_001d2$randcharacter)
    
    #take a subset
    dat1 <- subset(knz_001d, knz_001d$subplot=="A_1")
    
    #rename the subset
    dat2 <- dat1
    names(dat2) = c("sp", "yr", "sub", "abund")
    
    #make subplot a character
    dat3 <- dat1
    dat3$subplot <- as.character(dat3$subplot)
    
    #take a two replicate subset
    dat4 <- subset(knz_001d, subplot == "A_1"|subplot == "A_2")
    
    # test that calculation from turnover is correct and does not regress
    myresults <- turnover(knz_001d, time.var="year", species.var="species", 
                        abundance.var="abundance", 
                        replicate.var="subplot", 
                        metric="total")
    expect_equal(class(myresults), "data.frame")
    expect_equal(nrow(myresults), 460)
    expect_equal(myresults[460,1], 0.47368421, tolerance=0.00001)
    expect_equal(sum(myresults[,1]), 116.2359, tolerance=0.00001)
    
    # test that works regardless of order of the input replicates
    knz_001dreorder <- knz_001d[order(knz_001d$abundance, knz_001d$year, knz_001d$species),]
    myresults_reorder <- turnover(knz_001dreorder, time.var="year", species.var="species", abundance.var="abundance", 
                                replicate.var="subplot", 
                                metric="total")
    expect_equal(myresults, myresults_reorder)
    
    #test that works regardless of whether parameter is specified or just ordered
    myresults2 <- turnover(df=knz_001d, replicate.var="subplot", time.var = "year",
                         species.var = "species", abundance.var = "abundance",
                         metric="total")    
    expect_identical(myresults, myresults2)
    
    #test that total is the default metric
    myresults3 <- turnover(df=knz_001d, time.var="year", species="species", abundance="abundance", replicate.var="subplot")
    expect_identical(myresults, myresults3)
    
    #test that works with different column nmaes
    myresults4 <- turnover(df=knz_001d2, time.var="yr", species="sp", abundance="abund", replicate="sub", metric="total")
    expect_equal(sum(myresults$total), sum(myresults4$total))
    
    #test that works with different column orders if names specified
    myresults5 <- turnover(df=knz_001d2, time.var="yr", abundance="abund", replicate="sub", metric="total", species="sp")
    expect_equal(sum(myresults$turnover), sum(myresults5$turnover))
    
    #test that gives an error if abundance is a character or factor
    expect_error(turnover(df=knz_001d2, time.var="yr", species.var="sp", abundance.var="randcharacter", replicate.var="sub"))
    expect_error(turnover(df=knz_001d2, time.var="yr", species.var="sp", abundance.var="randfactor", replicate.var="sub"))
    
    # test that stops if more than one record for a species within a year, in one replicate
    dat4 = rbind(dat1, dat1[nrow(dat1),])
    expect_error(turnover(dat4, time.var="year", species.var="species", abundance.var="abundance", replicate.var="subplot"))
    
    dat5 = rbind(knz_001d, knz_001d[nrow(knz_001d),], knz_001d[1,])
    expect_error(turnover(dat5, time.var="year", species.var="species", abundance.var="abundance", replicate.var="subplot"))
    
    #Last test: manually calculate turnover for a mini-data set
    test.spp1 <- data.frame(species=letters[1:5])
    test.spp2 <- data.frame(species=letters[2:7])
    # total turnover: 4 shared species. 1 disappears, 2 appear. Total richness is 7
    
    # default for turnover_twoyears is total; partial matching also being tested
    expect_equal(turnover_twoyears(test.spp1, test.spp2, species.var = "species"), (2 + 1) / 7)
    
    expect_equal(turnover_twoyears(test.spp1, test.spp2, species.var = "species", metric = "tot"), (2 + 1) / 7)
    
    expect_equal(turnover_twoyears(test.spp1, test.spp2, species.var = "species", metric = "dis"), 1 / 7)
    
    expect_equal(turnover_twoyears(test.spp1, test.spp2, species.var = "species", metric = "app"), 2 / 7)
    
    # adding a test to make sure that the output actually is the metric specified
    tot <- turnover(knz_001d, time.var="year", species.var="species", abundance.var="abundance", replicate.var="subplot", metric = "total")
    expect_equal(names(tot)[1], "total")
    
    app <- turnover(knz_001d, time.var="year", species.var="species", abundance.var="abundance", replicate.var="subplot", metric = "appear")
    expect_equal(names(app)[1], "appearance")
    
    dis <- turnover(knz_001d, time.var="year", species.var="species", abundance.var="abundance", replicate.var="subplot", metric = "disap")
    expect_equal(names(dis)[1], "disappearance")
    
    # same, now for a single replicate
    
    tot <- turnover(dat1, time.var="year", species.var="species", abundance.var="abundance", metric = "total")
    expect_equal(names(tot)[1], "total")
    
    app <- turnover(dat1, time.var="year", species.var="species", abundance.var="abundance", replicate.var="subplot", metric = "appear")
    expect_equal(names(app)[1], "appearance")
    
    dis <- turnover(dat1, time.var="year", species.var="species", abundance.var="abundance", replicate.var="subplot", metric = "disap")
    expect_equal(names(dis)[1], "disappearance")
    
    # checking to make sure that when using a single replicate, should not be a problem if subplot is still specified
    tot.sub <- turnover(dat1, time.var="year", species.var="species", abundance.var="abundance", metric = "total", replicate.var = "subplot")
    expect_equal(tot[1:2], tot.sub[1:2])
    
    })

