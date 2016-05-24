library(nestedRanksTest)
data(woodpecker_multiyear)

context("nestedRanksTest, errors when using each interface")

adat = subset(woodpecker_multiyear, Species == "agrifolia")

test_that("Errors produced when specifying groups", {
    expect_error(nestedRanksTest(Distance ~ Year, data = adat),
                "invalid group specification in formula")
    expect_error(nestedRanksTest(Distance ~ Year | G, data = adat),
                "object 'G' not found")
    expect_error(nestedRanksTest(adat$Year, adat$Distance),
                "'groups' missing")
    expect_error(nestedRanksTest(Distance ~ Year | Granary, 
                                 groups = Granary, data = adat),
                "groups are specified with '|' in formula or with groups= argument, but not both")
    adat$Granary2 <- adat$Granary
    expect_error(nestedRanksTest(Distance ~ Year | Granary, 
                                 groups = Granary2, data = adat),
                "groups are specified with '|' in formula or with groups= argument, but not both")
})

test_that("Errors produced when other variable missing", {
    expect_error(nestedRanksTest(Year, Distance, data = adat),
                "object 'Year' not found")
})

test_that("Errors produced when wrong number of levels in treatment", {
    expect_error(nestedRanksTest(Distance ~ Year | Granary, 
                                 data = woodpecker_multiyear),
                 "must have exactly 2 levels")
    expect_error(nestedRanksTest(Distance ~ Year | Granary, 
                                 data = adat,
                                 subset = Year == "2007"),
                 "must have exactly 2 levels")
})

