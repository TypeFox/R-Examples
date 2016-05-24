library(nestedRanksTest)
data(woodpecker_multiyear)

context("nestedRanksTest, identity of results with each interface")

it <- 100

set.seed(42)
a.resf1 <- nestedRanksTest(Distance ~ Year | Granary, 
                         data = woodpecker_multiyear,
                         subset = Species == "agrifolia",
                         n.iter = it)

set.seed(42)
a.resf2 <- nestedRanksTest(Distance ~ Year | Granary, 
                         data = subset(woodpecker_multiyear,
                                       Species == "agrifolia"),
                         n.iter = it)

set.seed(42)
l.resf1 <- nestedRanksTest(Distance ~ Year | Granary, 
                         data = woodpecker_multiyear,
                         subset = Species == "lobata",
                         n.iter = it)

set.seed(42)
l.resf2 <- nestedRanksTest(Distance ~ Year | Granary, 
                         data = subset(woodpecker_multiyear,
                                       Species == "lobata"),
                         n.iter = it)

test_that("results identical regardless of subsetting", {
    expect_equal(a.resf1, a.resf2)
    expect_equal(l.resf1, l.resf2)
})

set.seed(42)
a.resf3 <- nestedRanksTest(Distance ~ Year | Granary, 
                         data = woodpecker_multiyear,
                         subset = Species == "agrifolia",
                         lightweight = TRUE,
                         n.iter = it)

a.resf4 <- a.resf2
a.resf4$null.distribution <- NULL

set.seed(42)
l.resf3 <- nestedRanksTest(Distance ~ Year | Granary, 
                         data = woodpecker_multiyear,
                         subset = Species == "lobata",
                         lightweight = TRUE,
                         n.iter = it)

l.resf4 <- l.resf2
l.resf4$null.distribution <- NULL


test_that("results identical with and without lightweight=", {
    expect_equal(a.resf3, a.resf4)
    expect_equal(l.resf3, l.resf4)
})

set.seed(42)
a.resf5 <- with(subset(woodpecker_multiyear, Species == "agrifolia"),
              nestedRanksTest(y = Distance, x = Year, groups = Granary,
                              n.iter = it))
set.seed(42)
a.resf6 <- with(subset(woodpecker_multiyear, Species == "agrifolia"),
              nestedRanksTest(Year, Distance, Granary, n.iter = it))
set.seed(42)
l.resf5 <- with(subset(woodpecker_multiyear, Species == "lobata"),
              nestedRanksTest(y = Distance, x = Year, groups = Granary,
                              n.iter = it))
set.seed(42)
l.resf6 <- with(subset(woodpecker_multiyear, Species == "lobata"),
              nestedRanksTest(Year, Distance, Granary, n.iter = it))
set.seed(42)
ldat <- subset(woodpecker_multiyear, Species == "lobata")
l.resf7 <- nestedRanksTest(ldat$Year, ldat$Distance, ldat$Granary, n.iter = it)
set.seed(42)
l.resf8 <- nestedRanksTest(Distance ~ Year, group = Granary,
                           data = woodpecker_multiyear,
                           subset = Species == "lobata",
                           n.iter = it)

test_that("results identical between formula and default interfaces", {
    expect_equal(a.resf1, a.resf5)
    expect_equal(a.resf1, a.resf6)
    expect_equal(l.resf1, l.resf5)
    expect_equal(l.resf1, l.resf6)
    expect_identical(l.resf6$statistic, l.resf7$statistic)
    expect_identical(l.resf6$p.value, l.resf7$p.value)
    expect_equal(l.resf1, l.resf8)
})

test_that("object identities correct regardless of interface", {
    expect_is(a.resf1, "htest_boot")
    expect_is(a.resf1, "htest")
    expect_is(a.resf5, "htest_boot")
    expect_is(a.resf5, "htest")
    expect_is(a.resf6, "htest_boot")
    expect_is(a.resf6, "htest")
})

