# test whether identical results are returned by different interfaces

library(nestedRanksTest)
data(woodpecker_multiyear)

run.tests <- function(test.species, it = 100) {
    set.seed(42)
    a <- nestedRanksTest(Distance ~ Year | Granary, 
                    data = woodpecker_multiyear,
                    subset = Species == test.species,
                    n.iter = it)
    print(a)
    b <- nestedRanksTest(Distance ~ Year | Granary, 
                    data = subset(woodpecker_multiyear,
                                  Species == test.species),
                    n.iter = it)
    print(b)
    c <- nestedRanksTest(Distance ~ Year | Granary, 
                    data = woodpecker_multiyear,
                    subset = Species == test.species,
                    lightweight = TRUE,
                    n.iter = it)
    print(c)
    d <- with(subset(woodpecker_multiyear, Species == test.species),
         nestedRanksTest(y = Distance, x = Year, groups = Granary,
                         n.iter = it))
    print(d)
    e <- with(subset(woodpecker_multiyear, Species == test.species),
         nestedRanksTest(Year, Distance, Granary, n.iter = it))
    print(e)
}

# run.tests("agrifolia", 100)

run.tests("lobata", 100)
