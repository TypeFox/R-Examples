## ---- results = "hide", echo=FALSE, message = FALSE----------------------
library(matchingR)

## ------------------------------------------------------------------------
# number of men and women, respectively
n = 100
# level of commonality
lambda = 0.5
# men's preferences
uM = lambda * matrix(runif(n), nrow = n, ncol = n) +
    (1 - lambda) * runif(n ^ 2)

## ------------------------------------------------------------------------
# womens's preferences
uW = lambda * matrix(runif(n), nrow = n, ncol = n) +
    (1 - lambda) * runif(n ^ 2)

## ------------------------------------------------------------------------
# market sizes
N = c(100, 500, 1000)

## ------------------------------------------------------------------------
# levels of commonality
Commonality = c(0.25, 0.5, 0.75)

## ------------------------------------------------------------------------
set.seed(1)

test_galeShapley.marriageMarket = function(n, lambda) {
    uM = lambda * matrix(runif(n), nrow = n, ncol = n) +
        (1 - lambda) * runif(n ^ 2)
    uW = lambda * matrix(runif(n), nrow = n, ncol = n) +
        (1 - lambda) * runif(n ^ 2)
    galeShapley.marriageMarket(uM, uW)
}

## ---- echo=FALSE---------------------------------------------------------
res1 = microbenchmark::microbenchmark(
    test_galeShapley.marriageMarket(N[1], Commonality[1]),
    test_galeShapley.marriageMarket(N[2], Commonality[1]),
    test_galeShapley.marriageMarket(N[3], Commonality[1]), times = 1, unit = "s")

res2 = microbenchmark::microbenchmark(
    test_galeShapley.marriageMarket(N[1], Commonality[2]),
    test_galeShapley.marriageMarket(N[2], Commonality[2]),
    test_galeShapley.marriageMarket(N[3], Commonality[2]), times = 1, unit = "s")

res3 = microbenchmark::microbenchmark(
    test_galeShapley.marriageMarket(N[1], Commonality[3]),
    test_galeShapley.marriageMarket(N[2], Commonality[3]),
    test_galeShapley.marriageMarket(N[3], Commonality[3]), times = 1, unit = "s")

## ---- echo=FALSE---------------------------------------------------------
table = matrix(c(summary(res1, unit = 's')$mean,
                 summary(res2, unit = 's')$mean,
                 summary(res3, unit = 's')$mean), ncol = length(Commonality), 
               dimnames = list("N" = gsub(" ", "", sprintf("N=%s", format(N, big.mark = ","))), 
                               "Commonality" = sprintf("lambda=%.2f" , Commonality)))
knitr::kable(table, digits = 6, caption = "Run time (in seconds) for the matching of men to women")

## ------------------------------------------------------------------------
u = matrix(rnorm(16), nrow = 4)
results = roommate(utils = u)

## ------------------------------------------------------------------------
test_roommate = function(n) {
    roommate(utils = matrix(rnorm( 4 * n ^ 2 ), nrow = 2 * n))
}

## ---- echo=FALSE---------------------------------------------------------
N = c(2, 3, 4, 5, 6, 7, 8, 9, 10, 
      11, 12, 13, 14, 15, 16, 17, 
      18, 19, 20)

res_irving = microbenchmark::microbenchmark(
    test_roommate(N[1]),
    test_roommate(N[2]),
    test_roommate(N[3]),
    test_roommate(N[4]),
    test_roommate(N[5]),
    test_roommate(N[6]),
    test_roommate(N[7]),
    test_roommate(N[8]),
    test_roommate(N[9]),
    test_roommate(N[10]),
    test_roommate(N[11]),
    test_roommate(N[12]),
    test_roommate(N[13]),
    test_roommate(N[14]),
    test_roommate(N[15]),
    test_roommate(N[16]),
    test_roommate(N[17]),
    test_roommate(N[18]),
    test_roommate(N[19]), times = 1, unit = "s")

solution_exists = rep(0, 19);
for (n in 1:19) {
    for (i in 1:50) {
        results = roommate(utils = matrix(rnorm( 4 * N[n] ^ 2 ), nrow = 2 * N[n]))
        if (!is.null(results)) {
            solution_exists[n] = solution_exists[n] + 1;
        }
    }
}

irving_original_prop = c(0.961, 0.931, 0.909, 0.868, 
                         0.871, 0.867, 0.857, 0.830, 
                         0.815, 0.808, 0.816, 0.788, 
                         0.750, 0.766, 0.755, 0.770, 
                         0.740, 0.725, 0.745, 0.775, 
                         0.740, 0.740, 0.730, 0.710, 
                         0.725, 0.670, 0.675, 0.690)

## ---- echo=FALSE---------------------------------------------------------
table = matrix(c(rep(50, 19), solution_exists, solution_exists/50, irving_original_prop[1:19], summary(res_irving)$mean),  
               ncol = 5, nrow = 19, 
               dimnames = list("N" = gsub(" ", "", sprintf("N = %s", format(N, big.mark = ","))), 
                               "t" = c("No. of instances", 
                                       "No. with solution", 
                                       "Proportion with solution", 
                                       "Proportion with solution (Irving 1985)",
                                       "Average cpu time")))
knitr::kable(table, digits = 6, caption = "Proportion of matchings which exist and run time (in seconds) for the stable roommate algorithm.")

## ------------------------------------------------------------------------
u = matrix(rnorm(16), nrow = 4)
results = toptrading(utils = u)

## ------------------------------------------------------------------------
test_toptrading = function(n) {
    toptrading(utils = matrix(rnorm( n ^ 2 ), nrow = n))
}

## ---- echo=FALSE---------------------------------------------------------
N = c(2, 4, 8, 16, 32, 64, 128, 256, 512, 1024)

res_irving = microbenchmark::microbenchmark(
    test_toptrading(N[1]),
    test_toptrading(N[2]),
    test_toptrading(N[3]),
    test_toptrading(N[4]),
    test_toptrading(N[5]),
    test_toptrading(N[6]),
    test_toptrading(N[7]),
    test_toptrading(N[8]),
    test_toptrading(N[9]),
    test_toptrading(N[10]), times = 1, unit = "s")

## ---- echo=FALSE---------------------------------------------------------
table = matrix(c(summary(res_irving)$mean),  
               ncol = 1, nrow = 10, 
               dimnames = list("N" = gsub(" ", "", sprintf("N = %s", format(N, big.mark = ","))), 
                               "t" = c("Average cpu time")))
knitr::kable(table, digits = 6, caption = "Run time (in s) for top trading cycle algorithm")

