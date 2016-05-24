library(ECctmc)

context("Simulating endpoint conditioned CTMC paths")

# This function merely executes a basic t test at a set of discrete times to
# check that the state at those times is not different when simulating using
# uniformization and modified rejection sampling. A full test comparing these
# algorithms to sample paths produced by the GillespieSSA package (endpoint
# conditioned simple rejection sampling) is available as an commented extended R
# script in the /tests directory of this package.
test_that("Modified rejection sampling and uniformization target the same distribution", {

        # set the seed, observation times, and parameters
        set.seed(52787)
        obstimes <- 0:5
        niter <- 500
        params <- c(mu = rnorm(1, 0.5, 1e-3), beta = rnorm(1, 0.5, 1e-3))

        # construct the rate matrix
        Q <- matrix(c(-params[1], params[1], params[2], -params[2]), nrow = 2, byrow = T)

        # function to determine the state at a sequence of times
        state_seq <- function(path, times) {
                path[findInterval(times, path[,1]), 2]
        }

        # sample the paths
        MR_paths <- ECctmc::sample_path(a=1, b=2, t0=0, t1=5, Q=Q, method = "mr", npaths=niter)
        Unif_paths <- ECctmc::sample_path(a=1, b=2, t0=0, t1=5, Q=Q, method = "unif", npaths=niter)

        # get the state at each observation time
        MR_res <- vapply(MR_paths, state_seq, FUN.VALUE = as.numeric(obstimes), times = obstimes)
        Unif_res <- vapply(Unif_paths, state_seq, FUN.VALUE = as.numeric(obstimes), times = obstimes)

        expect_true(all(MR_res[1,] == 1) & all(Unif_res[1,] == 1) & all(MR_res[6,]==2 & all(Unif_res[6,]==2)))

        for(k in 2:5) {
                expect_false(t.test(MR_res[k,], Unif_res[k,])$p.value < 0.1)
        }
})