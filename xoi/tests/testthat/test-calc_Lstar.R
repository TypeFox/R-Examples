
context("calc_Lstar")

test_that("calc_Lstar works", {

    # interference + large chr : Lstar == L
    expect_equal(calc_Lstar(500, 10, 0), 500)

    # p==1 is same as m==0
    expect_equal(calc_Lstar(100, 5, 1), calc_Lstar(100, 0, 0))

    # no interference case
    L <- seq(55, 205, by=25)
    tol <- (.Machine$double.eps)^(1/3)
    for(i in L) {
        Lstar <- calc_Lstar(i, 0, 0)
        expect_equal(i*2, 2*Lstar/(1-dpois(0, Lstar/50)), tolerance=tol)
    }

    # expected no. chiasmata under chi-square model
    expected_chiasmata <-
    function(Lstar, m=0, max_pts)
    {
        lambda <- Lstar/50*(m+1)
        if(missing(max_pts) || is.null(max_pts))
            max_pts <- qpois(1-1e-8, lambda)+10

        x <- 0:max_pts
        p <- dpois(x, lambda)
        xstar <- x %/% (m+1)
        xstar <- c(xstar, xstar+1)
        pstar <- (x %% (m+1))/(m+1)
        pstar <- c(p*(1-pstar), p*pstar)

        expected <- sum(pstar * xstar)
        expected/sum(pstar[xstar > 0])
    }

    # plain chi-square model (p=0)
    for(i in L) {
        Lstar <- calc_Lstar(i, 1, 0)
        expect_equal(i/50, expected_chiasmata(Lstar, 1), tolerance=tol)
    }

})

