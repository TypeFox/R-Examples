context('multivar')

## Generate test data without littering the environment with temporary
## variables
x <- NULL
y <- NULL
d <- NULL
o <- NULL
local({
    set.seed(123)
    f <- -1:1
    N <- 3
    T <- 2
    K <- 2
    rho <- .5
    beta <- c(.1, -.5)

    ## $x_i = 0.75 f + N(0, 1)$:
    x <- aperm(array(.75*f, dim=c(N, K, T+1)) + rnorm(N*K*(T+1)), 3:1)

    ## $y_{i,t} = \rho y_{i,t-1} + \beta x_{i,t} + f_i + N(0,1)$:
    y <- matrix(0, T+1, N)
    for (t in 1:(T+1)) {
        yy <- if (t>1) y[t-1, ] else 0
        y[t,] <- rho * yy + f  + beta %*% x[t,,] + rnorm(N)
    }

    x <<- x
    y <<- y
    d <<- data.frame(i=rep(1:3, each=3), t=rep(1:3, 3),
                     matrix(aperm(x, c(1, 3, 2)), N*(T+1), K,
                                    dimnames = list(NULL, paste0('x', seq(K)))),
                     y = c(y))

    set.seed(123)
    o <<- opm(y~x1+x2, d, n.samp = 10)
})


## Sanity check
test_that('data', {
    expect_equal(x, array(c(-1.31047564655221, -0.289083794010798,
                            -0.349228549405948, -0.679491608575424,
                            -1.19566197009996, 1.03691313680308,
                            -0.23017748948328, -1.26506123460653,
                            0.11068271594512, 0.129287735160946,
                            1.22408179743946, 0.497850478229239,
                            2.30870831414912, 0.0631471481064739,
                            0.194158865245925, 2.46506498688328,
                            1.10981382705736, -1.21661715662964),
                           dim = c(3, 2, 3)))
    expect_equal(y, matrix(c(-0.0899458588038237, -0.694025238411308,
                             -2.52543131039704, -0.560453024256735,
                             -2.04477798261599, -2.94693926957052,
                             -1.06948536801357, -0.812226112015961,
                             2.05939845332596),
                           nrow = 3, ncol = 3))
})


test_that('default', {
    set.seed(123)
    expect_equal(opm(x, y, n = 10),
                 structure(list(samples = list(rho = c(0.805, 0.425,
                                                       -0.181, 0.763,
                                                       0.581, 0.965,
                                                       0.0559999999999999, 0.533,
                                                       0.648, -0.086),
                                               sig2 = 1 / c(11.6175662271191,
                                                            0.919115793520271,
                                                            0.989369792715488,
                                                            0.00440918521157215,
                                                            1.35621264272223,
                                                            4.24440799476086,
                                                            0.0931249250470442,
                                                            1.56847551750474,
                                                            0.371922970545138,
                                                            5.00022577787265),
                                               beta = matrix(c(0.149966256134789, -1.07201170366086,
                                                               -3.32784640301039, 6.31976053552415,
                                                               0.399370400131068, -0.409833517459282,
                                                               0.942409502841064, -1.01613441888489,
                                                               -1.06148590924992, -1.70364397749032,
                                                               -0.746392848950218, -1.22143920569537,
                                                               -1.63866906134468, 7.01947247537598,
                                                               -0.536929048279301, -0.845295844019301,
                                                               -0.565341215754439, 0.0765501718682319,
                                                               -0.171463431207414, -1.13494915857112),
                                                             10, 2)),
                                call = quote(opm(x = x, y = y, n.samp = 10)),
                                time.indicators = FALSE),
                           class = 'opm'))
})


test_that('named matrix', {
    set.seed(123)
    colnames(x) <- c('x1', 'x2')
    expect_equal(opm(x, y, n = 10),
                 structure(list(samples = list(rho = c(0.805, 0.425,
                                                       -0.181, 0.763,
                                                       0.581, 0.965,
                                                       0.0559999999999999, 0.533,
                                                       0.648, -0.086),
                                               sig2 = 1 / c(11.6175662271191,
                                                            0.919115793520271,
                                                            0.989369792715488,
                                                            0.00440918521157215,
                                                            1.35621264272223,
                                                            4.24440799476086,
                                                            0.0931249250470442,
                                                            1.56847551750474,
                                                            0.371922970545138,
                                                            5.00022577787265),
                                               beta = matrix(c(0.149966256134789, -1.07201170366086,
                                                               -3.32784640301039, 6.31976053552415,
                                                               0.399370400131068, -0.409833517459282,
                                                               0.942409502841064, -1.01613441888489,
                                                               -1.06148590924992, -1.70364397749032,
                                                               -0.746392848950218, -1.22143920569537,
                                                               -1.63866906134468, 7.01947247537598,
                                                               -0.536929048279301, -0.845295844019301,
                                                               -0.565341215754439, 0.0765501718682319,
                                                               -0.171463431207414, -1.13494915857112),
                                                             10, 2, dimnames = list(NULL, c('x1', 'x2')))),
                                call = quote(opm(x = x, y = y, n.samp = 10)),
                                time.indicators = FALSE),
                           class = 'opm'))
})


test_that('formula', {
    expected <- structure(list(samples = list(rho = c(0.805, 0.425,
                                                      -0.181, 0.763,
                                                      0.581, 0.965,
                                                      0.0559999999999999, 0.533,
                                                      0.648, -0.086),
                                              sig2 = 1 / c(11.6175662271191,
                                                           0.919115793520271,
                                                           0.989369792715488,
                                                           0.00440918521157215,
                                                           1.35621264272223,
                                                           4.24440799476086,
                                                           0.0931249250470442,
                                                           1.56847551750474,
                                                           0.371922970545138,
                                                           5.00022577787265),
                                              beta = matrix(c(0.149966256134789, -1.07201170366086,
                                                              -3.32784640301039, 6.31976053552415,
                                                              0.399370400131068, -0.409833517459282,
                                                              0.942409502841064, -1.01613441888489,
                                                              -1.06148590924992, -1.70364397749032,
                                                              -0.746392848950218, -1.22143920569537,
                                                              -1.63866906134468, 7.01947247537598,
                                                              -0.536929048279301, -0.845295844019301,
                                                              -0.565341215754439, 0.0765501718682319,
                                                              -0.171463431207414, -1.13494915857112),
                                                            10, 2, dimnames = list(NULL, c('x1', 'x2')))),
                               call = quote(opm(x = y~x1+x2, data = d, n.samp = 10)),
                               time.indicators = FALSE,
                               index = c('i', 't'),
                               terms = structure(y~x1+x2,
                                                 variables = quote(list(y, x1, x2)),
                                                 factors = matrix(c(0, 1, 0, 0, 0, 1), 3, 2,
                                                                  dimnames=list(c('y', 'x1', 'x2'),
                                                                                c('x1', 'x2'))),
                                                 term.labels = c('x1', 'x2'),
                                                 order = c(1L, 1L),
                                                 intercept = 0L,
                                                 response = 1L,
                                                 .Environment = parent.frame(),
                                                 predvars = quote(list(y, x1, x2)),
                                                 dataClasses = c(y='numeric', x1='numeric', x2='numeric'),
                                                 class=c('terms', 'formula'))),
                          class = 'opm')

    ## default index is the first two columns of the data
    set.seed(123)
    expect_equal(opm(y~x1+x2, d,
                     n.samp = 10),
                 expected)

    d <- d[,c('y', 'x1', 't', 'i', 'x2')]

    ## specify index column by position
    set.seed(123)
    expected$call <- quote(opm(x = y~x1+x2, data = d,
                               index = 4:3, n.samp = 10))
    expect_equal(opm(y~x1+x2, d, index=4:3,
                     n.samp = 10),
                 expected)
   
    ## specify index column by name
    set.seed(123)
    expected$call <- quote(opm(x = y~x1+x2, data = d,
                               index = c('i', 't'), n.samp = 10))
    expect_equal(opm(y~x1+x2, d, index=c('i', 't'), n.samp = 10),
                 expected)
   
})


test_that('confint', {
    expect_equal(confint(o),
                 matrix(c(-0.159625, 0.111707292991574,
                          -2.96240085726837, -1.54479234382358,
                          0.929, 178.18554581309,
                          5.10985655317045, 5.45731495708674),
                        nrow = 4, ncol = 2,
                        dimnames = list(c("rho", "sig2", "beta.x1", "beta.x2"),
                                        c("2.5%", "97.5%"))))

    expect_equal(confint(o, 1:2, level = 0.68),
                 matrix(c(-0.0235200000000001, 0.215660754512765,
                          0.78652, 7.19646832807839),
                        nrow = 2, ncol = 2,
                        dimnames = list(c("rho", "sig2"), c("16%", "84%"))))
})


test_that('coef', {
    expect_equal(coef(o),
                 c(rho = 0.557,
                   sig2 = 0.874045960780248,
                   beta.x1 = -0.712983968172086,
                   beta.x2 = -0.655867032352328))

    expect_error(coef(o, probs = .6), "Arguments 'probs' and 'names' are not allowed")
    expect_error(coef(o, names = TRUE), "Arguments 'probs' and 'names' are not allowed")
})


test_that('print', {
    expect_equal(capture.output(o),
                 c('Call:',
                   'opm(x = y ~ x1 + x2, data = d, n.samp = 10)',
                   '',
                   'Coefficients:',
                   '                 mean (SD)       med           95-CI',
                   'rho        0.450900 (0.39)   0.55700   (-0.16, 0.93)',
                   'sig2     24.422159 (71.18)   0.87405  (0.11, 178.19)',
                   'beta.x1   -0.077945 (2.55)  -0.71298     (-3.0, 5.1)',
                   'beta.x2    0.023554 (2.51)  -0.65587     (-1.5, 5.5)',
                   ''))
})


test_that('summary', {
    expect_equal(capture.output(summary(o)),
                 c('Call:',
                   'opm(x = y ~ x1 + x2, data = d, n.samp = 10)',
                   '',
                   'Parameter estimates:',
                   '           <--95CI    <--68CI        med     68CI-->    95CI-->',
                   'rho       -0.15962   -0.02352    0.55700    0.786520     0.9290',
                   'sig2       0.11171    0.21566    0.87405    7.196468   178.1855',
                   'beta.x1   -2.96240   -1.42573   -0.71298    0.703472     5.1099',
                   'beta.x2   -1.54479   -1.18338   -0.65587   -0.032576     5.4573'))
})
