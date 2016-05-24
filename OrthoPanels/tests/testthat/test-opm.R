context('opm')

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
    beta <- .5
    rho <- .5

    ## $x_i = 0.75 f + N(0, 1)$:
    x <- aperm(array(.75*f + rnorm(N*(T+1)), dim = c(N, 1, T+1)), 3:1)

    ## $y_{i,t} = \rho y_{i,t-1} + \beta x_{i,t} + f_i + N(0,1)$:
    y <- matrix(rep(beta*x[1,,] + f + rnorm(N), each = T+1), T+1, N)
    for (t in seq_len(T)+1) {
        y[t,] <- rho*y[t-1,] + f + beta*x[t,,] + rnorm(N)
    }

    x <<- x
    y <<- y
    d <<- data.frame(i=rep(1:3, each=3), t=rep(1:3, 3), x=c(x), y = c(y))

    set.seed(123)
    o <<- opm(y~x, d, n.samp = 10)
})


## Sanity check
test_that('data', {
    expect_equal(x, array(c(-1.31047564655221, -0.679491608575424, -0.289083794010798, 
                            -0.23017748948328, 0.129287735160946, -1.26506123460653,
                            2.30870831414912, 2.46506498688328, 0.0631471481064739),
                           dim = c(3, 1, 3)))
    expect_equal(y, matrix(c(-2.10089979337606, 1.10899305269782, 2.51416798413193,
                             -1.98942425038169, 0.729823109874504, 2.93377535075353,
                             -0.352340885393167, 0.230231415863224, 0.531844092800363),
                           nrow = 3, ncol = 3, byrow=TRUE))
})


test_that('default', {
    set.seed(123)
    expect_equal(opm(x, y, n = 10),
                 structure(list(samples = list(rho = c(0.728, 0.574,
                                                       -0.181, 0.763,
                                                       0.877, 0.965,
                                                       0.0559999999999999,
                                                       0.781, 0.102,
                                                       -0.086),
                                               sig2 = 1 / c(3.50041867260549,
                                                            4.00944598251174,
                                                            1.74928182974987,
                                                            0.0471496595688863,
                                                            0.342753347371835,
                                                            0.0355468448260593,
                                                            0.578145844500398,
                                                            1.84521258726027,
                                                            1.89471612513104,
                                                            0.815747433844173),
                                               beta = as.matrix(c(1.06781404188239,
                                                                  1.57176166272772,
                                                                  0.984064976368265,
                                                                  2.10164684546106,
                                                                  0.663384108167381,
                                                                  -0.714948757610917,
                                                                  1.49999093333423,
                                                                  1.25150458784134,
                                                                  1.11183526859818,
                                                                  1.05848623306165))),
                                call = quote(opm(x = x, y = y, n.samp = 10)),
                                time.indicators = FALSE),
                           class = 'opm'))
})


test_that('formula', {
    d <- data.frame(i = rep(1:3, each=3),
                    t = rep(1:3, 3),
                    x = c(x), y = c(y))
    expected <- structure(list(samples = list(rho = c(0.728, 0.574,
                                                      -0.181, 0.763,
                                                      0.877, 0.965,
                                                      0.0559999999999999,
                                                      0.781, 0.102,
                                                      -0.086),
                                              sig2 = 1 / c(3.50041867260549,
                                                           4.00944598251174,
                                                           1.74928182974987,
                                                           0.0471496595688863,
                                                           0.342753347371835,
                                                           0.0355468448260593,
                                                           0.578145844500398,
                                                           1.84521258726027,
                                                           1.89471612513104,
                                                           0.815747433844173),
                                              beta = as.matrix(c(1.06781404188239,
                                                                 1.57176166272772,
                                                                 0.984064976368265,
                                                                 2.10164684546106,
                                                                 0.663384108167381,
                                                                 -0.714948757610917,
                                                                 1.49999093333423,
                                                                 1.25150458784134,
                                                                 1.11183526859818,
                                                                 1.05848623306165))),
                               call = quote(opm(x = y~x, data = d, n.samp = 10)),
                               time.indicators = FALSE,
                               index = c('i', 't'),
                               terms = structure(y~x,
                                                 variables = quote(list(y, x)),
                                                 factors = matrix(0:1, 2, 1,
                                                                  dimnames=list(c('y', 'x'),
                                                                                'x')),
                                                 term.labels = 'x',
                                                 order = 1L,
                                                 intercept = 1L,
                                                 response = 1L,
                                                 .Environment = parent.frame(),
                                                 predvars = quote(list(y, x)),
                                                 dataClasses = c(y='numeric', x='numeric'),
                                                 class=c('terms', 'formula'))),
                          class = 'opm')

    ## default index is the first two columns of the data
    set.seed(123)
    expect_equal(opm(y~x, d,
                     n.samp = 10),
                 expected)

    ## look for variables in the parent frame
    x <- c(x)
    y <- c(y)
    i <- rep(1:3, each=3)
    t <- rep(1:3, 3)
    set.seed(123)
    expected$call <- quote(opm(x = y~x, index = c('i', 't'),
                               n.samp = 10))
    expect_equal(opm(y~x, index = c('i', 't'),
                     n.samp = 10),
                 expected)

    d <- d[,c('y', 'x', 't', 'i')]

    ## specify index column by position
    set.seed(123)
    expected$call <- quote(opm(x = y~x, data = d,
                               index = 4:3, n.samp = 10))
    expect_equal(opm(y~x, d, index=4:3,
                     n.samp = 10),
                 expected)
   
    ## specify index column by name
    set.seed(123)
    expected$call <- quote(opm(x = y~x, data = d,
                               index = c('i', 't'), n.samp = 10))
    expect_equal(opm(y~x, d, index=c('i', 't'), n.samp = 10),
                 expected)
   
})


test_that('dummy vals', {
    set.seed(123)
    expect_equal(opm(x, y, n = 10, add = TRUE),
                 structure(list(samples = list(rho = c(-0.647,
                                                       -0.647,
                                                       -0.647,
                                                       0.763,
                                                       -0.651,
                                                       -0.635,
                                                       -0.647,
                                                       -0.648,
                                                       -0.647,
                                                       -0.086),
                                               sig2 = 1 / c(5781053178.63559,
                                                            382064945.310797,
                                                            297081166.091305,
                                                            0.00556021027046193,
                                                            202407.002118016,
                                                            79927.7878466404,
                                                            31958664.3506085,
                                                            4096230.88017378,
                                                            172152798.137179,
                                                            25.8707402214816),
                                               beta = matrix(c(1.37506444649828,
                                                               1.3750727761096,
                                                               1.37511878023697,
                                                               3.60053789546568,
                                                               1.37414275681303,
                                                               1.37772284729159,
                                                               1.37501825167012,
                                                               1.37574557193104,
                                                               1.37512079428095,
                                                               1.51620694338847,
                                                               1.17236073947895,
                                                               1.17240636180504,
                                                               1.17250411197871, 
                                                               -8.92051368795613,
                                                               1.1693793246162,
                                                               1.17516181067838,
                                                               1.17225716029063, 
                                                               1.17220803859993,
                                                               1.17238653820271,
                                                               1.35629081137943),
                                                             10, 2,
                                                             dimnames = list(NULL, c('1', 'tind.2')))),
                                call = quote(opm(x = x, y = y, n.samp = 10, add.time.indicators = TRUE)),
                                time.indicators = TRUE),
                           class = 'opm'))
})


test_that('formula dummy vals', {
    set.seed(123)
    expect_equal(opm(y ~ x, d, n = 10, add = TRUE),
                 structure(list(samples = list(rho = c(-0.647,
                                                       -0.647,
                                                       -0.647,
                                                       0.763,
                                                       -0.651,
                                                       -0.635,
                                                       -0.647,
                                                       -0.648,
                                                       -0.647,
                                                       -0.086),
                                               sig2 = 1 / c(5781053178.63559,
                                                            382064945.310797,
                                                            297081166.091305,
                                                            0.00556021027046193,
                                                            202407.002118016,
                                                            79927.7878466404,
                                                            31958664.3506085,
                                                            4096230.88017378,
                                                            172152798.137179,
                                                            25.8707402214816),
                                               beta = matrix(c(1.37506444649828,
                                                               1.3750727761096,
                                                               1.37511878023697,
                                                               3.60053789546568,
                                                               1.37414275681303,
                                                               1.37772284729159,
                                                               1.37501825167012,
                                                               1.37574557193104,
                                                               1.37512079428095,
                                                               1.51620694338847,
                                                               1.17236073947895,
                                                               1.17240636180504,
                                                               1.17250411197871, 
                                                               -8.92051368795613,
                                                               1.1693793246162,
                                                               1.17516181067838,
                                                               1.17225716029063, 
                                                               1.17220803859993,
                                                               1.17238653820271,
                                                               1.35629081137943),
                                                             10, 2,
                                                             dimnames = list(NULL, c('1', 'tind.2')))),
                                call = quote(opm(x = y ~ x, data = d, n.samp = 10, add = TRUE)),
                                time.indicators = TRUE,
                                index = c('i', 't'),
                                terms = structure(y~x,
                                                  variables = quote(list(y, x)),
                                                  factors = matrix(0:1, 2, 1,
                                                                   dimnames=list(c('y', 'x'),
                                                                                 'x')),
                                                  term.labels = 'x',
                                                  order = 1L,
                                                  intercept = 1L,
                                                  response = 1L,
                                                  .Environment = parent.frame(),
                                                  predvars = quote(list(y, x)),
                                                  dataClasses = c(y='numeric', x='numeric'),
                                                  class=c('terms', 'formula'))),
                           class = 'opm'))
})


test_that('confint', {
    expect_equal(confint(o),
                 matrix(c(-0.159625, 0.257571563455217, -0.4048238628108,
                          0.9452, 26.5742550057164, 1.98242267934606),
                        nrow = 3, ncol = 2,
                        dimnames = list(c("rho", "sig2", "beta"),
                                        c("2.5%", "97.5%"))))

    expect_equal(confint(o, 1:2, level = 0.68),
                 matrix(c(-0.0235200000000001, 0.39220562432823,
                          0.83476, 13.1607960528734),
                        nrow = 2, ncol = 2,
                        dimnames = list(c("rho", "sig2"), c("16%", "84%"))))
})


test_that('coef', {
    expect_equal(coef(o),
                 c(rho = 0.651,
                   sig2 = 0.898766398134737,
                   beta = 1.08982465524029))

    expect_error(coef(o, probs = .6), "Arguments 'probs' and 'names' are not allowed")
    expect_error(coef(o, names = TRUE), "Arguments 'probs' and 'names' are not allowed")
})


test_that('print', {
    expect_equal(capture.output(o),
                 c('Call:',
                   'opm(x = y ~ x, data = d, n.samp = 10)',
                   '',
                   'Coefficients:',
                   '           mean (SD)      med          95-CI',
                   'rho    0.4579 (0.44)  0.65100  (-0.16, 0.95)',
                   'sig2  5.7391 (10.14)  0.89877  (0.26, 26.57)',
                   'beta   1.0596 (0.74)  1.08982    (-0.4, 2.0)',
                   ''))
})


test_that('summary', {
    expect_equal(capture.output(summary(o)),
                 c('Call:',
                   'opm(x = y ~ x, data = d, n.samp = 10)',
                   '',
                   'Parameter estimates:',
                   '        <--95CI    <--68CI       med    68CI-->   95CI-->',
                   'rho    -0.15962   -0.02352   0.65100    0.83476    0.9452',
                   'sig2    0.25757    0.39221   0.89877   13.16080   26.5743',
                   'beta   -0.40482    0.80448   1.08982    1.54018    1.9824'))
})
