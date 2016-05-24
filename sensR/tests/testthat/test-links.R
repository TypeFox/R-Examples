## library(testthat)
context("Accuracy of binomial link functions")

##################################################################
## Tests of triangle family object:
delta.lim <- c(-Inf, -1, 0, 1e-3, 1e-2, .1, 1:20, 30, Inf)
trimat <- cbind(delta.lim,
                tri.linkinv=triangle()$linkinv(delta.lim),
                muEta=triangle()$mu.eta(delta.lim))
pc <- c(-Inf, -1, 0:3/10, 1/3, 1/3 + 1e-9, 1/3 + 1e-7,
        4:9/10, 1-1e-7, 1-1e-9, 1, 2, Inf)
trimat2 <- cbind(pc, triangle()$linkfun(pc))

trimat.expect <-
    structure(c(-Inf, -1, 0, 0.001, 0.01, 0.1, 1, 2, 3, 4, 5, 6,
                7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 30, Inf,
                0.333333333333333, 0.333333333333333, 0.333333333333333, 0.333333425221489,
                0.33334252221233, 0.334251449601924, 0.418046674796045, 0.604806932920651,
                0.781427598857822, 0.897660380694802, 0.958777656024149, 0.985694205716359,
                0.995733276661051, 0.998909165425241, 0.99976143708784, 0.999955443364009,
                0.999992902702296, 0.999999036983941, 0.999999889054218, 0.999999989273082,
                0.999999999198863, 0.999999999963858, 0.999999999999226, 0.999999999999993,
                1, 1, 1, 1, 0, 0, 0, 0.000183776267844556, 0.00183773235566254,
                0.0183470310702133, 0.155988959563108, 0.196685802264629, 0.148651046331582,
                0.0854611723784497, 0.0405422168900257, 0.0162170328596508, 0.00548902209724158,
                0.00157263230530614, 0.000381396482494682, 7.82967533087891e-05,
                1.36059355818558e-05, 2.00138505966741e-06, 2.49201402860651e-07,
                2.62656346212414e-08, 2.34338097972049e-09, 1.76976453975268e-10,
                1.13137261865862e-11, 6.12228333969946e-13, 2.80439234112079e-14,
                0, 0, 0), .Dim = c(28L, 3L), .Dimnames = list(NULL, c("delta.lim",
                                             "tri.linkinv", "muEta")))
trimat2.expect <-
    structure(c(-Inf, -1, 0, 0.1, 0.2, 0.3, 0.333333333333333, 0.333333334333333,
                0.333333433333333, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.9999999, 0.999999999,
                1, 2, Inf, 0, 0, 0, 0, 0, 0, 0, 0.0001220703125, 0.00101989522980266,
                0.879114569369712, 1.46626281186262, 1.97560154484988, 2.50372333631553,
                3.12859127911218, 4.02762596771379, 13.0460576380004, Inf, Inf,
                Inf, Inf), .Dim = c(20L, 2L), .Dimnames = list(NULL, c("pc",
                                              "")))

test_that("Triangle linkfun, linkinv and deriv return what we expect", {
    expect_equal(trimat, trimat.expect)
    expect_equal(trimat2, trimat2.expect)

    pcc <- triangle()$linkinv(14)
    expect_true(pcc < 1)
    expect_equal(triangle()$linkfun(pcc), 14)
    pcc <- triangle()$linkinv(14 + .1)
    expect_true(pcc <= 1)
    expect_true(triangle()$linkfun(pcc) == Inf)
    pcc <- triangle()$linkinv(15)
    expect_true(pcc <= 1)
    expect_true(triangle()$linkfun(pcc) == Inf)
})

##################################################################
#################################
## Tests of 3-AFC family object:

## Evaluating linkinv, mu.eta and linkfun for a range of values:
delta.lim <- c(-Inf, -1, 0, 1e-3, 1e-2, .1, 1:20, 30, Inf)
pc <- c(-Inf, -1, 0:3/10, 1/3, 1/3 + 1e-4, 1/3 + 1e-3,
        4:9/10, 1-1e-7, 1-1e-8, 1, 2, Inf)
threeMat <-
    cbind(delta.lim,
          linkinv= threeAFC()$linkinv(delta.lim),
          muEta= threeAFC()$mu.eta(delta.lim))
threeMat2 <- cbind(pc, threeAFC()$linkfun(pc))

threeMat.expect <-
    structure(c(-Inf, -1, 0, 0.001, 0.01, 0.1, 1, 2, 3, 4, 5, 6,
7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 30, Inf,
0.333333333333333, 0.333333333333333, 0.333333333333333, 0.333615473942033,
0.336158851977074, 0.361978125272024, 0.633702050336587, 0.865767175634742,
0.968795477969601, 0.995496522284447, 0.999599192120491, 0.999978025840482,
0.999999258053642, 0.999999984588949, 1, 1, 1, 1, 1, 1, 1, 1,
1, 1, 1, 1, 1, 1, 0, 0, 0, 0.282186609373901, 0.283006595460479,
0.290553761744979, 0.28931908430394, 0.164567688424708, 0.0529040173850969,
0.00980405393282727, 0.0010666911724293, 6.9128491660797e-05,
2.69395392971106e-06, 6.34564996659925e-08, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0), .Dim = c(28L, 3L), .Dimnames = list(
    NULL, c("delta.lim", "linkinv", "muEta")))
threeMat2.expect <-
    structure(c(-Inf, -1, 0, 0.1, 0.2, 0.3, 0.333333333333333, 0.333433333333333,
0.334333333333333, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.9999999, 0.99999999,
1, 2, Inf, 0, 0, 0, 0, 0, 0, 0, 0.000354412995510282, 0.00353718447587665,
0.228853730602745, 0.556523357673813, 0.885182389787633, 1.23795339238868,
1.65241120055872, 2.23019992157, 7.53291713420153, Inf, Inf,
Inf, Inf), .Dim = c(20L, 2L), .Dimnames = list(NULL, c("pc",
"")))


test_that("3-AFC linkfun, linkinv and deriv return what we expect", {
    expect_equal(threeMat, threeMat.expect)
    expect_equal(threeMat2, threeMat2.expect)

    ## Going smooth toward zero:
    x <- threeAFC()$mu.eta(8 - 10^(0:-7))
    expect_true(all(diff(x) < 0))

    ## Upper limit for linkinv:
    expect_true(1 - threeAFC()$linkinv(9-.1) > 0) ## 3.102338e-10
    expect_true(1 - threeAFC()$linkinv(9) == 0) # 0
    expect_true(1 - threeAFC()$linkinv(9+.1) == 0) # 0
    ## Lower limit for linkinv:
    expect_true((threeAFC()$linkinv(1e-9) - 1/3) > 0)# 1.791433e-10
    expect_true(threeAFC()$linkinv(1e-10) == 1/3) # 0

    ## The upper limit for linkfun:
    x <- threeAFC()$linkfun(1 - 1e-7) ## > Inf ~ 7.53
    ## dput(x)
    expect_equal(x, 7.53291713420153)
    x <- threeAFC()$linkfun(1 - 1.1e-8) ## < Inf ~ 8.0819
    expect_equal(x, 8.08189980744933)
    expect_true(threeAFC()$linkfun(1 - 1e-8) == Inf) ## Inf
    ## Lower limit for linkfun:
    expect_true(threeAFC()$linkfun(1/3 + 1e-5) > 0)# 6.103516e-05
    expect_true(threeAFC()$linkfun(1/3 + 1e-6)  == 0)# 0
})

##################################################################
## Tests of tetrad family object:

## Evaluating linkinv, mu.eta and linkfun for a range of values:
delta.lim <- c(-Inf, -1, 0, 1e-3, 1e-2, .1, 1:20, 30, Inf)
pc <- c(-Inf, -1, 0:3/10, 1/3, 1/3 + 1e-4, 1/3 + 1e-3,
        4:9/10, 1-1e-7, 1-1e-8, 1, 2, Inf)
tetradMat <- cbind(delta.lim,
            linkinv= tetrad()$linkinv(delta.lim),
            muEta= tetrad()$mu.eta(delta.lim))
tetradMat2 <- cbind(pc, tetrad()$linkfun(pc))

## dput(tetradMat)
tetradMat.expect <-
    structure(c(-Inf, -1, 0, 0.001, 0.01, 0.1, 1, 2, 3, 4, 5, 6,
7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 30, Inf,
0.333333333333333, 0.333333333333333, 0.333333333333333, 0.333333517315509,
0.333351710913841, 0.335168546455991, 0.493808427692137, 0.777667116845886,
0.942971619101258, 0.991341559084544, 0.999210674034481, 0.999956060315649,
0.999998497688177, 0.999999970097948, 1, 1, 1, 1, 1, 1, 1, 1,
1, 1, 1, 1, 1, 1, 0, 0, 0, 0.000367556882833047, 0.00367542471210065,
0.0366533045258264, 0.278493758283903, 0.243163256289553, 0.0926857799809961,
0.0185492303766036, 0.00208847744671707, 0.000137347787572007,
5.46406146954242e-06, 1.25374155477369e-07, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0), .Dim = c(28L, 3L), .Dimnames = list(
    NULL, c("delta.lim", "linkinv", "muEta")))
tetradMat2.expect <-
    structure(c(-Inf, -1, 0, 0.1, 0.2, 0.3, 0.333333333333333, 0.333433333333333,
0.334333333333333, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.9999999, 0.99999999,
1, 2, Inf, 0, 0, 0, 0, 0, 0, 0, 0.0233423654546829, 0.073791381662713,
0.618329709770766, 1.02212226454153, 1.36259707370522, 1.70450506671692,
2.09460011556416, 2.63258173217749, 7.70657597572176, Inf, Inf,
Inf, Inf), .Dim = c(20L, 2L), .Dimnames = list(NULL, c("pc",
"")))

test_that("Tetrad linkfun, linkinv and deriv return what we expect", {
    expect_equal(tetradMat, tetradMat.expect)
    expect_equal(tetradMat2, tetradMat2.expect)

    ## Going smooth toward zero:
    x <- tetrad()$mu.eta(8 - 10^(0:-7))
    expect_true(all(diff(x) < 0))

    ## Upper limit for linkinv:
    expect_true((1 - tetrad()$linkinv(9-.1)) > 0) # 5.737801e-10
    expect_true((1 - tetrad()$linkinv(9)) == 0) # 0
    expect_true((1 - tetrad()$linkinv(9+.1)) == 0) # 0

    ## Lower limit for linkinv:
    tetrad()$linkinv(1e-4) - 1/3 # 2.043666e-09
    tetrad()$linkinv(1e-5) - 1/3 # 2.242808e-10
    tetrad()$linkinv(1e-6) - 1/3 # 2.060869e-10
    expect_true((tetrad()$linkinv(1e-7) - 1/3) > 0) # 2.059049e-10
    expect_true((tetrad()$linkinv(1e-8) - 1/3) == 0) # 0

    ## The upper limit for linkfun:
    tetrad()$linkfun(1 - 1e-7) ## > Inf ~ 7.706576
    expect_true(tetrad()$linkfun(1 - 1.1e-8) < Inf) ## < Inf ~ 8.235212
    expect_true(tetrad()$linkfun(1 - 1e-8) == Inf) ## Inf

    ## Lower limit for linkfun:
    tetrad()$linkfun(1/3 + 1e-8) # 0.0002441406
    expect_true(tetrad()$linkfun(1/3 + 1e-9) > 0) # 6.103516e-05
    expect_true(tetrad()$linkfun(1/3 + 1e-10) == 0)# 0
    expect_true(tetrad()$linkfun(1/3) == 0) # 0

    ## mu.eta - upper limit:
    expect_true(tetrad()$mu.eta(8) > 0)
    expect_equal(tetrad()$mu.eta(9), 0)

    ## mu.eta - lower limit
    expect_equal(tetrad()$mu.eta(0), 0)
    expect_equal(tetrad()$mu.eta(1e-8), 0)
    expect_true(tetrad()$mu.eta(1e-7) > 0)
})


#################################
## Tests of duo-trio family object:

## Evaluating linkinv, mu.eta and linkfun for a range of values:
delta.lim <- c(-Inf, -1, 0, 1e-3, 1e-2, .1, 1:20, 30, Inf)
pc <- c(-Inf, -1, 0:5/10, 0.5 + 1e-9, 0.5 + 1e-7,
        6:9/10, 1-1e-9, 1-1e-10, 1, 2, Inf)
duotrioMat <- cbind(delta.lim,
                    DT.linkinv=duotrio()$linkinv(delta.lim),
                    DT.mu.eta=duotrio()$mu.eta(delta.lim))
duotrioMat2 <-
    cbind(pc, duotrio()$linkfun(pc))
## dput(duotrioMat2)
duotrioMat.expect <-
    structure(c(-Inf, -1, 0, 0.001, 0.01, 0.1, 1, 2, 3, 4, 5, 6,
7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 30, Inf,
0.5, 0.5, 0.5, 0.500000091888139, 0.500009188712827, 0.500917861363183,
0.582475444228876, 0.746820245546069, 0.876456704022415, 0.946665579848435,
0.979191495994195, 0.992836173545768, 0.997866267625024, 0.999454574711717,
0.999880718174701, 0.999977721453929, 0.999996451045831, 0.999999518321496,
0.999999944347472, 0.999999994530134, 0.999999999542935, 0.999999999967546,
0.999999999998043, 0.9999999999999, 0.999999999999996, 1, 1,
1, 0, 0, 0, 0.000183776257634716, 0.00183772214607225, 0.0183368416735807,
0.147617919351631, 0.159133715202119, 0.0974969682475734, 0.0473678937836294,
0.0207932285623696, 0.00814283165850389, 0.00274585514584807,
0.000786347863610548, 0.000190698693965787, 3.91483805719366e-05,
6.80296781148603e-06, 1.00069252989914e-06, 1.24600701430452e-07,
1.31328173106209e-08, 1.17169048986025e-09, 8.84882269876348e-11,
5.65686309329313e-12, 3.06114166984973e-13, 1.40219617056039e-14,
5.43690457381535e-16, 4.36262916005144e-34, 0), .Dim = c(28L,
3L), .Dimnames = list(NULL, c("delta.lim", "DT.linkinv", "DT.mu.eta"
)))
duotrioMat2.expect <-
    structure(c(-Inf, -1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.500000001,
0.5000001, 0.6, 0.7, 0.8, 0.9, 0.999999999, 0.9999999999, 1,
2, Inf, 0, 0, 0, 0, 0, 0, 0, 0, 0.0001220703125, 0.00101652663801013,
1.11520955683069, 1.71530153283467, 2.35486196368808, 3.26308047458777,
14.6915420837302, Inf, Inf, Inf, Inf), .Dim = c(19L, 2L), .Dimnames = list(
    NULL, c("pc", "")))

test_that("Duo-trio linkfun, linkinv and deriv return what we expect", {
    expect_equal(duotrioMat, duotrioMat.expect)
    expect_equal(duotrioMat2, duotrioMat2.expect)

    ## linkinv - upper limit:
    expect_true(1 - duotrio()$linkinv(19) > 0)
    expect_equal(1 - duotrio()$linkinv(20), 0)
    ## eta > 20 -> mu = 1

    ## linkinv - lower limit:
    expect_equal((duotrio()$linkinv(1e-8) - .5), 0)
    expect_true((duotrio()$linkinv(1e-6) - .5) > 0)

    ## mu.eta - upper limit:
    expect_true(duotrio()$mu.eta(94) > 0)
    expect_equal(duotrio()$mu.eta(95), 0)

    ## mu.eta - lower limit
    expect_equal(duotrio()$mu.eta(0), 0)
    expect_equal(duotrio()$mu.eta(1e-16), 0)
    expect_true(duotrio()$mu.eta(1e-15) > 0)

    ## linkfun - upper limit
    duotrio()$linkfun(1 - 1e-8)
    duotrio()$linkfun(1 - 1e-9)
    expect_true(duotrio()$linkfun(1 - 1.000001e-10) < Inf)
    expect_true(duotrio()$linkfun(1 - 1e-10) == Inf)

    ## linkfun - lower limit:
    expect_true(duotrio()$linkfun(.5) == 0)
    expect_true(duotrio()$linkfun(.5 + 1e-10)  == 0)
    expect_true(duotrio()$linkfun(.5 + 1.7e-10)  == 0)
    expect_true(duotrio()$linkfun(.5 + 1.8e-10) > 0)
    duotrio()$linkfun(.5 + 5e-10)
    duotrio()$linkfun(.5 + 1e-9)
    duotrio()$linkfun(.5 + 1e-8)
    duotrio()$linkfun(.6)
})

##################################################################
## Tests of the 2-AFC family object:

## Evaluating linkinv, mu.eta and linkfun for a range of values:
delta.lim <- c(-Inf, -1, 0, 1e-3, 1e-2, .1, 1:20, 30, Inf)
pc <- c(-Inf, -1, 0:5/10, 0.5 + 1e-9, 0.5 + 1e-7,
        6:9/10, 1-1e-9, 1-1e-10, 1, 2, Inf)
twoafcMat <-
    cbind(delta.lim,
          DT.linkinv=twoAFC()$linkinv(delta.lim),
          DT.mu.eta=twoAFC()$mu.eta(delta.lim))
twoafcMat2 <-
    cbind(pc, twoAFC()$linkfun(pc))
## dput(twoafcMat2)
twoafcMat.expect <-
    structure(c(-Inf, -1, 0, 0.001, 0.01, 0.1, 1, 2, 3, 4, 5, 6,
7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 30, Inf,
0.5, 0.5, 0.5, 0.500282094768266, 0.502820924410016, 0.528185988898508,
0.760249938906523, 0.921350396474857, 0.983052573237655, 0.997661132509476,
0.999796523991277, 0.999988954751501, 0.999999628450814, 0.999999992291371,
0.999999999901692, 0.999999999999231, 0.999999999999996, 1, 1,
1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0.282094791773878, 0.282094721250189,
0.282087739492238, 0.281390435606505, 0.219695644733861, 0.103776874355149,
0.0297325723059073, 0.00516674633852301, 0.000544571057588177,
3.48132629866869e-05, 1.3498566943462e-06, 3.17455866796663e-08,
4.52826473977171e-10, 3.91771663275432e-12, 2.05582901131572e-14,
6.54325309812312e-17, 1.26314500076565e-19, 1.47899073950077e-22,
1.05034134452874e-25, 4.5242669921399e-29, 1.18200346721402e-32,
1.87301835706147e-36, 1.80018888638499e-40, 1.04941405783859e-44,
5.42171444080732e-99, 0), .Dim = c(28L, 3L), .Dimnames = list(
    NULL, c("delta.lim", "DT.linkinv", "DT.mu.eta")))
twoafcMat2.expect <-
    structure(c(-Inf, -1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.500000001,
0.5000001, 0.6, 0.7, 0.8, 0.9, 0.999999999, 0.9999999999, 1,
2, Inf, 0, 0, 0, 0, 0, 0, 0, 0, 3.5449076015542e-09, 3.54490769994519e-07,
0.358286909242583, 0.741614317187116, 1.19023216289999, 1.81238760487365,
8.48218003161719, 8.99629456108862, Inf, Inf, Inf), .Dim = c(19L,
2L), .Dimnames = list(NULL, c("pc", "")))

test_that("2-AFC linkfun, linkinv and deriv return what we expect", {
    expect_equal(twoafcMat, twoafcMat.expect)
    expect_equal(twoafcMat2, twoafcMat2.expect)

    ## linkinv - upper limit:
    expect_true((1 - twoAFC()$linkinv(11)) > 0)
    expect_true((1 - twoAFC()$linkinv(12)) == 0)
    ## eta > 11 -> mu = 1

    ## linkinv - lower limit:
    expect_true((twoAFC()$linkinv(1e-16) - .5) == 0)
    expect_true((twoAFC()$linkinv(1e-15) - .5) > 0)

    ## mu.eta - upper limit:
    expect_true(twoAFC()$mu.eta(54) > 0)
    expect_true(twoAFC()$mu.eta(55) == 0)

    ## mu.eta - lower limit
    expect_equal(twoAFC()$mu.eta(0), 0.282094791773878)

    ## linkfun - upper limit
    expect_true(twoAFC()$linkfun(1 - 1e-16) < Inf)
    expect_true(twoAFC()$linkfun(1 - 1e-17) == Inf)

    ## linkfun - lower limit:
    expect_true(twoAFC()$linkfun(.5) == 0)
    expect_true(twoAFC()$linkfun(.5 + 1e-17) == 0)
    expect_true(twoAFC()$linkfun(.5 + 1e-16) > 0)
})
