

test.pcountOpen.null <- function()
{
  y <- matrix(c(
      3, 2, 1, 4,
      3, 4, 2, 1,
      0, 1, 2, 3,
      5, 3, 3, 4,
      2, 4, 3, 3), 5, 4, byrow=TRUE)
  siteCovs <- data.frame(x = c(0,2,3,4,1))
  obsCovs <- data.frame(o1 = 1:20)
  umf <- unmarkedFramePCO(y = y, siteCovs = siteCovs, obsCovs = obsCovs,
      numPrimary=4)

  fm1 <- pcountOpen(~1, ~1, ~1, ~1, data = umf, se=FALSE, K=10,
      starts=c(1, 0, 0, 7))
  checkEqualsNumeric(coef(fm1),
                     c(0.9565311, 0.2741022, 0.1352888, 7.0041290),
                     tol = 1e-5)

  fm2 <- pcountOpen(~1, ~1, ~1, ~1, data = umf, se=FALSE, fix="gamma",
                    K=10)
  checkEqualsNumeric(coef(fm2), c(1.8219354, 8.7416638, -0.2873611),
      tol = 1e-3)

  fm3 <- pcountOpen(~1, ~1, ~1, ~1, data = umf, se=FALSE, fix="omega",
                    K=30)
  checkEqualsNumeric(coef(fm3), c(1.8861091, -1.3102890, -0.4934883),
      tol = 1e-5)

}






test.pcountOpen.na <- function()
{

  y1 <- matrix(c(
      NA, 2, 1, 4,
      3, NA, 2, 1,
      0, 1, 2, NA,
      5, NA, 3, NA,
      NA, NA, 3, NA,
      2, NA, NA, NA), 6, 4, byrow=TRUE)
  siteCovs <- data.frame(x = c(0,2,3,4,1,1))
  obsCovs <- data.frame(o1 = 1:24)
  umf1 <- unmarkedFramePCO(y = y1, siteCovs = siteCovs, obsCovs = obsCovs,
      numPrimary=4)

  fm1 <- pcountOpen(~1, ~1, ~1, ~1, data = umf1, se=FALSE, K=10,
      starts=c(1.6, 0.24, 1.16, -0.268))
  checkEqualsNumeric(coef(fm1),
      c(0.9536852, 0.4611643, -0.8655834, 3.2154420), tol = 1e-5)

  y2 <- matrix(c(
      1, 2, 1, 4,
      3, 1, 2, 1,
      0, 1, 2, 1,
      5, 1, 3, 1,
      1, 1, 3, 1,
      2, 1, 1, 1), 6, 4, byrow=TRUE)
  oc <- y1 + -2:3

  siteCovs <- data.frame(x = c(0,2,3,4,1,1))
  obsCovs <- list(o1 = oc)
  ysc <- list(o2=oc)
  umf2 <- unmarkedFramePCO(y = y2, siteCovs = siteCovs, obsCovs = obsCovs,
      yearlySiteCovs=ysc, numPrimary=4)

  fm2.1 <- pcountOpen(~1, ~1, ~1, ~o1, data = umf2, se=FALSE, K=10,
      starts=c(1.4, -1.3, 1.8, -1.1, 0.7))
  checkEqualsNumeric(coef(fm2.1),
      c(1.239636, -2.085200, 1.770919, -0.602612, 1.255386),
                     tol = 1e-4)

  fm2.2 <- pcountOpen(~1, ~1, ~o2, ~1, data = umf2, se=FALSE, K=10,
      starts=c(1.2, -1, 2, 0, 0))
  checkEqualsNumeric(coef(fm2.2),
      c(1.3242059, 0.8439311, -2.8217070, -10.1414153, 0.1176959),
      tol = 1e-5)

  fm2.3 <- pcountOpen(~1, ~o2, ~1, ~1, data = umf2, se=FALSE, K=10,
      starts=c(1, 0, 0, -5, -1))
  checkEqualsNumeric(coef(fm2.3),
      c(0.7013386, 0.5277811, -0.2350951, -1.8346326, 4.7771974),
                     tol = 1e-2)

  y3 <- matrix(c(
      NA, 2, 1, 4,
      3, NA, 2, 1,
      0, 1, 2, NA,
      5, NA, 3, NA,
      NA, NA, 3, NA,
      NA, NA, NA, NA), 6, 4, byrow=TRUE)
  siteCovs <- data.frame(x = c(0,2,3,4,1,1))
  obsCovs <- data.frame(o1 = 1:24)
  umf3 <- unmarkedFramePCO(y = y3, siteCovs = siteCovs, obsCovs = obsCovs,
      numPrimary=4)

  fm3 <- pcountOpen(~1, ~1, ~1, ~1, data = umf3, se=FALSE, K=10,
      starts=c(1.5, 0, 1, 0))
  checkEqualsNumeric(coef(fm3),
                     c(0.9714456, 0.4481042, -0.8920462, 4.0379739 ),
                     tol = 1e-4)
  checkEquals(fm3@sitesRemoved, 6)



  y4 <- matrix(c(
      NA, 2, 1, 4,
      3, 1, 2, 1,
      0, 1, 2, 1,
      5, 1, 3, 1,
      1, 1, 3, 1,
      2, 1, 1, 1), 6, 4, byrow=TRUE)
  go4 <- matrix(c(
      NA, NA, NA, # remove y[1, 2:4]
      NA, 1, 2,   # remove y[2, 2]
      0, NA, 2,   # remove y[3, 3]. Creates an interior NA
      5, 1, NA,   # remove y[4, 4]. Creates an end NA
      1, NA, NA,  # remove y[5, 3:4]
      NA, NA, 1), 6, 3, byrow=TRUE)
  o2y <- matrix(c(
      1, 0, 0,
      0, 1, 0,
      0, 0, 1), 3, 3, byrow=TRUE)
  y4.na <- is.na(go4) %*% o2y
  y4.2 <- y4
  y4.2[,-1][y4.na>0] <- NA
  y4.2

  umf4 <- unmarkedFramePCO(y=y4, yearlySiteCovs=list(go4=cbind(go4, 1)),
      numPrimary=4)

  fm4.1 <- pcountOpen(~1, ~go4, ~1, ~1, umf4, se=FALSE,
      starts=c(.8, .5, -.3, -1.5, 6))
  checkEquals(fm4.1@sitesRemoved, 1)

  fm4.2 <- pcountOpen(~1, ~1, ~go4, ~1, umf4, se=FALSE,
      starts=c(.8, 0, 5, -5, 7))
  checkEquals(fm4.2@sitesRemoved, 1)



    # Now with secondary sampling periods

    y5 <- matrix(c(
        2,2,  2,2,  1,0,
        3,2,  1,1,  2,2,
        0,0,  1,1,  1,1,
        3,3,  2,0,  3,3), 4, 6, byrow=TRUE)

    umf5 <- unmarkedFramePCO(y=y5, numPrimary=3)
    fm5 <- pcountOpen(~1, ~1, ~1, ~1, umf5, se=FALSE, K=10)
    checkEqualsNumeric(coef(fm5),
        c(0.7269958, -0.3484145, 0.1494188, 1.9391898), tol=1e-5)

    y6 <- y5
    y6[1,1] <- y6[2,3:4] <- y6[3,5:6] <- y6[4,6] <- NA
    umf6 <- unmarkedFramePCO(y=y6, numPrimary=3)
    fm6 <- pcountOpen(~1, ~1, ~1, ~1, umf6, se=FALSE, K=10)
    checkEqualsNumeric(coef(fm6),
        c(0.7945817, -0.4340502, 0.5614526, 1.4161393), tol=1e-5)

    y7 <- y5
    oc7 <- y6 + -2:1
    umf7 <- unmarkedFramePCO(y=y7, obsCovs=list(oc=oc7), numPrimary=3)
    fm7 <- pcountOpen(~1, ~1, ~1, ~oc, umf7, se=FALSE, K=10)
    checkEqualsNumeric(coef(fm7),
        c(1.1986029, -9.9298367, 12.5760064, -0.8876606, 0.9525870),
                       tol=1e-4)

    y8 <- y5
    ysc8 <- matrix(1:3, #rnorm(12),
                   4, 3, byrow=TRUE)
    ysc8[2,1] <- NA # note this will fail if it's ysc8[1,1] <- NA
    umf8 <- unmarkedFramePCO(y=y8, yearlySiteCovs=list(ysc=ysc8),
                             numPrimary=3)
  ## !!!!!!!!!!!
    fm8 <- pcountOpen(~1, ~1, ~ysc, ~1, umf8, se=FALSE, K=10)
    checkEqualsNumeric(coef(fm8),
        c(0.7278796, -0.8770411, 0.9170578, 0.0399341, 1.8956210),
                       tol=1e-4)

}








test.pcountOpen.delta <- function()
{
    M <- 5
    T <- 4
    y <- matrix(c(
        NA, 2, 1, 4,
        3, NA, 2, 1,
        0, 1, NA, 3,
        5, 3, 3, NA,
        NA, 4, NA, NA), M, T, byrow=TRUE)
    if(!exists("formatDelta"))
        formatDelta <- unmarked:::formatDelta
    dates <- matrix(c(1,3,5,7), M, T, byrow=TRUE)
    delta <- formatDelta(dates, is.na(y))
    ans <- matrix(c(
        1, 2, 2, 2,
        1, 2, 4, 2,
        1, 2, 2, 4,
        1, 2, 2, 2,
        1, 2, 2, 2), M, T, byrow=TRUE)

    checkEquals(delta, ans)

    dates2 <- matrix(c(
      2, 4, 6, 8,
      1, 4, 6, 8,
      2, 4, 6, 8,
      1, 4, 6, 8,
      2, 4, 6, 8), M, T, byrow=TRUE)
    delta2 <- formatDelta(dates2, is.na(y))
    ans2 <- matrix(c(
        2, 3, 2, 2,
        1, 3, 5, 2,
        2, 2, 2, 4,
        1, 3, 2, 2,
        2, 3, 2, 2), M, T, byrow=TRUE)

    checkEquals(delta2, ans2)

    dates3 <- matrix(as.integer(c(
      2, NA, 6, 8,
      1, 4, 6, 8,
      2, 4, 6, 8,
      1, 4, 6, 8,
      2, 4, 6, 8)), M, T, byrow=TRUE)
    checkException(unmarkedFramePCO(y=y, primaryPeriod=dates3,
                                    numPrimary=4))

    dates4 <- dates2
#    dates4[is.na(y)] <- NA
    mode(dates4) <- "integer"
    delta4 <- formatDelta(dates4, is.na(y))
    umf <- unmarkedFramePCO(y=y, primaryPeriod=dates4, numPrimary=4)
    fm <- pcountOpen(~1, ~1, ~1, ~1, umf, K=10, starts=c(1.2, 0, 1.4, 1.2))
    checkEqualsNumeric(coef(fm),
        c(1.35779948, 0.11911809, -0.06946651, 5.78090618),
#        c(1.2543989, -0.5429887, 0.6715887, 5.6593500),
    tol = 1e-5)

    y5 <- matrix(c(
        1, NA, 1, 4,
        NA, 3, 2, 1,
        0, 1, 2, NA,
        NA, NA, 3, NA,
        NA, 4, NA, NA), M, T, byrow=TRUE)
    dates5 <- matrix(c(
        2, NA, 6, 8,
        NA, 4, 6, 8,
        2, 4, 6, NA,
        NA, NA, 6, NA,
        2, 4, 6, 8), M, T, byrow=TRUE)
    ans5 <- matrix(c(
        1, NA, 4, 2,
        NA, 2, 2, 2,
        1, 2, 2, NA,
        NA, NA, 4, NA, # 4 not 5 b/c primary period 1 is day 2
        1, 2, 2, 2), M, T, byrow=TRUE)
    delta5 <- formatDelta(dates5, is.na(y5))
    checkEquals(delta5, ans5)

    dates6 <- y6 <- matrix(c(2L, 1L), 1, 2)
    checkException(unmarkedFramePCO(y=y6, primaryPeriod=dates6,
                                    numPrimary=2))



}






test.pcountOpen.secondSamps <- function()
{
    y <- matrix(c(
        0,0,  2,2,  3,2,  2,2,
        2,2,  2,1,  3,2,  1,1,
        1,0,  1,1,  0,0,  0,0,
        0,0,  0,0,  0,0,  0,0), nrow=4, ncol=8, byrow=TRUE)

    sc <- data.frame(x1 = 1:4, x2 = c('A','A','B','B'))

    oc <- list(
        x3 = matrix(1:8, nrow=4, ncol=8, byrow=TRUE),
        x4 = matrix(letters[1:8], nrow=4, ncol=8, byrow=TRUE))

    ysc <- list(
        x5 = matrix(c(
            1,2,3,4,
            1,2,3,4,
            1,2,3,4,
            1,2,3,4), nrow=4, ncol=4, byrow=TRUE))

    umf1 <- unmarkedFramePCO(y=y, siteCovs=sc, obsCovs=oc,
        yearlySiteCovs=ysc, numPrimary=4)

    m1 <- pcountOpen(~1, ~1, ~1, ~1, umf1, K=10)
    checkEqualsNumeric(coef(m1),
        c(-0.2438797, -0.7838448, 0.5572557, 1.6925454), tol=1e-5)


    y2 <- y
    y2[1,1] <- NA
    umf2 <- unmarkedFramePCO(y=y2, siteCovs=sc, obsCovs=oc,
        yearlySiteCovs=ysc, numPrimary=4)

    m2 <- pcountOpen(~1, ~1, ~1, ~1, umf2, K=10)



}













test.pcountOpen.dynamics <- function()
{

    set.seed(3)
    M <- 20
    T <- 5
    lambda <- 4
    gamma <- .2
    omega <- 0.7
    p <- 0.7
    y <- N <- matrix(NA, M, T)
    S <- G <- matrix(NA, M, T-1)
    N[,1] <- rpois(M, lambda)
    for(t in 1:(T-1)) {
	S[,t] <- rbinom(M, N[,t], omega)
	G[,t] <- rpois(M, gamma*N[,t])
	N[,t+1] <- S[,t] + G[,t]
    }
    y[] <- rbinom(M*T, N, p)
    colMeans(y)
    umf <- unmarkedFramePCO(y = y, numPrimary=T)

    m1 <- pcountOpen(~1, ~1, ~1, ~1, umf, K=20, dynamics="autoreg")
    checkEqualsNumeric(coef(m1),
                       c(1.5457081, -1.2129776,  0.5668830,  0.4987492),
                       tolerance=1e-5)

    m2 <- pcountOpen(~1, ~1, ~1, ~1, umf, K=20, dynamics="notrend")
    checkEqualsNumeric(coef(m2),
                       c(1.2131713,  0.7301736,  1.1949289),
                       tolerance=1e-5)

    m3 <- pcountOpen(~1, ~1, ~1, ~1, umf, K=20, dynamics="trend")
    checkEqualsNumeric(coef(m3),
                       c(1.67211946, -0.06534021, 0.18287762),
                       tolerance=1e-5)

}
