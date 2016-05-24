ARMARECURSION1 <- function (iStart, iEnd, phi1, theta1,
    yzeroadj, innov, lny2adj, uadj)
.C("ARMARECURSION1", iStart = as.integer(iStart),
    iEnd = as.integer(iEnd), phi1 = as.double(phi1),
    theta1 = as.double(theta1), yzeroadj = as.double(yzeroadj),
    innov = as.double(innov), lny2adj = as.double(lny2adj),
    uadj = as.double(uadj), PACKAGE = "lgarch")

LGARCHSIM <- function (maxpq, nmaxpq, lnsigma2, phi, phisum, innov)
.C("LGARCHSIM", maxpq = as.integer(maxpq), nmaxpq = as.integer(nmaxpq),
  lnsigma2 = as.double(lnsigma2), phi = as.double(phi),
  phisum = as.double(phisum), innov = as.double(innov), PACKAGE = "lgarch")

VARMARECURSION1 <- function (iStart, n, m, mU, mY, mInnov,
    PHI, THETA, mYiszeroadj)
.C("VARMARECURSION1",  iStart = as.integer(iStart),
    n = as.integer(n), m = as.integer(m), mU = as.double(mU),
    mY = as.double(mY), mInnov = as.double(mInnov),
    PHI = as.double(PHI), THETA = as.double(THETA),
    mYiszeroadj = as.double(mYiszeroadj), PACKAGE = "lgarch")
