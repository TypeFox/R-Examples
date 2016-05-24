### Missing at random
MARnocov <- function(Y, D, Z) {
  R <- (!is.na(Y))*1
  ITT <- mean(Y[R == 1 & D == 1 & Z == 1])*mean(D[Z == 1]) +
    mean(Y[R == 1 & D == 0 & Z == 1])*(1-mean(D[Z == 1])) -
      mean(Y[R == 1 & D == 1 & Z == 0])*mean(D[Z == 0]) -
        mean(Y[R == 1 & D == 0 & Z == 0]) * (1-mean(D[Z == 0]))
  CACE <- ITT/(mean(D[Z==1])-mean(D[Z==0]))
  return(list(ITT = ITT, CACE = CACE))
}
