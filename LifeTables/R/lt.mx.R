lt.mx <- function (nmx, sex="female", age = c(0,1,seq(5,110,5)), nax=NULL){
  n <- c(diff(age), 999)
  if(is.null(nax)){
    nax <- 0.5*n
    if(n[2]==4){
      if(sex == "male"){
        if (nmx[1] >= 0.107) {
          nax[1] <- 0.33
          nax[2] <- 1.352
        } else {
          nax[1] <- 0.045 + 2.684 * nmx[1]
          nax[2] <- 1.651 - 2.816 * nmx[1]
        } # else
      } # if (sex == "male")
      if(sex == "female"){
        if (nmx[1] >= 0.107) {
          nax[1] <- 0.35
          nax[2] <- 1.361
        } else {
          nax[1] <- 0.053 + 2.8 * nmx[1]
          nax[2] <- 1.522 - 1.518 * nmx[1]
        } # else
      } #if (sex == "female")
    } #if(n[2]==4)
  } #if(is.null(nax))
  
  nqx <- (n * nmx)/(1 + (n - nax) * nmx)
  nqx <- c(nqx[-(length(nqx))], 1)
  for (i in 1:length(nqx)) {
    if (nqx[i] > 1) 
      nqx[i] <- 1
  }
  nage <- length(age)
 # nqx <- round(nqx, 4)
  npx <- 1 - nqx
  l0 = 1e+05
  lx <- round(cumprod(c(l0, npx)))
  ndx <- -diff(lx)
  lxpn <- lx[-1]
  nLx <- n * lxpn + ndx * nax
  Tx <- c(rev(cumsum(rev(nLx[-length(nLx)]))),0)
  lx <- lx[1:length(age)]
  ex <- Tx/lx
  lt <- cbind(Age = age, nax = c(round(nax[-length(nax)], 3),NA), nmx = round(nmx,4), nqx = round(nqx, 4), npx = round(npx, 4), ndx = ndx, 
              lx = lx, nLx = c(round(nLx[-length(nLx)]),NA), Tx = c(round(Tx[-length(Tx)]),NA), ex = c(round(ex[-length(ex)],2),NA))
  e0 <- lt[1, 10]
  lt.45q15 <- 1 - (lx[which(age==60)]/lx[which(age==15)])
  lt.5q0 <- 1 - (lx[which(age==5)]/lx[which(age==0)])
  lt.top.age <- min(which(nqx==1))
  lt <- lt[1:lt.top.age,]
  return(list(e0 = e0, lt.5q0 = lt.5q0, lt.45q15 = lt.45q15, 
              lt = lt))
}