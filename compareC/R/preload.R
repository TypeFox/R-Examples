.onLoad <- function(lib, pkg)
{
    library.dynam("compareC", pkg, lib)
#   if (!is.loaded("TauXX")) dyn.load(paste("survivalC", .Platform$dynlib.ext, sep=""), PACKAGE = "compareC")
#   if (!is.loaded("TauXX")) print("functions in the shared library cannot be loaded!")
}

estC <- function(timeX,statusX,scoreY) {

Tau_XX <- function(timeX,statusX) {
  .C("TauXX",as.double(timeX),
     as.integer(statusX),
     as.integer(length(timeX)),
     as.double(1),PACKAGE = "compareC"
  )[[4]]
}

Tau_XY <- function(timeX,statusX,scoreY) {
  .C("TauXY",as.double(timeX),
     as.integer(statusX),
     as.double(scoreY),
     as.integer(length(timeX)),
     as.double(1),PACKAGE = "compareC"
  )[[5]]
}
  
return(as.numeric((Tau_XY(timeX,statusX,scoreY)/Tau_XX(timeX,statusX)+1)/2))
}


vardiffC <- function(timeX,statusX,scoreY,scoreZ) {

Tau_XX <- function(timeX,statusX) {
  .C("TauXX",as.double(timeX),
     as.integer(statusX),
     as.integer(length(timeX)),
     as.double(1),PACKAGE = "compareC"
  )[[4]]
}

Tau_XY <- function(timeX,statusX,scoreY) {
  .C("TauXY",as.double(timeX),
     as.integer(statusX),
     as.double(scoreY),
     as.integer(length(timeX)),
     as.double(1),PACKAGE = "compareC"
  )[[5]]
}
  
Var.Tau_XX <- function(timeX,statusX) {
  .C("VarTauXX",as.double(timeX),
     as.integer(statusX),
     as.integer(length(timeX)),
     as.double(1),PACKAGE = "compareC"
  )[[4]]
}

Var.Tau_XY <- function(timeX,statusX,scoreY) {
  .C("VarTauXY",as.double(timeX),
     as.integer(statusX),
     as.double(scoreY),
     as.integer(length(timeX)),
     as.double(1),PACKAGE = "compareC"
  )[[5]]
}

Cov.Tau_XXXY <- function(timeX,statusX,scoreY) {
  .C("CovTauXXXY",as.double(timeX),
     as.integer(statusX),
     as.double(scoreY),
     as.integer(length(timeX)),
     as.double(1),PACKAGE = "compareC"
  )[[5]]
}

Cov.Tau_XYXZ <- function(timeX,statusX,scoreY,scoreZ) {
  .C("CovTauXYXZ",as.double(timeX),
     as.integer(statusX),
     as.double(scoreY),
     as.double(scoreZ),
     as.integer(length(timeX)),
     as.double(1),PACKAGE = "compareC"
  )[[6]]
}
  
  
t11 = Tau_XX(timeX,statusX)
t12 = Tau_XY(timeX,statusX,scoreY)
t13 = Tau_XY(timeX,statusX,scoreZ)
var.t11 = Var.Tau_XX(timeX,statusX)
var.t12 = Var.Tau_XY(timeX,statusX,scoreY)
var.t13 = Var.Tau_XY(timeX,statusX,scoreZ)
cov.t1112 = Cov.Tau_XXXY(timeX,statusX,scoreY)
cov.t1113 = Cov.Tau_XXXY(timeX,statusX,scoreZ)
cov.t1213 = Cov.Tau_XYXZ(timeX,statusX,scoreY,scoreZ)

est.varCxy = as.numeric(1/4*(c(1/t11,-t12/t11^2)%*%matrix(c(var.t12,  cov.t1112,cov.t1112,var.t11),2)%*%c(1/t11,-t12/t11^2)))
est.varCxz = as.numeric(1/4*(c(1/t11,-t13/t11^2)%*%matrix(c(var.t13,  cov.t1113,cov.t1113,var.t11),2)%*%c(1/t11,-t13/t11^2)))
est.cov =    as.numeric(1/4*(c(1/t11,-t12/t11^2)%*%matrix(c(cov.t1213,cov.t1113,cov.t1112,var.t11),2)%*%c(1/t11,-t13/t11^2)))
est.vardiff_c = est.varCxy + est.varCxz - 2*est.cov

return(list(est.vardiff_c=est.vardiff_c,est.varCxy=est.varCxy,est.varCxz=est.varCxz,est.cov=est.cov))
}



compareC <- function(timeX,statusX,scoreY,scoreZ) {
est.c = c(estC(timeX,statusX,scoreY),estC(timeX,statusX,scoreZ))
names(est.c) <- c("Cxy","Cxz")
est.diff_c = est.c[1] - est.c[2]
names(est.diff_c) <- NULL

tmpout = vardiffC(timeX,statusX,scoreY,scoreZ)
zscore = est.diff_c/sqrt(tmpout$est.vardiff_c)
pval = 2*(1-pnorm(abs(zscore)))

return(list(est.c=est.c,est.diff_c=est.diff_c,est.vardiff_c=tmpout$est.vardiff_c, 
            est.varCxy=tmpout$est.varCxy,est.varCxz=tmpout$est.varCxz,est.cov=tmpout$est.cov,
            zscore=zscore,pval=pval))
}
