fit.gam.sp1 <- function (Ry, RB, RrS, Rfamily)
.Call("fit_gam_sp_cpp", Ry, RB, RrS, Rfamily, PACKAGE = "robustgam")

