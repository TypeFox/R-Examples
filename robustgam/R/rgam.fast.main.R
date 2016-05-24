 rgam.fast.main <- function (Ry, RX, Rfamily, Rc, RB, RrS, Rexpect, Rw_fun, Rm_initial,
    Rbeta_old, Rcount_lim, Rw_count_lim, Rdisplay)
.Call("rgam_fast_main_cpp", Ry, RX, Rfamily, Rc, RB, RrS, Rexpect,
    Rw_fun, Rm_initial, Rbeta_old, Rcount_lim, Rw_count_lim,
    Rdisplay, PACKAGE = "robustgam")

