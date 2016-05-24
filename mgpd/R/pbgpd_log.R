pbgpd_log <-
function
(x, y,
mar1 = c(0, 1, 0.1), mar2 = c(0, 1, 0.1),
dep = 2, ...)
{
    mu1     = expression((x^(-alpha) + y^(-alpha))^(1/alpha))
    mu      = function(x, y) eval({x<-x; y<-y; mu1})
    param   = as.numeric(c(mar1, mar2, dep))
    mux     = param[1]; muy   = param[4]
    sigx    = param[2]; sigy  = param[5]
    gamx    = param[3]; gamy  = param[6]
    alpha   = param[7]
    Hxy     = NULL

    error = FALSE
    if (sigx < 0 | sigy < 0 | alpha < 1) error = TRUE
    if (!error)
    { Hxy   = NA
      tx    = (1 + gamx * (x - mux)/sigx)^(1/gamx)
      ty    = (1 + gamy * (y - muy)/sigy)^(1/gamy)
      tx0   = (1 + gamx * (-mux)/sigx)^(1/gamx)
      ty0   = (1 + gamy * (-muy)/sigy)^(1/gamy)
      c0    = -mu(tx0, ty0)
      Hxy   = 1/c0 * (mu(tx, ty) - mu(pmin(tx,rep(tx0,length(tx))), pmin(ty,rep(ty0,length(ty)))))
}     else stop("invalid parameter(s)")
      Hxy
}
