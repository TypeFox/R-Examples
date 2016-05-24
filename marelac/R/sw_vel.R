## -----------------------------------------------------------------------------
## Velocity of Sound in Seawater
## -----------------------------------------------------------------------------

sw_svel <- function (S = 35, t = 25, p = P-1.013253, P = 1.013253,
                   method = c("Gibbs","UNESCO")) {
  if (any (S < 0))
    stop ("Salinity should be >= 0")

   method <- match.arg(method)
  if (method=="UNESCO") {
    P  <- p  # P is used in code as synonym for hydroP
    T  <- t
    P2 <- P*P
    P3 <- P2*P
    Cw <- 1402.388 + (5.03711   + (-5.80852e-2 + (3.3420e-4 +
                                (-1.47800e-6 + 3.1464e-9*T)*T)*T)*T)*T   +
       + (0.153563 + (6.8982e-4 + (-8.1788e-6  +
                                (1.3621e-7 -6.1185e-10*T)*T)*T)*T)*P     +
       + (3.1260e-5 +(-1.7107e-6+ (2.5974e-8 +
                                (-2.5335e-10 + 1.0405e-12*T)*T)*T)*T)*P2 +
       + (-9.7729e-9+(3.8504e-10 -2.3643e-12*T)*T)*P3

    A <- 1.389      + (-1.262e-2  + (7.164e-5    +
                                (2.006e-6  -3.21e-8*T)*T)*T)*T           +
      + (9.4742e-5  + (-1.2580e-5 + (-6.4885e-8  +
                                 (1.0507e-8 -2.0122e-10*T)*T)*T)*T)*P    +
      + (-3.9064e-7 + (9.1041e-9  + (-1.6002e-10 + 7.988e-12*T)*T)*T)*P2 +
      + (1.100e-10  + (6.649e-12  -3.389e-13*T)*T)*P3

    B <- -1.922e-2 -4.42e-5*T + (7.3637e-5 + 1.7945e-7*T)*P

    D <- 1.727e-3 + -7.9836e-6*P
    svel <- Cw + A*S + B*S**1.5 + D*S**2
  } else {

    g_tt <- sw_gibbs(S = S, t = t, p = p, dS = 0, dt = 2, dp = 0)
    g_2  <- sw_gibbs(S = S, t = t, p = p, dS = 0, dt = 0, dp = 1)
    g_tp <- sw_gibbs(S = S, t = t, p = p, dS = 0, dt = 1, dp = 1)
    svel <- sw_gibbs(S = S, t = t, p = p, dS = 0, dt = 0, dp = 1) *
      sqrt(g_tt/(g_tp*g_tp - g_tt *
                 sw_gibbs(S = S, t = t, p = p, dS = 0, dt = 0, dp = 2)))
  }
  return (svel)
}
