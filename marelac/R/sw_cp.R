## -----------------------------------------------------------------------------
##  Heat Capacity of Seawater
## -----------------------------------------------------------------------------

sw_cp <- function (S = 35 ,
                   t = 25,
                   p = P-1.013253,
                   P = 1.013253,
                   method = c("Gibbs", "UNESCO")) {
  if (any (S < 0))
    stop ("Salinity should be >= 0")
  method <- match.arg(method)

  if (method =="UNESCO") {
    P     <- p  # hydrostatic pressure is called "P" here...
    S3_2  <- S*sqrt(S)

    ## eqn 26 p.32 Fofonoff et al.
    ## specific heat for P=0
    Cpst0 <- 4217.4    - 3.720283 *t + 0.1412855*t^2 -2.654387e-3*t^3 +
             2.093236e-5*t^4 +
             (-7.64357  + 0.1072763*t -1.38385e-3*t^2)*S +
             (0.1770383 -4.07718e-3*t +  5.148e-5*t^2)*S*sqrt(S)

    ## eqn 28 p.33
    ## pressure and temperature terms for S=0

    del_Cp0t0 <- (-4.9592e-1 + 1.45747e-2*t - 3.13885e-4*t^2 +
                 2.0357e-6*t^3 + 1.7168e-8*t^4)*P    +
                 (2.4931e-4 - 1.08645e-5*t + 2.87533e-7*t^2 - 4.0027e-9*t^3 +
                 2.2956e-11*t^4)*P^2 +
                 (- 5.422e-8 +  2.6380e-9*t - 6.5637e-11*t^2 + 6.136e-13*t^3)*P^3

    ## eqn 29 p.34
    ## pressure and temperature terms for S>0

    del_Cpstp <-((4.9247e-3  -1.28315e-4*t + 9.802e-7*t^2 +
                2.5941e-8*t^3 -2.9179e-10*t^4)*S +
                (-1.2331e-4 -  1.517e-6*t + 3.122e-8*t^2)*S3_2)*P +
                ((-2.9558e-6 +1.17054e-7*t -2.3905e-9*t^2 + 1.8448e-11*t^3)*S +
                9.971e-8*S3_2)*P^2 +
                ((5.540e-10 - 1.7682e-11*t + 3.513e-13*t^2)*S +
                -1.4300e-12*t*S3_2)*P^3

    ## specific heat
    cp <- Cpst0 + del_Cp0t0 + del_Cpstp
  } else {
    cp <- -(t + 273.15) * sw_gibbs(S, t, p, dS = 0, dt = 2, dp = 0)
  }
  return(cp)
}
