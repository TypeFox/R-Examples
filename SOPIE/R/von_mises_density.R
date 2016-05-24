von_mises_density <-
function (x_grid, c, kappa, mu=pi)
{
     fv_x <- ((1/(2 * pi * besselI(kappa, 0))) * exp(kappa * cos(x_grid - pi)) - (exp(1 * kappa * cos(c - mu)) / 
             (2 * pi * besselI(kappa, 0)))) / (pvonmises(circular((2 * pi) - c), circular(pi), kappa) - pvonmises(circular(c), circular(pi), kappa) - 
             (((pi - c) / (pi * besselI(kappa, 0))) * exp(1 * kappa * cos(c - mu))))
     return (fv_x)
}
