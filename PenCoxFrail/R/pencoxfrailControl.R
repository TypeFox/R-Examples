### control
pencoxfrailControl<-function(start = NULL, q_start = NULL, conv.eps = 1e-4, standardize = FALSE, center = FALSE,
                          smooth=list(nbasis = 6, penal = 0.1), ridge.pen = 1e-4,
                          print.iter = FALSE, max.iter = 100, c.app = 1e-6, zeta = 0.5, 
                          exact = 1e-2, xr = NULL,...)
{                       
  list(start = start, q_start = q_start, conv.eps = conv.eps, standardize = standardize, center = center,
       smooth = smooth, print.iter = print.iter, ridge.pen = ridge.pen,
       max.iter = max.iter, c.app =c.app, zeta = zeta, exact = exact, xr = xr,...)
}
