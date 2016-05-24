EMGaussMax <-
function(Z, W, yw, ptype, mutype, sdtype, sd.fixed=0.05, p) {
   ZW <- W %*% Z
   mu    <- EMMax.mu(yw, ZW, mutype)
   sigma <- EMMax.sd(yw, ZW, mu, sdtype, sd.fixed)
   p     <- EMMax.p(Z, diag(W), p, ptype)

   list(mu=mu,sigma=sigma,p=p)
   #so mu is optimized first, these are used to optimize sigma
   #the new p are derived from the current p's of the samples
}
