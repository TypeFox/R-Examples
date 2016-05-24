 `resFun` <-
function(par, volFreq, t, sigma_t, V) {
   ff <- fitFun(par, t, sigma_t, V) 
   volFreq - ff 
 }
`resFunChiSq` <-
function(par, volFreq, t, sigma_t, V, skel) {
  ff <- fitFun(relist(par,skel), t=t, sigma_t = sigma_t, V=V)
  gg <- log(volFreq / ff)
  x2 <- which(is.finite(gg))
  ## multinomial likelihood chi-square
  2*sum(volFreq[x2]*gg[x2])
}

