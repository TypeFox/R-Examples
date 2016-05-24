EMMax.mu <-
function(yw, z, mutype) {
  if (mutype==0) apply(z,2,wmean,yw) # estimate free mu's
  else EMMax.mu.ratiodist(z,yw, mutype) # yw arsin-sqrt transformed, but distances between mu's on (0,1) scale according to ratio's based on copynumbers
}
