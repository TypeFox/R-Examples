psi <-
function(a, f, mn, y, d0, wt) {
  0.5 * t(a) %*% (f - mn) - sum(wt * d0(f, y))
}
