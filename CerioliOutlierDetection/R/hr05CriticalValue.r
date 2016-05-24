hr05CriticalValue <-
# Hardin and Rocke 2005, page 938
# cut off value for MCD distances for
# observations not included in the MCD
# subset
#
# it is assumed that the MCD covariance
# estimate was adjusted by the consistency
# factor c, so it is not included here.
#
# if it is needed the Croux Haesbroeck code
# provides 1/c as c.alpha. Also see MCDcons
# in the robustbase package.
# 
# AUTHOR: Christopher G. Green
# DATE: 2011
#
# ARGUMENTS
#
# em - the degrees of freedom
# p.dim - dimension
# signif.alpha - significance level for the test
#
function(em,p.dim,signif.alpha) {
  mmm    <- em - p.dim + 1
  if ( mmm > 0 ) {
    # cat("mmm: ",mmm,"em: ",em,"v: ",v,"\n")
    ((p.dim * em)/mmm) * 
      qf( 1 - signif.alpha, df1 = p.dim, df2 = mmm )
  } else {
    warning("In hr05CriticalValue: computed degrees of freedom was non-positive.")
    NA
  }
}
