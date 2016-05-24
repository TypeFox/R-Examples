integration.constant <- function( density, w.lo, w.hi, ... ){
#
#   Return the scalar so that integral from 0 to w of underlying density
#   is 1.0
#
#   Input:
#   density = a function to compute integration constant for.
#       this function must be capable of evaluating values from 0 to w
#   w = upper limit of integral.
#   ... = additional arguments to density.  These vary by density function,
#       but generally are parameter values, series, expansion terms, etc.
#
#   Output:
#   a divisor scalar such that density / scalar integrates to 1.0. i.e.,
#   this output scalar is the integral of unscaled density from 0 to w.
#

density = match.fun(density)

seqx = seq(w.lo, w.hi, length=200)
seqy = density( dist=seqx, scale=FALSE, w.lo=w.lo, w.hi=w.hi,...)

#   Trapazoid rule
scaler= (seqx[2]-seqx[1]) * sum(seqy[-length(seqy)]+seqy[-1]) / 2

scaler

}
