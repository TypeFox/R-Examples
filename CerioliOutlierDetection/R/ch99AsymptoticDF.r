ch99AsymptoticDF <- 
#
# computes the asymptotic degrees of freedom and
# consistency constants as described in
# Croux and Haesbroeck (1999) and 
# Hardin and Rocke (2005)
#
# Christopher G. Green
# 2011
#
function( n.obs, p.dim, mcd.alpha=max.bdp.mcd.alpha(n.obs,p.dim) )
#
{

  # 2011-08-12
  # mcd.alpha = 1 case: estimator is the sample covariance
  # matrix, so c.alpha = 1 and m = n.obs - 1
  if ( mcd.alpha == 1 ) {
    c.alpha          <- 1.0
    m.hat.asy        <- n.obs - 1.0

  } else {

    # constants from Croux and Haesbroeck 1999 paper
    # 1 - alpha is equal to h/n, where h is the number
    # of observations used to compute the MCD and n
    # is the total number of observations

    # 2011-08-12 for n small relative to p this will not quite
    # be 1-alpha. just use mcd.alpha
    #one.minus.alpha  <- robustbase::h.alpha.n(mcd.alpha,n.obs,p.dim)/n.obs
    one.minus.alpha  <- mcd.alpha

    # G(t) = P(Z'Z <= t | F_{0,I})
    # 1-alpha quantile of the G distribution
    # For MVN sqrt(Z'Z) is chi square(p)
    # q.alpha here is sqrt(q.alpha) in CH99 paper
    # and works out to be 1 - alpha upper quantile of
    # chi square(p)
    q.alpha          <- qchisq(one.minus.alpha, df=p.dim)
    # probability of square of z[1] being less than q.alpha
    p.alpha          <- pchisq(q.alpha, df=p.dim+2)
    # consistency factor
    #c.alpha          <- one.minus.alpha / p.alpha
    #c.alpha          <- robustbase:::MCDcons(p.dim, mcd.alpha) 
	# 2014-07-28 use exported version now
    c.alpha          <- robustbase::.MCDcons(p.dim, mcd.alpha) 
    #robustbase:::MCDcnp2(v, n, mcd.alpha)
    # constants from theorem 1 of CH99 for the case of
    # a normal distribution
    c2               <- -0.5*p.alpha
    c3               <- ifelse(p.dim >= 2, -0.5*pchisq(q.alpha, df=p.dim+4), 0)
    #c4              <- 3*c3
    #b1              <- c.alpha * (c3 - c4 )/one.minus.alpha
    #b1               <- c.alpha * (-2. * c3) /one.minus.alpha
    # 2011-08-12
    # c.alpha/one.minus.alpha = 1/p.alpha
    b1               <- (-2. * c3)/p.alpha 
    # 2011-08-12
    # for one.minus.alpha near 1, part of this expression
    # can be numerically unstable (q.alpha goes to Inf, c2 + 0.5*one.minus.alpha 
    # goes to 0)
    #b2               <- 0.5 + (c.alpha/one.minus.alpha)*(c3 - 
    #            (q.alpha/p.dim)*(c2 + 0.5*one.minus.alpha))
    y1               <- q.alpha * (-p.alpha + one.minus.alpha)
    b2               <- 0.5 + (c3 - 0.5*y1/p.dim)/p.alpha
    #cat("y1 = ",y1,"b2 = ",b2,"\n")

    z                <- b1 - p.dim*b2
    # 2011-06-26 cgg fixed slight bug here in placement of parentheses
    # maybe??? think R's order of operations was doing the right thing
    # before. Better to be clear here.
    #v1               <- one.minus.alpha * b1 * b1 * 
    #              ( mcd.alpha * ( (c.alpha*q.alpha)/p.dim - 1. )^2.  - 1.) 
    # 2011-08-12 mcd.alpha * q.alpha is another unstable quantity
    # should also be 1 - mcd.alpha here, not mcd.alpha (probably more important
    # than the above bug!)
    y2               <- (1. - mcd.alpha) * (( c.alpha * q.alpha )/p.dim - 1.)^2.
    v1               <- one.minus.alpha * b1 * b1 * ( y2 - 1.) 
    v1               <- v1 - 2.*c3*c.alpha*c.alpha*(3.*z*z + (p.dim + 2.)*b2*(b1 + z))
    v2               <- n.obs*c.alpha*c.alpha*(b1*z*one.minus.alpha)^2.
    v                <- v1/v2
    m.hat.asy        <- 2. / (c.alpha*c.alpha*v)

  }

  list(c.alpha=c.alpha, m.hat.asy=m.hat.asy)
}
