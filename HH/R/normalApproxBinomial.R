normalApproxBinomial <- function(p0= if (number.vars==1) .5 else 0,
                                 p1=NA, p2=NA,
                                 p.hat=if (number.vars==1) .75 else 0,
                                 n=1,
                                 xlim=if (number.vars==1) c(0,1) else c(-1,1),
                                 ylim=c(0, 5),
                                 type=c("hypothesis","confidence"),
                                 alpha.left=if (type=="hypothesis") 0 else .025,
                                 alpha.right=if (type=="hypothesis") .05 else .025,
                                 xlab=if (number.vars==1)
                                        "w = p = population proportion"
                                      else
                                        "w = p[1] - p[2] :: population proportions", ...,
                                 number.vars=if (!is.na(p1) && !is.na(p2)) 2 else 1) {

  if (number.vars == 2) stop("This function is not yet working for the two-sample case.", call.=FALSE)
  type <- match.arg(type)

  sigma.p0 <- sqrt(p0*(1-p0)/n)
  sigma.p1 <- sqrt(p1*(1-p1)/n)
  sigma.p2 <- sqrt(p2*(1-p2)/n)
  s.p.hat <- sqrt(p.hat*(1-p.hat)/n)
  sigma.p1.p2 <- sqrt(sigma.p1^2 + sigma.p2^2)
##  z.calc <- (p.hat-p0)/sigma.p0



  if (number.vars == 1) {
  if (type=="hypothesis")
    NormalAndTplot(mean0=p0, mean1=p1, sd=sigma.p0*sqrt(n), n=n, xbar=p.hat,
                   alpha.left=alpha.left, alpha.right=alpha.right,
                   xlim=xlim, ylim=ylim,
                   xlab=xlab, ..., distribution.name="binomial")
  else ## "confidence"
    NormalAndTplot(mean0=NA, mean1=NA, sd=s.p.hat*sqrt(n), n=n, xbar=p.hat,
                   alpha.left=alpha.left, alpha.right=alpha.right,
                   xlim=xlim, ylim=ylim,
                   xlab=xlab, type="confidence", ..., distribution.name="binomial", NTmethod="binomial")
}
  else {## number.vars == 2
  if (type=="hypothesis")
    NormalAndTplot(mean0=0, mean1=p1-p2, sd=sigma.p1.p2*sqrt(n), n=n, xbar=p.hat,
                   alpha.left=alpha.left, alpha.right=alpha.right,
                   xlim=xlim, ylim=ylim,
                   xlab=xlab, ...) ##, distribution.name="binomial", number.vars=2)
  else ## "confidence"
    NormalAndTplot(mean0=NA, mean1=NA, sd=sigma.p1.p2*sqrt(n), n=n, xbar=p1-p2,
                   alpha.left=alpha.left, alpha.right=alpha.right,
                   xlim=xlim, ylim=ylim,
                   xlab=xlab, type="confidence", ...) ##, distribution.name="binomial", NTmethod="binomial", number.vars=2)
}

}
