grp = function(theta, bins=9, breaks=NULL, equal="count", type="means") {
  if (!is.null(dim(theta))) theta = theta[,1]
  if (is.null(breaks)) {
    if (equal=="count") {
      qu = seq(0, 1, 1/bins)
      breaks=quantile(theta, probs=qu, names=FALSE, na.rm=TRUE)
    } else {
      r = range(theta)
      st = (r[2]-r[1])/bins
      breaks = seq(r[1], r[2], by=st)
    }
  }
  grmemb = findInterval(theta, breaks, rightmost.closed=TRUE)
  counts = table(grmemb)
  ref = switch(type,
    mids = {
      dd = diff(breaks,1)
      md = mean(dd)/2
      c(breaks[1]-md, breaks[-length(breaks)]+dd/2, breaks[length(breaks)]+md)
    },
    meds = tapply(theta, grmemb, median),
    means= tapply(theta, grmemb, mean)
  )
  return(list(breaks=breaks,counts=counts,ref=ref,grmemb=grmemb))
}
# compute an item fit statistic with df and pval.
# Optionally plot the IRF with residuals shown


#' Test item fit
#' 
#' Returns a statistic of item fit together with its degrees of freedom and
#' p-value. Optionally produces a plot.
#' 
#' Given a long test, say 20 items or more, a large-test statistic of item fit
#' could be constructed by dividing examinees into groups of similar ability,
#' and comparing the observed proportion of correct answers in each group with
#' the expected proportion under the proposed model. Different statistics have
#' been proposed for this purpose.
#' 
#' The chi-squared statistic
#' \deqn{X^2=\sum_g(N_g\frac{(p_g-\pi_g)^2}{\pi_g(1-\pi_g)},} where \eqn{N_g}
#' is the number of examinees in group \eqn{g}, \eqn{p_g=r_g/N_g}, \eqn{r_g} is
#' the number of correct responses to the item in group \eqn{g}, and
#' \eqn{\pi_g} is the IRF of the proposed model for the median ability in group
#' \eqn{g}, is attributed by Embretson & Reise to R. D. Bock, although the
#' article they cite does not actually mention it. The statistic is the sum of
#' the squares of quantities that are often called "Pearson residuals" in the
#' literature on categorical data analysis.
#' 
#' BILOG uses the likelihood-ratio statistic
#' \deqn{X^2=2\sum_g\left[r_g\log\frac{p_g}{\pi_g} +
#' (N_g-r_g)\log\frac{(1-p_g)}{(1-\pi_g)}\right],} where \eqn{\pi_g} is now the
#' IRF for the mean ability in group \eqn{g}, and all other symbols are as
#' above.
#' 
#' Both statistics are assumed to follow the chi-squared distribution with
#' degrees of freedom equal to the number of groups minus the number of
#' parameters of the model (eg 2 in the case of the 2PL model). The first
#' statistic is obtained in \code{itf} with \code{stat="chi"}, and the second
#' with \code{stat="lr"} (or not specifying \code{stat} at all).
#' 
#' In the real world we can only work with estimates of ability, not with
#' ability itself. \code{irtoys} allows use of any suitable ability measure 
#' via the argument \code{theta}. If \code{theta} is not specified, \code{itf} 
#' will compute EAP estimates of ability, group them in 9 groups having 
#' approximately the same number of cases, and use the means of the ability
#' eatimates in each group. This is the approximate behaviour of BILOG. 
#' 
#' If the test has less than 20 items, \code{itf} will issue a warning.
#' For tests of 10 items or less, BILOG has a special statistic of fit, which
#' can be found in the BILOG output. Also of interest is the fit in 2- and
#' 3-way marginal tables in package \code{ltm}.
#' 
#' @param resp A matrix of responses: persons as rows, items as columns,
#' entries are either 0 or 1, no missing data
#' @param ip Item parameters: a matrix with one row per item, and three
#' columns: [,1] item discrimination \eqn{a}, [,2] item difficulty \eqn{b}, and
#' [,3] asymptote \eqn{c}.
#' @param item A single number pointing to the item (column of \code{resp}, row
#' of \code{ip}), for which fit is to be tested
#' @param stat The statistic to be computed, either \code{"chi"} or
#' \code{"lr"}. Default is \code{"lr"}. See details below.
#' @param theta A vector containing some viable estimate of the latent variable
#' for the same persons whose responses are given in \code{resp}. If not given
#' (and \code{group} is also missing), EAP estimates will be computed from
#' \code{resp} and \code{ip}.
#' @param standardize Standardize the distribution of ability estimates?
#' @param mu Mean of the standardized distribution of ability estimates
#' @param sigma Standard deviation of the standardized distribution of ability
#' estimates
#' @param bins Desired number of bins (default is 9)
#' @param breaks A vector of cutpoints. Overrides \code{bins} if present.
#' @param equal Either \code{"width"} for bins of equal width, or
#' \code{"count"} for bins with roughly counts of observations. Default is
#' \code{"quant"}
#' @param type The points at which \code{itf} will evaluate the IRF. One of
#' \code{"mids"} (the mid-point of each bin), \code{"meds"} (the median of the
#' values in the bin), or \code{"means"} (the mean of the values in the bin).
#' Default is \code{"means"}.
#' @param do.plot Whether to do a plot
#' @param main The title of the plot if one is desired
#' @return A vector of three numbers: \item{Statistic}{The value of the statistic of item fit}
#' \item{DF}{The degrees of freedom} \item{P-value}{The p-value}
#' @author Ivailo Partchev
#' @seealso \code{\link{eap}}, \code{\link{qrs}}
#' @references S. E. Embretson and S. P. Reise (2000), Item Response Theory for
#' Psychologists, Lawrence Erlbaum Associates, Mahwah, NJ
#' 
#' M. F. Zimowski, E. Muraki, R. J. Mislevy and R. D. Bock (1996), BILOG--MG.
#' Multiple-Group IRT Analysis and Test Maintenance for Binary Items, SSI
#' Scientific Software International, Chicago, IL
#' @keywords models
#' @export
#' @examples
#' 
#' fit   <- itf(resp=Scored, ip=Scored2pl$est, item=7)
#' 
itf = function(resp, ip, item, stat = "lr", theta,
  standardize=TRUE, mu=0, sigma=1, 
  bins = 9, breaks = NULL, equal = "count", type = "means",
  do.plot=TRUE, main="Item fit") {
  if (!(item %in% 1:nrow(ip))) stop("bad item number")
  if (nrow(ip)<20) warning("item fit statistic computed for a test of less than 20 items")
  if (missing(theta)) theta = eap(resp, ip, normal.qu())
  if (!is.null(dim(theta))) theta = theta[, 1]
  pa = ip[item,]
  iresp = resp[,item]
  oo = !is.na(iresp)
  iresp = iresp[oo]
  theta = theta[oo]
  if (standardize) theta = scale(theta)*sigma + mu
  groups = grp(theta, bins=bins, breaks=breaks, equal=equal, type=type)
  ep = as.vector(irf(ip=pa, x=groups$ref)$f)
  nn = groups$counts
  gr = groups$grmemb * iresp
  rg = tabulate(gr[gr>0], nbins=length(nn))
  op = rg / nn
  npar = 3 - (all(ip[,3]==0)) - (var(ip[,1])==0)
  if (stat=="chi") tst = sum((op-ep)^2/(ep*(1-ep))*nn)
  else {
    cmp = rg*log(op/ep)+(nn-rg)*log((1-op)/(1-ep))
    tst = 2*sum(cmp[is.finite(cmp)])
  }
  dfr = sum(nn>0) - npar
  pval=pchisq(tst,dfr,lower.tail=FALSE)
  if (do.plot) {
    ifit = paste("Q = ",round(tst,2),"    d.f. = ",dfr,"    P = ",round(pval,4),sep="")
    plot(c(-4,4),c(0,1),main="",xlab="Ability",ylab="Proportion right",type="n")
    title(main=main,mtext(ifit,line=0.5,cex=0.7))
    points(groups$ref, op)
    plot(irf(pa),add=TRUE)
  }
  result = c(tst,dfr,pval)
  names(result) = c("Statistic","DF","P-value")
  result
}


#' The Z3 appropriateness index
#' 
#' Computes the Z3 appropriateness index, a measure of person fit in IRT models
#' 
#' 
#' @param resp A matrix of responses: persons as rows, items as columns,
#' entries are either 0 or 1, no missing data
#' @param ip Item parameters: a matrix with one row per item, and three
#' columns: [,1] item discrimination \eqn{a}, [,2] item difficulty \eqn{b}, and
#' [,3] asymptote \eqn{c}.
#' @param theta A measure of ability, typically produced with \code{mlebme}, 
#' \code{wle} etc. If missing, ML estimates will be computed automatically.
#' @return A vector of length equal to the number of rows in \code{resp},
#' containing the appropriateness indices
#' @author Ivailo Partchev
#' @references Drasgow, F., Levine, M. V., & Williams, E. A. (1985).
#' Appropriateness measurement with polychotomous item response models and
#' standardized indices. British Journal of Mathematical and Statistical
#' Psychology, 38, 67--80
#' @keywords models
#' @export
#' @examples
#' 
#' api(Scored, Scored2pl$est)
#' 
api = function(resp, ip, theta){
  if(missing(theta)) theta = mlebme(resp, ip)
  if (!is.null(dim(theta))) theta = theta[, 1]
  p = irf(ip, theta)$f
	p = pmax(p, .00001); pr = pmin(p, .99999)
  if(is.null(dim(p))) p=matrix(p,ncol=length(p))
  q = 1 - p
  lp= log(p)
  lq= log(q)
  l = apply(lp*resp + lq*(1-resp), 1, sum)
  m = apply(p*lp+q*lq, 1, sum)
  v = apply(p*q*(log(p/q))^2, 1, sum)
  (l-m)/sqrt(v)
}
