#' Item response function
#' 
#' Returns the item response function of the 3PL (1PL, 2PL) model, the i.e. the
#' probabilities defined by
#' \deqn{P(U_{ij}=1|\theta_i,a_j,b_j,c_j)=c_j+(1-c_j)\frac{\displaystyle\exp(a_j(\theta_i-b_j))}{1+\displaystyle\exp(a_j(\theta_i-b_j))}}
#' where \eqn{U_{ij}} is a binary response given by person \eqn{i} to item
#' \eqn{j}, \eqn{\theta_i} is the value of the latent variable ("ability") for
#' person \eqn{i}, \eqn{a_j} is the discrimination parameter for item \eqn{j},
#' \eqn{b_j} is the difficulty parameter for item \eqn{j}, \eqn{c_j} is the
#' asymptote for item \eqn{j}. Some authors call the IRF "the item characteristic curve".
#' 
#' In the 2PL model (\code{model="2PL"}), all asymptotes \eqn{c_j} are 0. In
#' the 1PL model (\code{model="1PL"}), all asymptotes \eqn{c_j} are 0 and the
#' discriminations \eqn{a_j} are equal for all items (and sometimes to 1).
#' 
#' A common use of this function would be to obtain a plot of the IRF.
#' 
#' 
#' @param ip Item parameters: a matrix with one row per item, and three
#' columns: [,1] item discrimination \eqn{a}, [,2] item difficulty \eqn{b}, and
#' [,3] asymptote \eqn{c}.
#' @param x The values of the latent variable (\eqn{\theta} in the equation
#' above), at which the IRF will be evaluated. If not given, 99 values spaced
#' evenly between -4 and +4 will be used, handy for plotting.
#' @return A list of: \item{x}{A copy of the argument \code{x}} \item{f}{A
#' matrix containing the IRF values: persons (values of (\code{x}) as rows and
#' items as columns}
#' @author Ivailo Partchev
#' @seealso \code{\link{plot.irf}}
#' @keywords models
#' @export
#' @examples
#' 
#' plot(irf(Scored2pl$est[1,]))
#' 
irf = function(ip,x=NULL) {
 if (is.null(x)) 
    x = seq(-4, 4, length = 101)
  if (is.null(dim(ip))) 
    dim(ip) = c(1, 3)
  f = sweep(outer(x, ip[,2], "-"), 2, ip[,1], "*")
  f = 1 / (1 + exp(-f))
  if (any(ip[,3]!=0)) 
    f = sweep(sweep(f, 2, 1-ip[,3], "*"), 2, ip[,3], "+")
  r = list(x = x, f = f)
  class(r) = "irf"
  return(r)
}



#' A plot method for item response functions
#' 
#' Useful for plotting item response functions. The \code{x} argument of
#' \code{irf} should better be left out unless something special is required.
#' 
#' 
#' @param x An object produced by function \code{irf}
#' @param add When \code{add=T}, the IRF is added to a plot, otherwise a new
#' plot is started. Default is F.
#' @param main The main title of the plot, given that \code{add=F}.
#' @param co The colour of the IRF curve. Default is 1 for black. Use
#' \code{co=NA} to plot each IRF in a different colour.
#' @param label When \code{label=T}, individual curves will be labeled with the
#' item number.
#' @param ... Any additional plotting parameters
#' @author Ivailo Partchev
#' @seealso \code{\link{irf}}
#' @keywords models
#' @method plot irf
#' @S3method plot irf
#' @examples
#' 
#' # plot IRF for all items in red, label with item number
#' plot(irf(Scored2pl$est), co="red", label=TRUE)
#' # plot IRF for items 2, 3, and 7 in different colours
#' plot(irf(Scored2pl$est[c(2,3,7),]), co=NA)
#' 
plot.irf = function(x,  
    add=FALSE, main="Item response function", co=1, label=FALSE, ...) {
  if (!add) plot(c(min(x$x), max(x$x)), c(0,1), ty="n", xlab="Ability",
  ylab="Probability of a correct response", main=main)
  invisible(lapply(1:ncol(x$f), function(i) {
    if (is.na(co)) co = i
    lines(x$x, x$f[,i], lw=2, co=co)
  }))
  if (label) invisible(lapply(1:ncol(x$f), function(i) {
    lx = sample(1:length(x$x), 1)
    points(x$x[lx], x$f[lx,i], co="white", cex=1.6, pch=19)
    text(x$x[lx], x$f[lx,i],i, co=1, cex=.6)
  })) 
}

#' Test response function
#' 
#' Returns the test response function (TRF) of the 3PL (1PL, 2PL) model. The
#' TRF is the sum of the item response functions (IRF) in a test, and
#' represents the expected score as a function of the latent variable
#' \eqn{\theta}.
#' 
#' A common use of this function would be to obtain a plot of the TRF.
#' 
#' 
#' @param ip Item parameters: a matrix with one row per item, and three
#' columns: [,1] item discrimination \eqn{a}, [,2] item difficulty \eqn{b}, and
#' [,3] asymptote \eqn{c}.
#' @param x The values of the latent variable (\eqn{\theta} in the equation
#' above), at which the IRF will be evaluated. If not given, 99 values spaced
#' evenly between -4 and +4 will be used, handy for plotting.
#' @return A list of: \item{x}{A copy of the argument \code{x}} \item{f}{A
#' vector containing the TRF values}
#' @author Ivailo Partchev
#' @seealso \code{\link{plot.trf}}, \code{\link{irf}}
#' @keywords models
#' @export
#' @examples
#' 
#' plot(trf(Scored2pl$est))
#' 
trf = function(ip,x=NULL) {
  i = irf(ip,x)
  if (is.null(dim(i$f))) dim(i$f) = c(length(i$x),length(i$f))
  f = apply(i$f,1,sum)
  r = list(x=i$x, f=f, ni=ncol(i$f))
  class(r) = "trf"
  return(r)
}



#' A plot method for test response functions
#' 
#' Useful for plotting test response functions. The \code{x} argument of
#' \code{trf} should better be left out unless something special is required.
#' 
#' 
#' @param x An object produced by function \code{trf}
#' @param add When \code{add=T}, the IRF is added to a plot, otherwise a new
#' plot is started. Default is F.
#' @param main The main title of the plot, given that \code{add=F}.
#' @param co The colour of the TRF curve. Default is 1 for black. Use
#' \code{co=NA} to plot each TRF in a different colour.
#' @param ... Any additional plotting parameters
#' @author Ivailo Partchev
#' @seealso \code{\link{trf}}
#' @keywords models
#' @method plot trf
#' @S3method plot trf
#' @examples
#' 
#' plot(trf(Scored2pl$est))
#' 
plot.trf = function(x, 
    add=FALSE, main="Test response function", co=1, ...) {
  if (is.na(co)) co = 1
  if (!add) plot(c(min(x$x), max(x$x)), c(0,x$ni), ty="n", xlab="Ability",
  ylab="Expected score", main=main)
  lines(x$x, x$f, lw=2, co=co)
}

#' Item information function
#' 
#' The item information function (IIF) for the 3PL model can be computed as
#' \deqn{I(\theta) =
#' a^2\frac{Q(\theta)}{P(\theta)}\left[\frac{P(\theta)-c}{1-c}\right]^2,} where
#' \eqn{\theta} is the value of the latent variable for a person, \eqn{a} is
#' the discrimination parameter for the item, \eqn{P} is the IRF for the person
#' and item, and \eqn{Q=1-P}. For the 1PL and 2PL models, the expression
#' reduces to \eqn{a^2PQ}.
#' 
#' A common use of this function would be to obtain a plot of the IIF.
#' 
#' 
#' @param ip Item parameters: a matrix with one row per item, and three
#' columns: [,1] item discrimination \eqn{a}, [,2] item difficulty \eqn{b}, and
#' [,3] asymptote \eqn{c}.
#' @param x The values of the latent variable (\eqn{\theta} in the equation
#' above), at which the IIF will be evaluated. If not given, 99 values spaced
#' evenly between -4 and +4 will be used, handy for plotting.
#' @return A list of: \item{x}{A copy of the argument \code{x}} \item{f}{A
#' matrix containing the IIF values: persons (values of (\code{x}) as rows and
#' items as columns}
#' @author Ivailo Partchev
#' @seealso \code{\link{plot.iif}}, \code{\link{irf}}
#' @keywords models
#' @export
#' @examples
#' 
#' plot(iif(Scored2pl$est[1:3,]))
#' 
iif = function(ip, x=NULL) {
  if (is.null(x)) 
    x = seq(-4, 4, length = 101)
  if (is.null(dim(ip))) 
    dim(ip) = c(1, 3)
  p = irf(ip, x)$f  
  if (any(ip[, 3] != 0)) {
    ip[,3] = 0
    q = irf(ip, x)$f
    f = q^2*(1-p)/p
  } else 
    f = p*(1-p)
  f = sweep(f, 2, ip[,1]^2, "*")  
  r = list(x = x, f = f)
  class(r) = "iif"
  return(r)
}



#' A plot method for item information functions
#' 
#' Useful for plotting item information functions. The \code{x} argument of
#' \code{iif} should better be left out unless something special is required.
#' 
#' 
#' @param x An object produced by function \code{iif}
#' @param add When \code{add=T}, the IIF is added to a plot, otherwise a new
#' plot is started. Default is F.
#' @param main The main title of the plot, given that \code{add=F}.
#' @param co The colour of the IIF curve. Default is 1 for black. Use
#' \code{co=NA} to plot each IIF in a different colour.
#' @param label When \code{label=T}, individual curves will be labeled with the
#' item number.
#' @param ... Any additional plotting parameters
#' @author Ivailo Partchev
#' @seealso \code{\link{iif}}
#' @keywords models
#' @method plot iif
#' @S3method plot iif
#' @examples
#' 
#' # plot IIF for all items in red, label with item number
#' plot(iif(Scored2pl$est), co="red", label=TRUE)
#' # plot IIF for items 2, 3, and 7 in different colours
#' plot(iif(Scored2pl$est[c(2,3,7),]), co=NA)
#' 
plot.iif = function(x,  
  add=FALSE, main="Item information function", co=1, label=FALSE, ...) {
  if (!add) plot(c(min(x$x), max(x$x)), c(0,max(x$f)), ty="n", xlab="Ability",
  ylab="Item information",main=main)
  invisible(lapply(1:ncol(x$f), function(i) {
    if (is.na(co)) co = i
    lines(x$x, x$f[,i], lw=2, co=co)
    }))
    if (label) invisible(lapply(1:ncol(x$f), function(i) {
      lx = sample(1:length(x$x),1)
      points(x$x[lx], x$f[lx,i], co="white", cex=1.6, pch=19)
      text(x$x[lx], x$f[lx,i], i, co=1, cex=.6)
    })) 
}


#' Test information function
#' 
#' Returns the test information function (TIF) of the 3PL (1PL, 2PL) model. The
#' TIF is the sum of the item information functions (IIF) in a test, and
#' indicates the precision of measurement that can be achieved with the test at
#' any value of the latent variable, bein inversely related to measurement
#' variance.
#' 
#' A common use of this function would be to obtain a plot of the TIF.
#' 
#' 
#' @param ip Item parameters: a matrix with one row per item, and three
#' columns: [,1] item discrimination \eqn{a}, [,2] item difficulty \eqn{b}, and
#' [,3] asymptote \eqn{c}.
#' @param x The values of the latent variable (\eqn{\theta} in the equation
#' above), at which the TIF will be evaluated. If not given, 99 values spaced
#' evenly between -4 and +4 will be used, handy for plotting.
#' @return A list of: \item{x}{A copy of the argument \code{x}} \item{f}{A
#' vector containing the TIF values}
#' @author Ivailo Partchev
#' @seealso \code{\link{plot.tif}}, \code{\link{iif}}
#' @keywords models
#' @export
#' @examples
#' 
#' plot(trf(Scored2pl$est))
#' 
tif = function(ip, x=NULL) {
  i = iif(ip, x)
  if (is.null(dim(i$f))) dim(i$f) = c(length(i$x),length(i$f))
  f = apply(i$f, 1, sum)
  r = list(x=i$x, f=f, ni=ncol(i$f))
  class(r) = "tif"
  return(r)
}



#' A plot method for test information functions
#' 
#' Useful for plotting test information functions. The \code{x} argument of
#' \code{tif} should better be left out unless something special is required.
#' 
#' 
#' @param x An object produced by function \code{tif}
#' @param add When \code{add=T}, the TIF is added to a plot, otherwise a new
#' plot is started. Default is F.
#' @param main The main title of the plot, given that \code{add=F}.
#' @param co The colour of the TIF curve. Default is 1 for black. Use
#' \code{co=NA} to plot each TIF in a different colour.
#' @param ... Any additional plotting parameters
#' @author Ivailo Partchev
#' @seealso \code{\link{tif}}
#' @keywords models
#' @method plot tif
#' @S3method plot tif
#' @examples
#' 
#' plot(tif(Scored2pl$est))
#' 
plot.tif = function(x, add=FALSE, main="Test information function", co=1, ...) {
  if (is.na(co)) co = 1
  if (!add) plot(c(min(x$x), max(x$x)), c(0,max(x$f)), ty="n", xlab="Ability",
  ylab="Information", main=main)
  lines(x$x, x$f, lw=2, co=co)
}


#' True scores with standard errors
#' 
#' Computes the IRT true scores (test response function at the estimated
#' ability) and an estimate of their standard error via the delta theorem,
#' treating item parameters as known).
#' 
#' 
#' @param ip Item parameters: a matrix with one row per item, and three
#' columns: [,1] item discrimination \eqn{a}, [,2] item difficulty \eqn{b}, and
#' [,3] asymptote \eqn{c}.
#' @param theta An object containing ability estimates, as output by function
#' \code{mlebme} or \code{eap}
#' @return A matrix with the true scores in column 1, and their standard errors
#' of measurement (SEM) in column 2
#' @author Ivailo Partchev
#' @seealso \code{\link{mlebme}}, \code{\link{eap}}, \code{\link{trf}}
#' @keywords models
#' @export
#' @examples
#' 
#' th <- mlebme(resp=Scored, ip=Scored2pl$est)
#' tsc(Scored2pl$est, th)
#' 
tsc = function(ip, theta){
   p = irf(ip, theta[,1])$f
   if (is.null(dim(ip))) dim(ip) = c(1,3)
   if (is.null(dim(p)))  p = matrix(p, ncol=1)
   sc = apply(p, 1, sum)
   aq = sweep(1-p, 2, ip[,1], "*")
   if(any(ip[,3]!=0)) {
      p = sweep(p, 2, ip[,3],  "-")
      p = sweep(p, 2, 1-ip[,3],"/")
   }
   jb = apply(aq*p, 1, sum)
   se = jb*theta[,2]
   cbind(sc,se)
}

#' Plot observed and predicted scores against ability
#' 
#' Produces a plot of IRT true scores (test response function at the estimated
#' ability) with a confidence band (plus/minus standard error). The observed
#' sum scores are shown in red.
#' 
#' 
#' @param resp A matrix of binary responses to a test, with persons as rows and
#' items as columns.
#' @param ip Item parameters: a matrix with one row per item, and three
#' columns: [,1] item discrimination \eqn{a}, [,2] item difficulty \eqn{b}, and
#' [,3] asymptote \eqn{c}.
#' @param theta An object containing ability estimates, as output by function
#' \code{mlebme} or \code{eap}.  If \code{NULL}, MLE will be estimated from
#' \code{resp} and \code{ip}.
#' @return None
#' @author Ivailo Partchev
#' @seealso \code{\link{mlebme}}, \code{\link{eap}}, \code{\link{tsc}},
#' \code{\link{trf}}
#' @keywords models
#' @export
#' @examples
#' 
#' scp(Scored, Scored2pl$est)
#' 
scp = function(resp, ip, theta=NULL) {
  if (is.null(theta)) theta=mlebme(resp,ip)
	or = order(theta[,1])
	theta = theta[or,]
  ts = tsc(ip, theta)
	os = apply(resp, 1, sum, na.rm=TRUE)
	os = os[or]
	plot(theta[,1], ts[,1], type="l", lwd = 2, xlab="Estimated ability", ylab="Score",
		main="Observed and predicted scores")
	points(theta[,1], ts[,1]-ts[,2], type="l")
	points(theta[,1], ts[,1]+ts[,2], type="l")
	points(theta[,1], os, co="red")
}  
