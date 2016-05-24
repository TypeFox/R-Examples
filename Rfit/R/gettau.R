#' Estimate of the scale parameter tau
#' 
#' An estimate of the scale parameter tau is needed for the standard errors of
#' the coefficents in rank-based regression.
#' 
#' This is the confidence interval type estimate of the scale parameter tau
#' developed my Koul, Sievers, and McKean (1987). This estimate is also
#' discussed in Section 3.7.1 of Hettmansperger and McKean (1998). One of these
#' function is called in rfit.  The default is to use the faster FORTRAN
#' version. The R version can be more precise in small samples, but also can be
#' much slower especially when sample sizes are large.
#' 
#' @aliases gettau gettauF0
#' @param ehat full model residuals
#' @param p number of regression coefficents
#' @param scores object of class scores, defaults to Wilcoxon scores
#' @param delta confidence level
#' @param hparm Joe's hparm
#' @return Length one numeric object.
#' @author Joseph McKean, John Kloke
#' @seealso \code{\link{rfit}}
#' @references Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust
#' Nonparametric Statistical Methods, 2nd ed.}, New York: Chapman-Hall.
#' 
#' Koul, H.L., Sievers, G.L., and McKean, J.W. (1987) An esimator of the scale
#' parameter for the rank analysis of linear models under general score
#' functions, \emph{Scandinavian Journal of Statistics}, 14, 131-141.
#' @export gettau
gettau = function(ehat,p=0,scores=Rfit::wscores, delta = 0.8, hparm = 2, ...){

     n = length(ehat)
     asc = getScores(scores, 1:n/(n+1))
     ascpr = getScoresDeriv(scores, 1:n/(n+1) )
     temph = hstarreadyscr(ehat,asc,ascpr)
     abdord = temph$absdifford
     wtord = temph$wtsord
     const = temph$cn
     templ = looptau(delta,abdord,wtord,const,n)
     tn =  templ$quan/sqrt(n)
     pn = hstar(abdord,wtord,const,n,tn)
     tauscr = ((asc[n] - asc[1])*pn)/(2*tn)
     tauscr = sqrt(n/(n-p))*(1/tauscr)
     w = rep(0,n)
     stan = (ehat-median(ehat))/mad(ehat)
     w[abs(stan) < hparm] = 1
     hubcor = sum(w)/n
     if(hubcor < .000001){hubcor = .000001}
     fincor = 1 + (((p)/n)*((1-hubcor)/hubcor))
     tauscr = fincor*tauscr
     tauscr
#     list(tauscr=tauscr,asc=asc,ascpr=ascpr)
}
