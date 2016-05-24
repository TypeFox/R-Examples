#kwayr = function(levs,data,ci=F,ciconf=.95,cimat=c(0)){

# jk 06-19-2011: removed confidence interval 


#' Internal Functions for K-Way analysis of variance
#' 
#' These are internal functions used to construct the robust anova table.
#' 
#' 
#' @aliases kwayr cellx khmat pasteColsRfit redmod subsets
#' @param levs vector of levels corresponding to each of the factors
#' @param data data matrix in the form y, factor 1,..., factor k
#' @param X n x k matrix where the columns represent the levels of the k
#' factors.
#' @param levsind Ask Joe
#' @param permh Ask Joe
#' @param x n x k matrix where the columns represent the levels of the k
#' factors.
#' @param xmat n x p full model design matrix
#' @param amat Ask Joe
#' @param k Ask Joe
#' @param sep Seperator used in pasteColsRfit
#' @note Renamed pasteCols of library plotrix written by Jim Lemon et. al. June
#' 2011 under GPL 2
#' @author Joseph McKean, John Kloke
#' @seealso \code{\link{raov}}
#' @references Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust
#' Nonparametric Statistical Methods, 2nd ed.}, New York: Chapman-Hall.
#' 
#' Hocking, R. R. (1985), \emph{The Analysis of Linear Models}, Monterey,
#' California: Brooks/Cole.
#' @export kwayr
kwayr = function(levs,data) {
#
#   Input:
#     levs = vector of levels corresponding to the factors A, B, C, etc.
#     data = matrix whose first column is the response variable then the succeeding columns
#            are the level of the first factor, second factor, etc.
#            For example, the Data 
#              for a three way is of the form:
#              Y_ijkl i j k 
#            NOTE: The data must be presorted so that k runs the fastest, than j, than i
#     CONTRASTS = T if confidence intervals are to be obtained.
#     CONTRASTSCONF is the confidence coefficient for the confidence intervals for
#                   the contrasts.
#     CONTRASTSMAT is the matrix of contrasts, one row for each contrast.
#

    nf = length(levs)
    y = data[,1]
    levsind = data[,2:(nf+1)]
    n = length(y)
    xfull = cellx(data[,2:(nf+1)])
    listtests = subsets(nf)
    p = length(xfull[1,])
     xuse = xfull
     fitF = rfit(y~xfull-1)
     iflagq = 0

    drf = disp(fitF$betahat, xuse, fitF$y, fitF$scores)
    ehatr = fitF$resid
    dfull = n - p
    nt = length(listtests[,1])
    tab2 = matrix(rep(0,nt*5),ncol=5)
    tab3 = matrix(rep(0,nt*5),ncol=5)

    for(i in 1:nt){
      permh = listtests[i,]
      hmat = khmat(levsind,permh)
      q = length(hmat[,1])
#      see = cbind(rep(i,q),hmat)
#      write(t(see),ncol=13,append=T,file="allh.dat")
      xred = redmod(xfull,hmat)
      fitr = rfit(y~xred-1)
      drr = disp(fitr$betahat, xred, fitr$y, fitr$scores)
      rd = drr - drf
      ft = (rd/q)/(fitF$tauhat/2)
      pv = 1 - pf(ft,q,dfull)
      tab2[i,] = c(q,rd,(rd/q),ft,pv)
#      tab3[i,] = c(dfull,0,(fitF$tauhat/2),drr,drf)
    }
	list(tab=tab2,fit=fitF)
#         list(listtests=listtests,tab2=tab2,tab3=tab3,ci=ci)
 
}
