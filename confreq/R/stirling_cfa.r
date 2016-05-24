#' @title  Approximation to the binomial using Stirling's Formula
#' @export stirling_cfa
#' @description Calculates the binomial aproximation using stirling's formula (Version of function: V 1.0 - November 2013)
#' @author R.W. Alexandrowicz
#' @param observed a integer vector with observed freqencies
#' @param expected a vector giving the expected frequencies. expected can be set to \code{expected=NULL} if an vector of cell probabilities is given in argument \code{p}. 
#' @param n number of trials (scalar) default is \code{n = sum(observed) }.
#' @param p a vector of cell probabilities. If p is not NULL the argument \code{expected} is ignored and this vector p of cell probabilities is used for calculatio instead of expected counts
#' @param cum a logical - computation of cumulative density. If \code{cum=TRUE} (default) computes tail probability. If \code{cum=FALSE} computes prob. only for one cell (i.e. execute stircore only).
#' @param verb logical - verbose results: If \code{verb=TRUE} (default) builds a results table. If \code{verb=FALSE} returns vector of cell p-values only.
#' @details \itemize{ \item{Vector p must be of same length as observed _or_ p may be a scalar (e.g. in case of the zero-order CFA).}
#' \item {The routine autoselects the upper or lower tail: 
#' \itemize{ \item{if obs > exp then sum obs:n} \item{ else              sum 0:obs}}}
#' \item {The stirling approximation cannot be evaluated if the observed frequency is 0 or n. Therefore, the proposal of A. von Eye (20xx) is adopted, taking the sum up to 1 or n-1, respectively.}}
#' @references 
#' von Eye, A. (2002). \emph{Configural Frequency Analysis. Methods, Models, and Applications.} Mahwah, NJ, LEA.

  stirling_cfa = function(observed,expected=NULL,n=sum(observed),p=NULL,cum=T,verb=T) {
    ## ergänzungen JHH
    if( (length(p)==0) & (length(expected)==0) ) stop("either expected counts or cell probabilities must be given")
    if( (length(p)!=0) ){p<-p; cat("vector of cell probabilities is used instead of expected counts", "\n")}
    if(length(p)==0){p<-expected/n}
    
    ## ENDE ergänzungen JHH

     stircore_x = function(n,observed,p) {   # core function (working horse)
        sqrt(n/(2*pi*observed*(n-observed))) * (n*p/observed)^observed * (n*(1-p)/(n-observed))^(n-observed)  
     }                                # end of core function

     stircore = function(n,observed,e) {     # core function (working horse)
        sqrt(n/(2*pi*observed*(n-observed))) * (e/observed)^observed * ((n-e)/(n-observed))^(n-observed)
     }                                # end of core function

     stopifnot(is.vector(n) & length(n)==1)  # n must be scalar
  
     nc = length(observed)                   # number of cells (= different pattern)
     if (length(p)==1) p = rep(p,nc)  # expand p if scalar (for loop)

     if(length(expected)==0){e <- n*p} else {e <- expected} # exp freq. JHH: only when no exp. counts are given 
     
     up = observed >= e               # upper or lower tail prob
     lim = rep(1,nc)                  # set lower limits
     lim[up] = n-1                    # set upper limits
     
     pval = NULL                      # init result
     if (cum) {                       # tail area (like pbinom)
        # for (i in 1:nc) pval[i] = sum(stircore(n,o[i]:lim[i],p[i]))
          for (i in 1:nc) pval[i] = sum(stircore(n,observed[i]:lim[i],e[i]))
          if (verb) {                 # verbose output: build result table
             res= cbind(observed,e,up,lim,pval)   # build table
             colnames(res) = c("obs.freq","exp.freq","upper tail","limit","p-value")
             }                        # end verbose
          else res = pval             # vector of p-values only
         }                            # end if (cum density)
     else {                           # only density (like dbinom)
          dens = stircore(n,observed,e)      # p(K=k|n,p)
          if (verb) {                 # verbose output
             res= cbind(observed,e,dens)     # build table
             colnames(res) = c("obs.freq","exp.freq","density")
             }                        # end verbose
          else res = dens             # vector of densities only
          }                           # end else (density only)
     return(res)                      # return result table
     
  } # end of function
  
