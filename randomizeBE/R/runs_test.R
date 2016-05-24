###############################################################################
# p.value of the 2-sided runs test 
# 
# Author: dlabes
###############################################################################
# adapted from runs.test package lawstat
# Author(s): Wallace Hui, Yulia R. Gel, Joseph L. Gastwirth, Weiwen Miao
# and runs.test package tseries
# Author(s): A. Trapletti
# only the 2-sided p-value is of concern for our purposes
# in package randomizeBE

runs.pvalue <- function(y, pmethod=c("exact", "normal", "cc"))
{
  pmethod <- match.arg(pmethod)
  y <- na.omit(y)
  # do not dichotomize if already +1, -1 or some othe binary coding
  yuniq <- unique(y)
  if (length(yuniq)==2){
    s <- as.numeric(y)
    s[yuniq[1]==y] <- -1
    s[yuniq[2]==y] <- +1
  } else {
    med <- median(y, na.rm = TRUE)
# lawstat handling of values == med
# gives zero values if the first value(s) == med !!!
# annoying p-values for 3 sequences  
#    for (k in 2:length(y)) {
#      if ((y[k] == med) & (y[k - 1] < med)) {
#        y[k] = y[k - 1]
#      }
#      else if ((y[k] == med) & (y[k - 1] > med)) {
#        y[k] = y[k - 1]
#      }
#    }
    # -1 and +1
    s  <- sign(y - med)
    s[s==0] <- +1 # values==med were classified as +1
  }
  #print(s)
  n  <- length(s)
  # runs: see runs.test from package tseries
  R  <- 1 + sum(as.numeric(s[-1] != s[-n]))
  # number of above/below
  n1 <- sum(s == +1)
  n2 <- sum(s == -1)
  n  <- n1 + n2
  #cat("runs:",R," ns:",n1,n2,"\n") 
  E  <- 1 + 2*n1*n2/n
  s2 <- (2*n1*n2*(2*n1*n2 - n))/(n^2 * (n - 1))
  # with continuity corr. ? see SPSS npar tests
  # or SAS support: http://support.sas.com/kb/33/092.html
  ccf <- 0
  if (pmethod=="cc" | pmethod=="exact"){
    ccf <- ifelse((R - E) < 0, +0.5, -0.5)
  }
  statistic <- (R - E + ccf)/sqrt(s2)

  if((n>30 & n1>12 & n2>12) | pmethod != "exact"){
    # asymptotic normal
    pvalue <- 2 * pnorm(-abs(statistic))
  } else {
    pvalue <- pruns.exact(R, n1, n2, tail="2-sided")
  }
  pvalue
}

# exact conditional cumulative distribution
pruns.exact <- function(r, n1, n2, tail=c("2-sided", "lower", "upper"))
{
  tail <- match.arg(tail)
  
  n <- n1+n2
  
  if(r<=1) stop("Number of runs must be >1")
  if(r>n) stop("Number of runs must be <(n1+n2")
  if(n1<1 | n2<1) return(0) #??? is not random!

  E <- 1 + 2*n1*n2/n
  #cat("E:",E,"\n")
  
  denom <- choose(n,n1)
  # how long should we make the r vector?
  # in very unsymmetric cases only a few elements of
  # pp = density have values > 0 if rmax=n1+n2
  # number of runs possible: 2*m if n=m, 2*m+1 if m<n
  rmax <- ifelse(n1==n2, 2*n1, 2*min(n1,n2)+1)
  rv <- 2:rmax
  pp <- druns_nom(rv, n1, n2)
  #print(pp)
  # pL is p(R<=r) -> left/lower tail
  pL <- sum(pp[rv<=r])/denom 
  #pU is p(R>=r) -> right/upper tail
  pU <- 1 - sum(pp[rv<=(r-1)])/denom
  #print(abs(rv-E)>=abs(r-E))
  # Equn. 4.7 of the SPSS documentation
  p2 <- sum(pp[abs(rv-E)>=abs(r-E)])/denom
  #cat("pL:",pL,"\n")
  #cat("pU:",pU,"\n")
  #cat("p2:",p2,"\n")
  # Next is the rule from:
  # Gibbons "Nonparametric Methods for Quantitative Analysis"
  # 0.5 is to avoid p>1 if both pL and pU are >0.5
  p2min <- 2*min(c(pL, pU, 0.5))
  # we are using the SPSS approach wich takes into account the 
  # unsymmetric form of the distribution if n1 << n2
  if(tail=="2-sided") return(p2)
  if(tail=="lower")  return(pL)
  if(tail=="upper") return(pU)
}

# function for calculating the denominator of the runs distribution
druns_nom <- function(r, n1, n2){
# vector form, i.e. if r is a vector
#  k    <- r
#  pp   <- r
#  even <- 2*r%/%2==r
#  k[even]   <- r[even]/2
#  pp[even]  <- 2*choose(n1-1,k[even]-1)*choose(n2-1,k[even]-1)
#  k[!even]  <- (r[!even]-1)/2
# the odd case seems not vectorize propperly!  
#  pp[!even] <-   choose(n1-1,k[!even]-1)*choose(n2-1,k[!even]) 
#               + choose(n1-1,k[!even])*choose(n2-1,k[!even]-1)
  #print(r);print(k);print(pp)               
# original code
  pp <- vector(mode="numeric",length=length(r))
  for (i in seq_along(r)){
    if (2*r[i]%/%2==r[i]){
      # even 2*k
      k <- r[i]/2
      pp[i] <- 2*choose(n1-1, k-1)*choose(n2-1, k-1)
    } else {
      # odd 2*k+1
      k <- (r[i]-1)/2
      pp[i] <- choose(n1-1,k-1) * choose(n2-1,k) +
               choose(n1-1,k)   * choose(n2-1,k-1)
    }
  }  
  pp
}

# SPSS example
# x <- c(31, 23, 36, 43, 51, 44, 12, 26, 43, 75,  2,  3, 15, 18, 78, 24, 13, 27, 86, 61, 13,  7,  6,  8)
# should give exact    0.3008894
# asymptotic (with cc) 0.2966895
# SPSS example small dataset
# x <- c(1, 1, 1, 1, 0, 0, 0, 0, 1, 1)
# exact p:  0.071 (3 runs, n1=4, n2=6)
# asymp. p: 0.106
