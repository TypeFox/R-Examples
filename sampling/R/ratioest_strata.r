ratioest_strata<-function(y,x,TX_strata,pik,strata,description=FALSE)
{
if (missing(x) | missing(y) | missing(pik) | missing(strata)) stop("incomplete input")
if(!is.vector(x)) x=as.vector(x)
if(length(y)!=length(x)) stop("x and y have different sizes")
str <- function(st, h, n) .C("str", as.double(st), as.integer(h), as.integer(n), s = double(n), PACKAGE = "sampling")$s
sample.size = length(y)
h = unique(strata)
s1 = 0
for (i in 1:length(h)) {
s=str(strata, h[i], sample.size)
ys=y[s==1]
xs=x[s==1]
r=ratioest(ys,xs,TX_strata[h[i]],pik[s==1])
s1 = s1 + r
if(description)
 {cat("Stratum ",h[i],", the ratio estimator is:",r,"\n")
  cat("Number of units:",sum(s),"\n")
    }  
  }
if(description) 
  cat("The ratio estimator is:\n")
s1
}


