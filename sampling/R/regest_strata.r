regest_strata<-function(formula,weights,Tx_strata,strata,pikl,sigma=rep(1,length(weights)),description=FALSE)
{
cl <- match.call() 
mf <- match.call(expand.dots = FALSE)
m <- match(c("formula", "weights"), names(mf), 0)
mf <- mf[c(1, m)]
mf$drop.unused.levels <- TRUE
mf[[1]] <- as.name("model.frame")
mf <- eval(mf, parent.frame())
mt <- attr(mf, "terms")
y <- model.response(mf, "numeric")
w <- as.vector(model.weights(mf))
x <- model.matrix(mt, mf, contrasts)
str <- function(st, h, n) .C("str", as.double(st), as.integer(h), as.integer(n), s = double(n), PACKAGE = "sampling")$s
sample.size = length(y)
h = unique(strata)
s1 = 0
for (i in 1:length(h)) {
s=str(strata, h[i], sample.size)
ys=y[s==1]
xs=x[s==1,]
r=regest(ys~xs-1,Tx=Tx_strata[h[i]],weights=weights[s==1],pikl=pikl[s==1,s==1],n=length(s[s==1]),sigma[s==1])
est=r$regest
s1 = s1 + est
if(description)
 {cat("Stratum ",h[i],", the regression estimator is:",est,"\n")
  cat("Number of units:",sum(s),"\n")
  cat("Beta coefficient(s):", r$coefficients,"\n")
  cat("Std. error:", r$std_error,"\n")
  cat("t-value:", r$t_value, "\n")
  cat("p_value:",r$p_value,"\n")  
  cat("cov_matrix:\n")
  print(r$cov_matrix)
    }  
  }
if(description) 
  cat("The regression estimator is:\n")
s1
}

