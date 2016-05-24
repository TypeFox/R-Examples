VuongClarke <- function(obj1, obj2, sig.lev = 0.05){

l1 <- obj1$logLik
l2 <- obj2$logLik

if(l1==l2) stop("The two competing models have identical log-likelihoods!")

p1 <- obj1$t.edf
p2 <- obj2$t.edf
n  <- obj1$n
n1 <- obj2$n

if(n != n1) stop("The two competing models have different sample sizes.")

li12 <- obj1$fit$l.par - obj2$fit$l.par
w  <- sqrt(var(li12)*(n-1))
vt <- (l1-l2-(p1-p2)/2*log(n))/w  
li12b <- li12-(p1-p2)/(2*n)*log(n) 
b     <- sum(li12b>0) 

if(abs(vt) <= qnorm(1 - sig.lev/2)) dvt <- " Vuong's test: it is not possible to discriminate between the two models."
if(vt >  qnorm(1 - sig.lev/2))      dvt <- " Vuong's test: Model 1 is preferred over Model 2."
if(vt < -qnorm(1 - sig.lev/2))      dvt <- " Vuong's test: Model 2 is preferred over Model 1."

 
db <- "Clarke's test: it is not possible to discriminate between the two models."   

if(b >= n/2){ pvalue <- 2 * (1 - pbinom(b - 1, n, 0.5))
              if(pvalue <= sig.lev) db <- "Clarke's test: Model 1 is preferred over Model 2."
            }
if(b < n/2) { pvalue <- 2 * (pbinom(b, n, 0.5))
              if(pvalue <= sig.lev) db <- "Clarke's test: Model 2 is preferred over Model 1."
            }

cat("\n",dvt,"\n",sep="")
cat(db,"\n\n",sep="")

}






