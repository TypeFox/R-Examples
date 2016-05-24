FiellerC <-
function (Data) {
lmResults <- lm(y ~ x1 + x2 + x1*x2,data=Data) # fit lm() 
# extract coefficients
coefficient <- summary(lmResults)$coefficients 
# extract variance covariance matrix
covariance <- vcov(lmResults)
B2 <- coefficient[3,1] # estimation of B2
B3 <- coefficient[4,1] # estimation of B3
COV22 <- covariance[3,3] # variance of B2
COV33 <- covariance[4,4] # variance of B3
COV23 <- covariance[3,4] # covariance between B2 and B3

k <- 3.84 # 95% percentile for chi square with df=1
# Quadratic Equation: a*x^2 + b*x + c = 0
a <- B3^2 - k*COV33
b <- 2*(B2*B3-k*COV23)
c <- B2^2 - k*COV22
delta <- b^2-4*a*c # determinant of quadratic equation
# select only a>0 (convex) and delta>0 (real solution)
#,which give confidence interval of the form (LowCI,UpperCI). 
if (a>0 && delta>0) {
LowCI <- (-b-sqrt(b^2-4*a*c))/(2*a)
UpperCI <- (-b+sqrt(b^2-4*a*c))/(2*a)
}
else {
return(-1) # function returns -1 to indicate the case where
                       # the above condition is not satisfied.  
}
# collect value into a list
results <- list(LowCI = LowCI, UpperCI = UpperCI)
return(results)
}
