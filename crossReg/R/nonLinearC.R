nonLinearC <-
function (Data,startingValue) {
catchError <- tryCatch (
nlsResults <- nls(y ~ A0 + B1*(x1-C)+B2*((x1-C)*x2)
,start=startingValue,data=Data),
error = function(e) e
)
# return -1 if nls() return error
if(inherits(catchError, "error")) return(-1) 
coefficient <- summary(nlsResults)$coefficients # extract coefficients
C_Hat <- coefficient[4,1] 
SE <- coefficient[4,2]
# Use the standard error (SE) from nonlinear regression 
# to construct a confidence interval under the assumption that
      # the sampling distribution of C_Hat is N(C,SE^2)
LowCI <- C_Hat-1.96*SE # lower bound of CI
UpperCI <- C_Hat+1.96*SE # upper bound of CI
# collect values into a list
results <- list(C_Hat = C_Hat, SE = SE, LowCI = LowCI, UpperCI =     
      UpperCI) 
print(coefficient)
return(results) # return list
}
