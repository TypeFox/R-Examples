`print.saws` <-
function (x, digits = NULL, ...) 
{
    if (is.null(digits)){ 
        digits <- options()$digits
    } else options(digits = digits)
    cat("\n", "SAWS: Small sample Adjusted Wald-type test using Sandwich variance")
    cat("\n")
    cat("\nOriginal Call:\n")
    dput(x$originalCall)
     cat("\n")

    method.message<-switch(x$method,
        dm="dm (model based variance, chi squared test)",
        d1="d1 (usual unadjusted sandwich test)",
        d2="d2 (unadjusted sandwich variance, F test, df=dhat)",
        d3="d3 (unadjusted sandwich variance, F test, df=dtilde)",
        d4="d4 (adjusted sandwich variance, F test, df=dhat_H)",
        d5="d5 (adjusted sandwich variance, F test, df=dtilde_H)")
    cat("\n","SAWS method=",method.message)
  #   cat("\n","(see Fay and Graubard, 2001, Biometrics, 1198-1206)")
     cat("\n")
     cat("\n")


        
    r<-dim(x$test)[[1]]
    ## remember if x$test is r X p, then length(coef)=r
    ## and x$coefficients = t(test) %*% beta
    ## where beta=original coefficients
    ## so lenght(x$coefficients)= r
    ## previously I had, p<-length(x$coefficients)
    ## Now (April 19, 2012) it is corrected
    p<-dim(x$test)[[2]]
    output<-matrix(NA,r,4)
    output[,1]<-x$coefficients
    output[,2:3]<-x$conf.int
    output[,4]<-x$p.value 

    dimnames(output)<-list(dimnames(x$test)[[1]],c("Estimate",
        paste("Lower",100*attr(x$conf.int,"conf.level"),"% CL"),
        paste("Upper",100*attr(x$conf.int,"conf.level"),"% CL"),
        "2-sided p-value"))

     if (r==p && all(x$test==diag(p))){
        print(output)
        cat("--- p-values associated with tests that estimates are different from zero")
        cat("\n")
    } else {
        cat("\n","Test Matrix:")
        cat("\n")
        print(x$test)
        cat("\n","beta0:")
        cat("\n")
        cat(x$beta0)
        cat("\n")
        print(output)
        cat("--- p-values and confidence limits associated with: test %*% (Estimate-beta0)=0")    
    }
    if (x$method=="dm"){
        cat("--- scale=",x$scale, ", where Var(Y) is estimated by scale times variance function of mean")
    }
    invisible(x)
}

