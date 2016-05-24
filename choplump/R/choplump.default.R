`choplump.default` <-
function(x, y, alternative = c("two.sided", "less", "greater"), 
    use.ranks=TRUE, exact = NULL, method=NULL, methodRule=methodRule1, 
    methodRuleParms=c(10^4), nMC=10^4-1,seed=1234321,
    printNumCalcs=TRUE,...){

    W<-c(x,y)
    if (any(W<0)) stop("choplump test not defined for negative responses")
    Z<-c(rep(0,length(x)),rep(1,length(y)))
    if (!any(W==0)) warning("No W values equal 0")

    M<-length(W[W!=0])
    N<-length(W)
    n1<-sum(Z)
    k<- N-M
    num.calcs<- sum( choose(M,max(0,n1-k):min(n1,M)) )
    if (is.null(method)) method<-methodRule1(W,Z,exact, methodRuleParms)
    method.OK<-(method=="approx" | method=="exact" | method=="exactMC")
    if (!method.OK) stop("method not one of: 'approx', 'exact', or 'exactMC'")
    if (method=="exact" | method=="exactMC") exact<-TRUE

    if (printNumCalcs){
        if (method=="exact"){
            cat(paste("calculating exact test...\nrequires ",num.calcs," evaluations of test statistic\n"))
            flush.console()
        } else if (method=="exactMC"){
            cat(paste("calculating exact test by Monte Carlo...\nrequires ",nMC," evaluations of test statistic\n"))
            flush.console()
        }
    }
    p.values<-switch(method,
        approx=choplumpApprox(W,Z,use.ranks),
        exact=choplumpExact(W,Z,use.ranks),
        exactMC=choplumpExactMC(W,Z,use.ranks,nMC,seed))
 
    alternative <- match.arg(alternative)
    PVAL <- switch(alternative,two.sided=p.values["p.2sided"],greater=p.values["p.upper"],
        less=p.values["p.lower"]) 
    if (use.ranks & method=="exact") METHOD<-"Exact Choplump Rank Test"
    else if (use.ranks & method=="exactMC") METHOD<-"Exact Choplump Rank Test by Monte Carlo"
    else if (use.ranks & method=="approx") METHOD<-"Asymptotic Choplump Rank Test"
    else if (!use.ranks & method=="exact") METHOD<-"Exact Choplump Difference in Means Test"
    else if (!use.ranks & method=="exactMC") METHOD<-"Exact Choplump Difference in Means Test by Monte Carlo"
    else if (!use.ranks & method=="approx") METHOD<-"Asymptotic Choplump Difference in Means Test"
    xname<-deparse(substitute(x))
    yname<-deparse(substitute(y))
    if (length(xname)>1) xname<-c("x")
    if (length(yname)>1) yname<-c("y")
   
    DNAME <- paste(xname, "and", yname)
    
    OUT <- list(statistic = NULL, parameter = NULL, p.value = as.numeric(PVAL), 
        null.value = NULL, alternative = alternative, method = METHOD, 
        data.name = DNAME, p.values=p.values)
    class(OUT) <- "htest"
    return(OUT)
}

