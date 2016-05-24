test.theta <-
function(x, y, alternative = c("two.sided", "less", "greater"),
       theta = 1, B = 1000, conf.level = 0.95){
    alternative <- match.arg(alternative)
    DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(y))) 
    x <- x[!is.na(x)]
    m <- length(x)
    if (m < 1) 
        stop("not enough 'x' data")
    PVAL <- NULL
    y <- y[!is.na(y)]
    n <- length(y)
    if (n < 1) 
        stop("not enough 'y' data")
    METHOD <- "Test for Proportional Odds Rate theta"
    N<-m+n
    lambda<-m/N
    Fx<-ecdf(x)
    theta0<-dd.est(x,y)#newton.theta(Fx(y))
    cat(theta0,"\n")
    p0<-phi(N, theta0, lambda)/N
    res<-try(mrle.sporm(x, y, theta0, p0), TRUE)
    theta.hat<-res$theta
    cat(theta.hat,"\n")
    STATISTIC<-theta.hat
    names(STATISTIC)<-"Estimated theta"
    Theta.b<-NULL
    sim<-function(){
        u<-runif(m)
        v<-runif(n)
        v<-v/(theta+(1-theta)*v)
        theta0<-dd.est(u,v)#newton.theta(ecdf(u)(v))
        p0<-phi(N, theta0, lambda)/N
        res<-try(mrle.sporm(u, v, theta0, p0), TRUE)
        theta.b<-res$theta
        theta.b 
    }
    Theta.b <- lapply(1:B, function(i) try(sim(), TRUE))
    Theta.b <-unlist(Theta.b[sapply(Theta.b, function(x) !inherits(x, "try-error"))])
    PVAL <- switch(alternative, two.sided = mean(abs(Theta.b-theta)>abs(theta.hat-theta)), 
                greater = mean(Theta.b>theta.hat), less = mean(Theta.b<theta.hat))
    nm_alternative <- switch(alternative, two.sided = "true difference in theta values is not equal to 0", 
        less = "true difference in theta values is less than 0", 
        greater = "true difference in theta values is greater than 0")
 
    RVAL <- list(statistic = STATISTIC, p.value = PVAL, alternative = nm_alternative, 
        method = METHOD, data.name = DNAME)
    class(RVAL) <- "htest"
    return(RVAL)
}
