power.mcnemar.test <- function(n = NULL, paid = NULL, psi = NULL, sig.level = 0.05, power = NULL,
                               alternative = c("two.sided", "one.sided"),
                               method = c("normal", "exact")) {

    if (sum(sapply(list(n, paid, psi, power, sig.level), is.null)) != 1) 
        stop("exactly one of 'n', 'paid', 'psi', 'power', and 'sig.level' must be NULL")
    if (!is.null(sig.level) && !is.numeric(sig.level) || any(0 > sig.level | sig.level > 1)) 
        stop("'sig.level' must be numeric in [0, 1]")
    alternative <- match.arg(alternative)
    method <- match.arg(method)
    tside <- switch(alternative, one.sided = 1, two.sided = 2)

        
    f <- function(n, paid, psi, sig.level, power) {
        pd <- (1 + psi)*paid
        d <- (psi-1)*paid
        r <- ceiling(log(sig.level)/log(.5))
        n <- ceiling(n)

        power <- sum(
            
         sapply(r:n, function(x) { 
            hhh <- sapply(0:(which.min(pbinom(0:x, size=x, prob=.5) <= sig.level)-2),
                   function(y) {
                lgamma(n+1) - lgamma(n-x+1) - lgamma(y+1) - lgamma(x-y+1) + (n-x)*log(1-pd) + y*log((d+pd)/2) + (x-y)*log((pd-d)/2)
            } )
#            cat("====\n")
#            print(res)
#            print(hhh)
#            print(exp(hhh))
            sum(exp(hhh))
        })
        )
        power
    }

    if (method=="normal") {
        p.body <- quote( pnorm (
            (sqrt(n * paid) * (psi-1) - qnorm(sig.level/tside, lower.tail=FALSE)*sqrt(psi+1)) / sqrt((psi+1) - paid*(psi-1)^2)))
    } else {
        p.body <- quote( f(n, paid, psi, sig.level, power) )
    }

    
    



    

#    print(f(n))


#    if (is.null(power)) {
#        power <- uniroot(function(power) eval(n.body) - n, c(0.001, 1-1e-10))$root        
#    } else if (is.null(n)) {
#        n <- eval(n.body)
#    } else if (is.null(paid)) 
#        paid <- uniroot(function(paid) eval(n.body) - n, c(0, 1-p10-1e-10))$root
#    else if (is.null(psi)) 
#        psi <- uniroot(function(psi) eval(n.body) - n, c(1e-10, 10))$root
#    else if (is.null(sig.level)) 
#        sig.level <- uniroot(function(sig.level) eval(n.body) - 
#            n, c(1e-10, 1 - 1e-10))$root

    if (is.null(power)) {
        power <- eval(p.body)
    } else if (is.null(n)) {
        n <- uniroot(function(n) eval(p.body) - power, c(ceiling(log(sig.level)/log(.5)), 1e+07))$root        
    } else if (is.null(paid)) 
        paid <- uniroot(function(paid) eval(p.body) - power, c(0, 1-psi-1e-10))$root
    else if (is.null(psi)) 
        psi <- uniroot(function(psi) eval(p.body) - power, c(1e-10, 1/paid-1-1e-10))$root
    else if (is.null(sig.level)) 
        sig.level <- uniroot(function(sig.level) eval(p.body) - 
            power, c(1e-10, 1 - 1e-10))$root
    else stop("internal error", domain = NA)
    NOTE <- "n is number of pairs"
    METHOD <- paste("McNemar paired comparison of proportions", ifelse(method=="normal", "approximate", "exact") ,"power calculation")
    structure(list(n = n, paid = paid, psi = psi, sig.level = sig.level, 
        power = power, alternative = alternative, note = NOTE, 
        method = METHOD), class = "power.htest")

}
