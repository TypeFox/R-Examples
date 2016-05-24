"neill.test" <- function(
object, grouping, method = c("c-finest", "finest", "percentiles"), breakp = NULL, display = TRUE)
{
    method <- match.arg(method)
    
    x <- object$"dataList"$"dose"
    noCluster <- floor(length(x)/2) 
    
    if (missing(grouping))
    {
        if (method == "finest")
        {      
            lenx <- length(x)
            grouping <- floor((1+1:lenx)/2)
            grouping[lenx] <- grouping[lenx-1]
        } 
        if (method == "c-finest")
        { 
            for (i in noCluster:(length(coef(object)) + 1))
            {
                grouping <- cutree(hclust(dist(x)), k = i)
                if (all(tapply(x, grouping, length) > 1)) {break}
            }
        }
        if (method == "percentiles")
        {
            cutVar <- c(-Inf, quantile(x, c(0.2, 0.4, 0.6, 0.8)), Inf)
            grouping <- cut(x, cutVar)
        }
   
        if (!is.null(breakp))
        {
            cutVar <- c(-Inf, breakp, Inf)
            grouping <- cut(x, cutVar)   
        }
    }
    if (display)
    {
        neill.default(object, grouping, anovaDisplay = TRUE)
    } else {
        neill.default(object, grouping, anovaDisplay = FALSE)[1, 2]
    }
}

"neill.default" <- function(object, grouping, anovaDisplay)
{
    M <- length(unique(grouping))
    N <- nrow(eval(object$data))
    denDF <- N - M
    
    ## Checking the number of groups
    if (denDF <= 0)  # (N <= M) 
    {
        stop("Too many groups in 'grouping'")
    }
    p <- N - df.residual(object)
    numDF <- M - p
    
    if (numDF <= 0)  # (M <= p)
    {
        stop("Too few groups in 'grouping'")
    }

    if (anovaDisplay)
    {
        ## Print information on grouping:
#    tapply(polcurve$dist, grouping, function(x){paste(as.character(range(x)), collapse = "-")})
# The above line should in the clustering function
        cat("Grouping used\n\n")
        grTable <- table(grouping)
        dimnames(grTable) <- list(dimnames(grTable)$grouping)  
        # to avoid 'grouping' is being displayed
        print(grTable)
        cat("\n")    
    }
    
    ## Calculating the test statistic
    resVec <- residuals(object)
    resAver0 <- tapply(resVec, grouping, mean)
    resAver <- rep(resAver0, tapply(grouping, grouping, length))

    resDiff <- resVec - resAver
    F <- (denDF/numDF)*(sum(resAver*resAver)/(sum(resDiff*resDiff)))
    p <- pf(F ,numDF, denDF, lower.tail = FALSE)

    if (anovaDisplay)
    {
        ## Formatting for 'print.anova'
        dataFra <- data.frame(F, p)
        dimnames(dataFra) <- list("", c("F value", "p value"))
        structure(dataFra, heading = "Neill's lack-of-fit test\n", class = c("anova", "data.frame"))    
    } else {
        return(matrix(c(F, p), 1, 2))
    }
}

if (FALSE)
{

"sim.neill.test" <- function(noSim, noObs = 20, noRep = NULL, seed = 20070327, 
true = c("llogistic", "weibull", "expdecay"), ...)
{
    true <- match.arg(true)
    
    set.seed(seed)

    expdecay <-
    function(x, b = 4, d = 2.15)
    {
        d * exp(-b * x)
    }    
    llogist <- 
    function(x, c = 0.08, d = 2.15, b = 1.84, e = 0.2)
    {
        c+(d-c)/(1+(x/e)^b)
    }
    weibull <-
    function(x, c = 0.08, d = 2.15, b = 1, e = 0.2)
    {
        c+(d-c)*exp(-exp(b*(log(x)-log(e))))
    }
    genFct <- switch(true, "expdecay" = expdecay, "llogistic" = llogist, "weibull" = weibull) 

    if (!is.null(noRep))
    {
        xFct <- function(){rep( (1:noObs)/(noObs+1), rep(noRep, noObs))}
        errFct <- function(){rnorm(noObs * noRep, 0, sqrt(0.0054))}
    } else {
        xFct <- function(){sort(runif(noObs, 0, 1))}
        errFct <- function(){rnorm(noObs, 0, sqrt(0.0054))}
    }
        
    pVec <- rep(NA, noSim)
    for (i in 1:noSim)
    {
#        x <- sort(runif(noObs, 0, 1))
#        y <- genFct(x) + rnorm(noObs, 0, sqrt(0.0054))   
        x <- xFct()
        y <- genFct(x) + errFct() 
        fit <- multdrc(y ~ x)  # fitting log-logistic model
    
        pVec[i] <- neill.wrap(fit, x, ...)
    }
    pVec
}

    ### Design without replicates
    ## Simulation using finest clustering
    
    ## Examining the power of the test (two alternatives)
    pmNeill <- matrix(NA, 1000, 7)
    colnames(pmNeill) <- as.character(c(10, 25, 50, 75, 100, 250, 500))

    pmNeill[, 1] <- sim.neill.test(1000, 10, seed=21, true = "weibull")
    pmNeill[, 2] <- sim.neill.test(1000, 25, seed=22, true = "weibull")
    pmNeill[, 3] <- sim.neill.test(1000, 50, seed=23, true = "weibull")
    pmNeill[, 4] <- sim.neill.test(1000, 75, seed=24, true = "weibull")
    pmNeill[, 5] <- sim.neill.test(1000, 100, seed=25, true = "weibull")
    pmNeill[, 6] <- sim.neill.test(1000, 250, seed=26, true = "weibull")
    pmNeill[, 7] <- sim.neill.test(1000, 500, seed=27, true = "weibull")
    pFct <- function(x, alpha = 0.05) {sum(x < alpha)/1000}    
    apply(pmNeill, 2, pFct)

    pmNeill.2 <- matrix(NA, 1000, 7)
    colnames(pmNeill.2) <- as.character(c(10, 25, 50, 75, 100, 250, 500))

#    pmNeill.2[, 1] <- sim.neill.test(1000, 10, seed=11121, true = "expdecay")
    pmNeill.2[, 2] <- sim.neill.test(1000, 25, seed=11122, true = "expdecay")
    pmNeill.2[, 3] <- sim.neill.test(1000, 50, seed=11123, true = "expdecay")
    pmNeill.2[, 4] <- sim.neill.test(1000, 75, seed=11124, true = "expdecay")
    pmNeill.2[, 5] <- sim.neill.test(1000, 100, seed=11125, true = "expdecay")
    pmNeill.2[, 6] <- sim.neill.test(1000, 250, seed=11126, true = "expdecay")
#    pmNeill[, 7] <- sim.neill.test(1000, 500, seed=11127, true = "expdecay")
    pFct <- function(x, alpha = 0.05) {sum(x < alpha)/1000}    
    apply(pmNeill.2, 2, pFct)
    
    
    ## Examining the level of the test
    pmNeillt <- matrix(NA, 1000, 7)
    colnames(pmNeillt) <- as.character(c(10, 25, 50, 75, 100, 250, 500))
    
    pmNeillt[, 1] <- sim.neill.test(1000, 10, seed=31)
    pmNeillt[, 2] <- sim.neill.test(1000, 25, seed=32)
    pmNeillt[, 3] <- sim.neill.test(1000, 50, seed=33)
    pmNeillt[, 4] <- sim.neill.test(1000, 75, seed=34)
    pmNeillt[, 5] <- sim.neill.test(1000, 100, seed=35)
    pmNeillt[, 6] <- sim.neill.test(1000, 250, seed=36)
    pmNeillt[, 7] <- sim.neill.test(1000, 500, seed=37)
    apply(pmNeillt, 2, pFct)
    

    ## Simulation using fivenum percentiles
    
    ## Examining the power of the test
    pmNeill2 <- matrix(NA, 1000, 7)
    colnames(pmNeill2) <- as.character(c(10, 25, 50, 75, 100, 250, 500))

    pmNeill2[, 1] <- sim.neill.test(1000, 10, seed=421, true = "weibull", method = "percentiles")
    pmNeill2[, 2] <- sim.neill.test(1000, 25, seed=422, true = "weibull", method = "percentiles")
    pmNeill2[, 3] <- sim.neill.test(1000, 50, seed=423, true = "weibull", method = "percentiles")
    pmNeill2[, 4] <- sim.neill.test(1000, 75, seed=424, true = "weibull", method = "percentiles")
    pmNeill2[, 5] <- sim.neill.test(1000, 100, seed=425, true = "weibull", method = "percentiles")
    pmNeill2[, 6] <- sim.neill.test(1000, 250, seed=426, true = "weibull", method = "percentiles")
    pmNeill2[, 7] <- sim.neill.test(1000, 500, seed=427, true = "weibull", method = "percentiles")
    pFct <- function(x, alpha = 0.05) {sum(x < alpha)/1000}    
    apply(pmNeill2, 2, pFct)

    ## Examining the level of the test
    pmNeill2t <- matrix(NA, 1000, 7)
    colnames(pmNeill2t) <- as.character(c(10, 25, 50, 75, 100, 250, 500))
    
    pmNeill2t[, 1] <- sim.neill.test(1000, 10, seed=431, method = "percentiles")
    pmNeill2t[, 2] <- sim.neill.test(1000, 25, seed=432, method = "percentiles")
    pmNeill2t[, 3] <- sim.neill.test(1000, 50, seed=433, method = "percentiles")
    pmNeill2t[, 4] <- sim.neill.test(1000, 75, seed=434, method = "percentiles")
    pmNeill2t[, 5] <- sim.neill.test(1000, 100, seed=435, method = "percentiles")
    pmNeill2t[, 6] <- sim.neill.test(1000, 250, seed=436, method = "percentiles")
    pmNeill2t[, 7] <- sim.neill.test(1000, 500, seed=437, method = "percentiles")
    apply(pmNeill2t, 2, pFct)    
    
    ## Simulation using pairs grouping
    
    ## Examining the power of the test (two alternatives)
    pmNeill3 <- matrix(NA, 1000, 7)
    colnames(pmNeill3) <- as.character(c(10, 25, 50, 75, 100, 250, 500))

    pmNeill3[, 1] <- sim.neill.test(1000, 10, seed=521, true = "weibull", method = "finest")
    pmNeill3[, 2] <- sim.neill.test(1000, 25, seed=522, true = "weibull", method = "finest")
    pmNeill3[, 3] <- sim.neill.test(1000, 50, seed=523, true = "weibull", method = "finest")
    pmNeill3[, 4] <- sim.neill.test(1000, 75, seed=524, true = "weibull", method = "finest")
    pmNeill3[, 5] <- sim.neill.test(1000, 100, seed=525, true = "weibull", method = "finest")
    pmNeill3[, 6] <- sim.neill.test(1000, 250, seed=526, true = "weibull", method = "finest")
    pmNeill3[, 7] <- sim.neill.test(1000, 500, seed=527, true = "weibull", method = "finest")
    pFct <- function(x, alpha = 0.05) {sum(x < alpha)/1000}    
    apply(pmNeill3, 2, pFct)

    pmNeill3.2 <- matrix(NA, 1000, 7)
    colnames(pmNeill3.2) <- as.character(c(10, 25, 50, 75, 100, 250, 500))

#    pmNeill3.2[, 1] <- sim.neill.test(1000, 10, seed=521, true = "expdecay", method = "finest")
    pmNeill3.2[, 2] <- sim.neill.test(1000, 25, seed=522, true = "expdecay", method = "finest")
    pmNeill3.2[, 3] <- sim.neill.test(1000, 50, seed=523, true = "expdecay", method = "finest")
    pmNeill3.2[, 4] <- sim.neill.test(1000, 75, seed=524, true = "expdecay", method = "finest")
    pmNeill3.2[, 5] <- sim.neill.test(1000, 100, seed=525, true = "expdecay", method = "finest")
    pmNeill3.2[, 6] <- sim.neill.test(1000, 250, seed=526, true = "expdecay", method = "finest")
#    pmNeill3[, 7] <- sim.neill.test(1000, 500, seed=527, true = "expdecay", method = "finest")
    pFct <- function(x, alpha = 0.05) {sum(x < alpha)/1000}    
    apply(pmNeill3.2, 2, pFct)

    ## Examining the level of the test
    pmNeill3t <- matrix(NA, 1000, 7)
    colnames(pmNeill3t) <- as.character(c(10, 25, 50, 75, 100, 250, 500))
    
    pmNeill3t[, 1] <- sim.neill.test(1000, 10, seed=531, method = "finest")
    pmNeill3t[, 2] <- sim.neill.test(1000, 25, seed=532, method = "finest")
    pmNeill3t[, 3] <- sim.neill.test(1000, 50, seed=533, method = "finest")
    pmNeill3t[, 4] <- sim.neill.test(1000, 75, seed=534, method = "finest")
    pmNeill3t[, 5] <- sim.neill.test(1000, 100, seed=535, method = "finest")
    pmNeill3t[, 6] <- sim.neill.test(1000, 250, seed=536, method = "finest")
    pmNeill3t[, 7] <- sim.neill.test(1000, 500, seed=537, method = "finest")
    apply(pmNeill3t, 2, pFct)    
        

    ### Design with replicates
    ## Simulation using finest clustering
    
    ## Examining the power of the test
    pmNeill.rep <- matrix(NA, 1000, 6)
    colnames(pmNeill.rep) <- as.character(c(1, 2, 3, 5, 10, 20))

    pmNeill.rep[, 1] <- sim.neill.test(1000, 10, 1, seed=200821, true = "weibull")
    pmNeill.rep[, 2] <- sim.neill.test(1000, 10, 2, seed=200822, true = "weibull")
    pmNeill.rep[, 3] <- sim.neill.test(1000, 10, 3, seed=200823, true = "weibull")
    pmNeill.rep[, 4] <- sim.neill.test(1000, 10, 5, seed=200824, true = "weibull")
    pmNeill.rep[, 5] <- sim.neill.test(1000, 10, 10, seed=200825, true = "weibull")
    pmNeill.rep[, 6] <- sim.neill.test(1000, 10, 20, seed=200826, true = "weibull")
    pFct <- function(x, alpha = 0.05) {sum(x < alpha)/1000}    
    apply(pmNeill.rep, 2, pFct)

    ## Examining the level of the test
    pmNeillt.rep <- matrix(NA, 1000, 6)
    colnames(pmNeillt.rep) <- as.character(c(1, 2, 3, 5, 10, 20))
    
    pmNeillt.rep[, 1] <- sim.neill.test(1000, 10, 1, seed=200831)
    pmNeillt.rep[, 2] <- sim.neill.test(1000, 10, 2, seed=200832)
    pmNeillt.rep[, 3] <- sim.neill.test(1000, 10, 3, seed=200833)
    pmNeillt.rep[, 4] <- sim.neill.test(1000, 10, 5, seed=200834)
    pmNeillt.rep[, 5] <- sim.neill.test(1000, 10, 10, seed=200835)
    pmNeillt.rep[, 6] <- sim.neill.test(1000, 10, 20, seed=200836)
    apply(pmNeillt.rep, 2, pFct)
    

    ## Simulation using fivenum percentiles
    
    ## Examining the power of the test
    pmNeill2.rep <- matrix(NA, 1000, 6)
    colnames(pmNeill2.rep) <- as.character(c(1, 2, 3, 5, 10, 20))

    pmNeill2.rep[, 1] <- sim.neill.test(1000, 10, 1, seed=2008421, true = "weibull", method = "percentiles")
    pmNeill2.rep[, 2] <- sim.neill.test(1000, 10, 2, seed=2008422, true = "weibull", method = "percentiles")
    pmNeill2.rep[, 3] <- sim.neill.test(1000, 10, 3, seed=2008423, true = "weibull", method = "percentiles")
    pmNeill2.rep[, 4] <- sim.neill.test(1000, 10, 5, seed=2008424, true = "weibull", method = "percentiles")
    pmNeill2.rep[, 5] <- sim.neill.test(1000, 10, 10, seed=2008425, true = "weibull", method = "percentiles")
    pmNeill2.rep[, 6] <- sim.neill.test(1000, 10, 20, seed=2008426, true = "weibull", method = "percentiles")
    pFct <- function(x, alpha = 0.05) {sum(x < alpha)/1000}    
    apply(pmNeill2.rep, 2, pFct)

    ## Examining the level of the test
    pmNeill2t.rep <- matrix(NA, 1000, 6)
    colnames(pmNeill2t.rep) <- as.character(c(1, 2, 3, 5, 10, 20))
    
    pmNeill2t.rep[, 1] <- sim.neill.test(1000, 10, 1, seed=2008431, method = "percentiles")
    pmNeill2t.rep[, 2] <- sim.neill.test(1000, 10, 2, seed=2008432, method = "percentiles")
    pmNeill2t.rep[, 3] <- sim.neill.test(1000, 10, 3, seed=2008433, method = "percentiles")
    pmNeill2t.rep[, 4] <- sim.neill.test(1000, 10, 5, seed=2008434, method = "percentiles")
    pmNeill2t.rep[, 5] <- sim.neill.test(1000, 10, 10, seed=2008435, method = "percentiles")
    pmNeill2t.rep[, 6] <- sim.neill.test(1000, 10, 20, seed=2008436, method = "percentiles")
    apply(pmNeill2t.rep, 2, pFct)    
    
    ## Simulation using pairs grouping
    
    ## Examining the power of the test
    pmNeill3.rep <- matrix(NA, 1000, 6)
    colnames(pmNeill3) <- as.character(c(1, 2, 3, 5, 10, 20))

    pmNeill3.rep[, 1] <- sim.neill.test(1000, 10, 1, seed=2008521, true = "weibull", method = "finest")
    pmNeill3.rep[, 2] <- sim.neill.test(1000, 10, 2, seed=2008522, true = "weibull", method = "finest")
    pmNeill3.rep[, 3] <- sim.neill.test(1000, 10, 3, seed=2008523, true = "weibull", method = "finest")
    pmNeill3.rep[, 4] <- sim.neill.test(1000, 10, 5, seed=2008524, true = "weibull", method = "finest")
    pmNeill3.rep[, 5] <- sim.neill.test(1000, 10, 10, seed=2008525, true = "weibull", method = "finest")
    pmNeill3.rep[, 6] <- sim.neill.test(1000, 10, 20, seed=2008526, true = "weibull", method = "finest")
    pFct <- function(x, alpha = 0.05) {sum(x < alpha)/1000}    
    apply(pmNeill3.rep, 2, pFct)

    ## Examining the level of the test
    pmNeill3t.rep <- matrix(NA, 1000, 6)
    colnames(pmNeill3t.rep) <- as.character(c(1, 2, 3, 5, 10, 20))
    
    pmNeill3t.rep[, 1] <- sim.neill.test(1000, 10, 1, seed=2008531, method = "finest")
    pmNeill3t.rep[, 2] <- sim.neill.test(1000, 10, 2, seed=2008532, method = "finest")
    pmNeill3t.rep[, 3] <- sim.neill.test(1000, 10, 3, seed=2008533, method = "finest")
    pmNeill3t.rep[, 4] <- sim.neill.test(1000, 10, 5, seed=2008534, method = "finest")
    pmNeill3t.rep[, 5] <- sim.neill.test(1000, 10, 10, seed=2008535, method = "finest")
    pmNeill3t.rep[, 6] <- sim.neill.test(1000, 10, 20, seed=2008536, method = "finest")
    apply(pmNeill3t.rep, 2, pFct)    
        

    
    
    ## Simulation using fixed grouping

    ## Examining the level of the test
    pmNeillt2 <- matrix(NA, 1000, 7)
    colnames(pmNeillt2) <- as.character(c(10, 25, 50, 75, 100, 250, 500))
        
    pmNeillt2[, 1] <- sim.neill.test(1000, 10, seed=200704201, 
    breaks=c(-Inf, quantiles(x, c(.2,.4,0.6,0.8)), Inf))  # too small
    pmNeillt2[, 1] <- sim.neill.test(1000, 10, seed=200704202, breaks=c(.2,.4,0.6,0.8))  # too small
    pmNeillt2[, 2] <- sim.neill.test(1000, 25, seed=200704202, breaks=c(.2,.4,0.6,0.8))  # too small
    pmNeillt2[, 3] <- sim.neill.test(1000, 50, seed=200704203, breaks=c(.2,.4,0.6,0.8))        
    pmNeillt2[, 4] <- sim.neill.test(1000, 75, seed=200704204, breaks=c(.2,.4,0.6,0.8))    
    pmNeillt2[, 5] <- sim.neill.test(1000, 100, seed=200704205, breaks=c(.2,.4,0.6,0.8))    
    pmNeillt2[, 6] <- sim.neill.test(1000, 250, seed=200704206, breaks=c(.2,.4,0.6,0.8))        
    apply(pmNeillt2, 2, pFct)
}
