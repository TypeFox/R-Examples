AUC.test <- function(pred1, lab1, pred2, lab2, conf.level = 0.95, paired = FALSE){
    if(missing(pred1) || missing(lab1))
        stop("'pred1' and 'lab1' have to be specified")

    if(is.factor(lab1)){ 
        lab <- paste(levels(factor(lab1)), collapse = " and ")
    }else{
        lab <- "0 and 1"
    }
    lab1 <- as.integer(factor(lab1))-1
    tab1 <- table(lab1)
    auc1 <- AUC(pred1, group = lab1)
    se1 <- .seAUC(auc1, npos = tab1[2], nneg = tab1[1])
    alpha <- 1-conf.level
    cint <- qnorm(1-alpha/2)
    res1 <- c(auc1, se1, auc1 - cint*se1, auc1 + cint*se1)
    names(res1) <- c("AUC", "SE", "low CI", "up CI")

    if(missing(pred2) && missing(lab2)){
        res <- wilcox.test(pred1[lab1 == 0], pred1[lab1 == 1], conf.level = conf.level)
        res$null.value <- 0.5
        names(res$null.value) <- "AUC"
        res$data.name <- lab
        return(list(Variable1 = res1, Test = res))
    }
    if(missing(pred2) || missing(lab2)){
        stop("For comparison of AUCs you have to specify 'pred2' and 'lab2'!")
    }
    
    auc2 <- AUC(pred2, group = lab2)
    lab2 <- as.integer(factor(lab2))-1
    tab2 <- table(lab2)
    se2 <- .seAUC(auc2, npos = tab2[2], nneg = tab2[1])
    res2 <- c(auc2, se2, auc2 - cint*se2, auc2 + cint*se2)
    names(res2) <- c("AUC", "SE", "low CI", "up CI")

    if(paired){
        stop("not yet implemented")
    }else{
        r <- 0
    }
    
    sed <- sqrt(se1^2 + se2^2 - 2*r*se1*se2)
    
    zstat <- (auc1-auc2)/sed
    names(zstat) <- "z"
    pval <- 2*pnorm(abs(zstat), lower.tail = FALSE)
    estimate <- auc1-auc2
    names(estimate) <- "Difference in AUC"
    cint <- estimate + c(-cint, cint)*sed
    attr(cint, "conf.level") <- conf.level
    method <- "Hanley and McNeil test for two AUCs"
    dname <- paste(deparse(substitute(pred1)), "and", deparse(substitute(pred2)))
    nullval <- 0
    names(nullval) <- "difference in AUC"
    
    rval <- list(statistic = zstat, p.value = pval, 
                 conf.int = cint, 
                 estimate = estimate,
                 null.value = nullval, 
                 alternative = "two.sided", 
                 method = method, 
                 data.name = dname)
    class(rval) <- "htest"

    list(Variable1 = res1, Variable2 = res2, Test = rval)
}

.seAUC <- function(AUC, npos, nneg){
    a <- AUC
    q1 <- a/(2-a)
    q2 <- (2*a^2)/(1+a)
    sqrt((a*(1-a)+(npos-1)*(q1-a^2)+(nneg-1)*(q2-a^2))/(npos*nneg))
}
