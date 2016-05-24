`nested.npmanova` <-
function(formula,data,method="euc",permutations=100,warnings=FALSE){
    randomize=function(data,toplev,lowlev){
        newdata <- data
        orig.levs <- levels(droplevels(data[,lowlev]))
        nl <- length(orig.levs)
        new.levs <- orig.levs[sample(nl)]
        for (i in 1:nl) {
            subs1 <- data[, lowlev] == orig.levs[i]
            subs2 <- data[, lowlev] == new.levs[i]
            newtops <- data[subs2,toplev]
            newtops <- newtops[1]
            newtops <- rep(newtops, sum(subs1))
            newdata[subs1,toplev] <- newtops
        }
        return(newdata)
    }
    randomize2=function(data,strata){
        newdata <- data
        orig.levs <- levels(droplevels(data[,strata]))
        nl <- length(orig.levs)
        for (i in 1:nl) {
            subs <- data[, strata] == orig.levs[i]
            nsub <- sum(subs==T)
            subdata <- data[subs,]
            newdata[subs,] <- subdata[sample(nsub),]
        }
        return(newdata)
    }
    ow <- options("warn")
    if (warnings==FALSE) {options(warn=0)}
    formula <- as.formula(formula)
    .BiodiversityR <- new.env()
    environment(formula) <- .BiodiversityR
    if (length(all.vars(formula)) > 3) 
        stop(paste("function only works with one main and one nested factor"))
    x <- eval(as.name((all.vars(formula)[1])))
    if (inherits(x, "dist")) {
        distmatrix <- as.matrix(x)
    }else{
        distmatrix <- as.matrix(vegdist(x, method = method))
    }
    SStot <- sum(distmatrix^2)/(2*nrow(distmatrix))
    cat("Total sum of squares of distance matrix:", SStot, "\n")
    resp <- all.vars(formula)[1]
    toplev <- all.vars(formula)[2]
    lowlev <- all.vars(formula)[3]
    data1 <- data
    assign("data1",data1,envir=.BiodiversityR) 
    adonis1 <- adonis(formula,data1,permutations=2,method=method)
    adonis1 <- data.frame(adonis1$aov.tab)
    anovadat <- adonis1[1:3, c(1:4,6)]
    df1 <- anovadat[1,1]
    df2 <- anovadat[2,1]
    df3 <- nrow(distmatrix)-df1-df2-1
    sstop <- anovadat[1,2]
    sslow <- anovadat[2,2]
    ssres <- anovadat[3,2]
    vartot <- adonis1[4,2]
    Ftop <- anovadat[1,3] <- (sstop/df1)/(sslow/df2)
    Flow <- anovadat[2,3] <- (sslow/df2)/(ssres/df3)
    counter <- 1
    for (i in 1:permutations) {
        data2 <- randomize(data,toplev,lowlev)
        assign("data2", data2, envir=.BiodiversityR)
        adonis2 <- adonis(formula,data=data2,method=method,permutations=2)
        adonis2 <- data.frame(adonis2$aov.tab)
        Frand <- (adonis2[1,2]/df1)/(adonis2[2,2]/df2)
        if (Frand >= Ftop) {counter <- counter+1}
    }    
    signi <- counter/(permutations+1)
    anovadat[1,4] <- anovadat[2,4] <- permutations
    anovadat[1,5] <- signi
    counter <- 1
    for (i in 1:permutations) {
        data2 <- randomize2(data,toplev)
        assign("data2", data2, envir=.BiodiversityR)
        adonis2 <- adonis(formula,data=data2,method=method,permutations=2)
        adonis2 <- data.frame(adonis2$aov.tab)
        Frand <- (adonis2[2,2]/df2)/(adonis2[3,2]/df3)
        if (Frand >= Flow) {counter <- counter+1}
    } 
    remove("data1", envir=.BiodiversityR)
    remove("data2", envir=.BiodiversityR)   
    signi <- counter/(permutations+1)
    anovadat[2,5] <- signi
    colnames(anovadat) <- c("Df", "SumsofSquares", "F", "N.Perm", "Pr(>F)")
    mod <- paste("Nested anova for", lowlev, "nested within", toplev, "\n")
    head <- paste("Total sum of squares for non-parametric manova:", vartot, "\n")
    options(ow)
    structure(anovadat, heading = c(head, mod), Random.seed = NA, 
        class = c("anova.cca", "anova", "data.frame"))
}
