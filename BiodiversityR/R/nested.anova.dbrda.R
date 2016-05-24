`nested.anova.dbrda` <-
function (formula, data, method = "euc", add = FALSE, permutations = 100, 
    warnings = FALSE) 
{
    randomize = function(data, toplev, lowlev) {
        newdata <- data
        orig.levs <- levels(droplevels(data[, lowlev]))
        nl <- length(orig.levs)
        new.levs <- orig.levs[sample(nl)]
        for (i in 1:nl) {
            subs1 <- data[, lowlev] == orig.levs[i]
            subs2 <- data[, lowlev] == new.levs[i]
            newtops <- data[subs2, toplev]
            newtops <- newtops[1]
            newtops <- rep(newtops, sum(subs1))
            newdata[subs1, toplev] <- newtops
        }
        return(newdata)
    }
    randomize2 = function(data, strata) {
        newdata <- data
        orig.levs <- levels(droplevels(data[, strata]))
        nl <- length(orig.levs)
        for (i in 1:nl) {
            subs <- data[, strata] == orig.levs[i]
            nsub <- sum(subs == T)
            subdata <- data[subs, ]
            newdata[subs, ] <- subdata[sample(nsub), ]
        }
        return(newdata)
    }
    ow <- options("warn")
    if (warnings == FALSE) {
        options(warn = 0)
    }
    formula <- as.formula(formula)
    if (length(all.vars(formula)) > 3) 
        stop(paste("function only works with one main and one nested factor"))
    x <- eval(as.name((all.vars(formula)[1])))
    if (inherits(x, "dist")) {
        distmatrix <- as.matrix(x)
    }
    else {
        distmatrix <- as.matrix(vegdist(x, method = method))
    }
    SStot <- sum(distmatrix^2)/(2 * nrow(distmatrix))
    cat("Total sum of squares of distance matrix:", SStot, "\n")
    resp <- all.vars(formula)[1]
    toplev <- all.vars(formula)[2]
    lowlev <- all.vars(formula)[3]
    .BiodiversityR <- new.env()
    environment(formula) <- .BiodiversityR
    data1 <- data
    assign("data1", data, envir=.BiodiversityR)
    assign("data1", data1, envir=.BiodiversityR)
    METHODS <- c("manhattan", "euclidean", "canberra", "bray", 
        "kulczynski", "gower", "morisita", "horn", "mountford", 
        "jaccard", "raup", "binomial", "chao")
    methodid <- pmatch(method, METHODS)
    method <- METHODS[methodid]
    if (add == TRUE) {
        if (method == "manhattan") {
            model <- capscale(formula, data1, distance = "manhattan", 
                add = TRUE)
        }
        if (method == "euclidean") {
            model <- capscale(formula, data1, distance = "euclidean")
        }
        if (method == "canberra") {
            model <- capscale(formula, data1, distance = "canberra", 
                add = TRUE)
        }
        if (method == "bray") {
            model <- capscale(formula, data1, distance = "bray", 
                add = TRUE)
        }
        if (method == "kulczynski") {
            model <- capscale(formula, data1, distance = "kulczynski", 
                add = TRUE)
        }
        if (method == "gower") {
            model <- capscale(formula, data1, distance = "gower", 
                add = TRUE)
        }
        if (method == "morisita") {
            model <- capscale(formula, data1, distance = "morisita", 
                add = TRUE)
        }
        if (method == "horn") {
            model <- capscale(formula, data1, distance = "horn", 
                add = TRUE)
        }
        if (method == "mountford") {
            model <- capscale(formula, data1, distance = "mountford", 
                add = TRUE)
        }
        if (method == "jaccard") {
            model <- capscale(formula, data1, distance = "jaccard", 
                add = TRUE)
        }
        if (method == "raup") {
            model <- capscale(formula, data1, distance = "raup", 
                add = TRUE)
        }
        if (method == "binomial") {
            model <- capscale(formula, data1, distance = "binomial", 
                add = TRUE)
        }
        if (method == "chao") {
            model <- capscale(formula, data1, distance = "chao", 
                add = TRUE)
        }
    }
    else {
        if (method == "manhattan") {
            model <- capscale(formula, data1, distance = "manhattan", 
                add = FALSE)
        }
        if (method == "euclidean") {
            model <- capscale(formula, data1, distance = "euclidean")
        }
        if (method == "canberra") {
            model <- capscale(formula, data1, distance = "canberra", 
                add = FALSE)
        }
        if (method == "bray") {
            model <- capscale(formula, data1, distance = "bray", 
                add = FALSE)
        }
        if (method == "kulczynski") {
            model <- capscale(formula, data1, distance = "kulczynski", 
                add = FALSE)
        }
        if (method == "gower") {
            model <- capscale(formula, data1, distance = "gower", 
                add = FALSE)
        }
        if (method == "morisita") {
            model <- capscale(formula, data1, distance = "morisita", 
                add = FALSE)
        }
        if (method == "horn") {
            model <- capscale(formula, data1, distance = "horn", 
                add = FALSE)
        }
        if (method == "mountford") {
            model <- capscale(formula, data1, distance = "mountford", 
                add = FALSE)
        }
        if (method == "jaccard") {
            model <- capscale(formula, data1, distance = "jaccard", 
                add = FALSE)
        }
        if (method == "raup") {
            model <- capscale(formula, data1, distance = "raup", 
                add = FALSE)
        }
        if (method == "binomial") {
            model <- capscale(formula, data1, distance = "binomial", 
                add = FALSE)
        }
        if (method == "chao") {
            model <- capscale(formula, data1, distance = "chao", 
                add = FALSE)
        }
    }
# remember the data
    model$call$data <- data1
#
    anovares <- anova(model, perm = 2, by="terms")
    anovadat <- data.frame(anovares)
    adjust <- nrow(model$CCA$u) - 1
    if (pmatch("mean", model$inertia, nomatch = -1) > 0) {
        anovadat[, 2] <- anovadat[, 2] * adjust
        model$tot.chi <- model$tot.chi * adjust
    }else{
        anovadat[, 2] <- anovadat[, 2] / adjust
        model$tot.chi <- model$tot.chi / adjust
    }
    df1 <- anovadat[1, 1]
    df2 <- anovadat[2, 1]
    anovadat[3, 1] <- df3 <- nrow(distmatrix) - df1 - df2 - 1
    anovadat[2, 3] <- anovadat[2, 2]/df3
    formula1 <- as.formula(paste(resp, "~", lowlev, "+Condition(", 
        toplev, ")"))
    environment(formula1) <- .BiodiversityR
    model1 <- capscale(formula1, data = data1, distance = method, 
        add = add)
    Ftop <- (model1$pCCA$tot.chi/df1)/(model1$CCA$tot.chi/df2)
    anovadat[1, 3] <- Ftop
    counter <- 1
    for (i in 1:permutations) {
        data2 <- randomize(data, toplev, lowlev)
        assign("data2", data2, envir=.BiodiversityR)
        Ordinationperm <- capscale(formula1, data = data2, 
            distance = method, add = add)
        randomF <- (Ordinationperm$pCCA$tot.chi/df1)/(Ordinationperm$CCA$tot.chi/df2)
        if (randomF >= Ftop) {
            counter <- counter + 1
        }
    }
    signi <- counter/(permutations + 1)
    anovadat[1, 4] <- anovadat[2, 4] <- permutations
    anovadat[1, 5] <- signi
    Flow <- (model1$CCA$tot.chi/df2)/(model1$CA$tot.chi/df3)
    anovadat[2, 3] <- Flow
    counter <- 1
    for (i in 1:permutations) {
        data2 <- randomize2(data, toplev)
        assign("data2", data2, envir=.BiodiversityR)
        Ordinationperm <- capscale(formula1, data = data2, 
            distance = method, add = add)
        randomF <- (Ordinationperm$CCA$tot.chi/df2)/(Ordinationperm$CA$tot.chi/df3)
        if (randomF >= Flow) {
            counter <- counter + 1
        }
    }
    remove("data1", envir=.BiodiversityR)
    remove("data2", envir=.BiodiversityR)
    signi <- counter/(permutations + 1)
    anovadat[2, 5] <- signi
    colnames(anovadat) <- c("Df", "SumsOfSquares", "F", "N.Perm", "Pr(>F)")
    mod <- paste("Nested anova for", lowlev, "nested within", 
        toplev, "\n")
    head <- paste("Total sum of squares of distance-based redundancy analysis:", 
        c(model$tot.chi), "\n")
    options(ow)
    structure(anovadat, heading = c(head, mod), Random.seed = NA, 
        class = c("anova.cca", "anova", "data.frame"))
}

# test
# model.test <- nested.anova.dbrda(warcom~rift.valley+popshort, data=warenv, method="jac", permutations=5)


