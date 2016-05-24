startingValues <-
function (phen, genZ, reference = "noia", max.level = 2, max.dom = 2, 
    fast = FALSE, e.unique = FALSE, start.algo = "linear", bilinear.steps = NULL, 
    ...) 
{
    if (is.null(start.algo)) {
        start.algo <- "nostart"
    }
    if (start.algo == "nostart") {
        start.values <- startingValuesReg(NULL, max.level, max.dom, 
            e.unique, nloc = (ncol(genZ)/3))
    }
    else if (start.algo == "linear") {
        linear <- linearRegression(phen = phen, genZ = genZ, 
            reference = reference, max.level = max.level, max.dom = max.dom, 
            fast = fast)
        start.values <- startingValuesReg(linear, max.level, 
            max.dom, e.unique)
    }
    else if (start.algo == "multilinear") {
        mlinear <- multilinearRegression(phen = phen, genZ = genZ, 
            reference = reference, max.level = 2, max.dom = 0, 
            fast = fast, e.unique = e.unique, start.algo = "linear", 
            start.values = NULL, ...)
        start.values <- startingValuesReg(mlinear, max.level, 
            max.dom, e.unique)
    }
    else if (start.algo == "subset") {
        subset.size <- as.integer(length(phen)/2)
        subsample <- sample(1:length(phen), subset.size)
        mlinear <- multilinearRegression(phen = phen[subsample], 
            genZ = genZ[subsample, ], reference = reference, 
            max.level = max.level, max.dom = max.dom, fast = fast, 
            e.unique = e.unique, start.algo = "linear", ...)
        start.values <- startingValuesReg(mlinear, max.level, 
            max.dom, e.unique)
    }
    else if (start.algo == "bilinear") {
        linear <- linearRegression(phen = phen, genZ = genZ, 
            reference = reference, max.level = max.level, max.dom = max.dom, 
            fast = fast)
        mdom <- 1
        if (max.dom < 1) 
            mdom <- max.dom
        st1 <- startingValuesReg(linear, max.level = 1, max.dom = mdom, 
            e.unique = e.unique)
        mdom <- 2
        if (max.dom < 2) 
            mdom <- max.dom
        mlev <- 2
        if (max.level < 2) 
            mlev <- max.level
        st2 <- startingValuesReg(linear, max.level = mlev, max.dom = mdom, 
            e.unique = e.unique)
        st2[names(st1)] <- NULL
        form <- formulaMultilinear(nloc = linear$nloc, max.level = mlev, 
            max.dom = mdom, e.unique = e.unique)
        X <- as.data.frame(linear$x)
        phen <- linear$phen
        temp <- list(marginal = st1, interactions = st2)
        st <- 1
        repeat {
            temp2 <- bilinearStep(form, X, phen, temp$marginal, 
                temp$interactions, ...)
            nvar.left <- sum(abs(1 - (unlist(temp)/unlist(temp2))) > 
                nls.control()$tol)
            conv.index <- sum(abs(1 - (unlist(temp)/unlist(temp2))))
            if (bilinear.steps > 2) 
                cat(nvar.left, "variables remaining; convergence index:", 
                  conv.index, "\n")
            if (nvar.left == 0) {
                temp <- temp2
                break
            }
            temp <- temp2
            st <- st + 1
            if (!is.null(bilinear.steps) && st > bilinear.steps) {
                break
            }
        }
        start.values <- as.list(merge(temp$marginal, temp$interactions))
    }
    else {
        stop("Starting values algorithm", start.algo, "is unknown")
    }
    return(start.values)
}
