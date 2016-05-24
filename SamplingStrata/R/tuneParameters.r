# ----------------------------------------------------
# Function to choose suitable values of the parameters
# for the function "optimizeStrata"
# Author: Giulio Barcaroli
# Date: 4 January 2012
# ----------------------------------------------------
tuneParameters <- function (noptim, nsampl, frame, errors = errors, strata = strata, 
    cens = NULL, strcens = FALSE, alldomains = FALSE, dom = 1, 
    initialStrata, addStrataFactor, minnumstr, iter, pops, mut_chance, 
    elitism_rate) 
{
	frame <- frame[frame$domainvalue == dom, ]
	colnames(frame) <- toupper(colnames(frame))
	checkInput(errors,strata,frame)
    if (length(initialStrata) != noptim) 
        stop("Number of optimisations does not equal values in initialStrata parameter")
    if (length(addStrataFactor) != noptim) 
        stop("Number of optimisations does not equal values in addStrataFactor parameter")
    if (length(minnumstr) != noptim) 
        stop("Number of optimisations does not equal values in minnumstr parameter")
    if (length(iter) != noptim) 
        stop("Number of optimisations does not equal values in iter parameter")
    if (length(pops) != noptim) 
        stop("Number of optimisations does not equal values in pops parameter")
    if (length(mut_chance) != noptim) 
        stop("Number of optimisations does not equal values in mut_chance parameter")
    if (length(elitism_rate) != noptim) 
        stop("Number of optimisations does not equal values in elitism_rate parameter")
    numY <- length(grep("Y", toupper(colnames(frame))))
    for (i in 1:numY) {
        stmt <- paste("Y", i, " <- sum(frame$Y", i, ")", sep = "")
        eval(parse(text = stmt))
    }
    estim <- array(0, c(noptim, nsampl, numY))
    nsimul <- c(1:noptim)
    nstrata <- rep(0, noptim)
    cost <- rep(0, noptim)
    for (i in 1:numY) {
        stmt <- paste("CV", i, " <- rep(0,noptim)", sep = "")
        eval(parse(text = stmt))
    }
    est <- NULL
    for (i in 1:numY) {
        if (i < numY) 
            stmt <- paste("est <- paste(est,'CV", i, "=CV", i, 
                ",',sep='')", sep = "")
        if (i == numY) 
            stmt <- paste("est <- paste(est,'CV", i, "=CV", i, 
                "',sep='')", sep = "")
        eval(parse(text = stmt))
    }
    stmt <- paste("simula <- list(nsimul=nsimul,nstrati=nstrata,cost=cost,", 
        est, ")", sep = "")
    eval(parse(text = stmt))
    simula <- as.data.frame(simula)
    nopt <- rep(0, (noptim * nsampl))
    nsamp <- rep(0, (noptim * nsampl))
    for (i in 1:numY) {
        stmt <- paste("diff", i, " <- rep(0,noptim*nsampl)", 
            sep = "")
        eval(parse(text = stmt))
    }
    est <- NULL
    for (i in 1:numY) {
        if (i < numY) 
            stmt <- paste("est <- paste(est,'diff", i, "=diff", 
                i, ",',sep='')", sep = "")
        if (i == numY) 
            stmt <- paste("est <- paste(est,'diff", i, "=diff", 
                i, "',sep='')", sep = "")
        eval(parse(text = stmt))
    }
    stmt <- paste("differ <- list(nopt=nopt,nsamp=nsamp,", est, 
        ")", sep = "")
    eval(parse(text = stmt))
    differ <- as.data.frame(differ)
    file.remove("iteration.txt")
    file.remove("tuning.txt")
    for (i in (1:noptim)) {
        solution <- optimizeStrata(errors = errors, strata = strata, 
            cens = NULL, strcens = FALSE, alldomains = FALSE, 
            dom = dom, initialStrata = initialStrata[i], addStrataFactor = addStrataFactor[i], 
            minnumstr = minnumstr[i], iter = iter[i], pops = pops[i], 
            mut_chance = mut_chance[i], elitism_rate = elitism_rate[i], 
            realAllocation = TRUE, highvalue = 1e+08, suggestions = NULL, writeFiles=TRUE)
        file.remove("outstrata.txt")
        eval(parse(text = stmt))
        stmt <- paste("file.remove('strata_dom", dom, "_iter", 
            i, ".txt')", sep = "")
        eval(parse(text = stmt))
        stmt <- paste("file.remove('solution_dom", dom, "_iter", 
            i, ".txt')", sep = "")
        eval(parse(text = stmt))
        stmt <- paste("file.remove('results_dom", dom, "_iter", 
            i, ".txt')", sep = "")
        eval(parse(text = stmt))
        stmt <- paste("file.remove('plot_dom", dom, "_iter", 
            i, ".pdf')", sep = "")
        eval(parse(text = stmt))
        stmt <- paste("file.copy('outstrata", dom, ".txt','outstrata.txt')", 
            sep = "")
        eval(parse(text = stmt))
        stmt <- paste("file.copy('outstrata", dom, ".txt','strata_dom", 
            dom, "_iter", i, ".txt')", sep = "")
        eval(parse(text = stmt))
        stmt <- paste("file.copy('solution", dom, ".txt','solution_dom", 
            dom, "_iter", i, ".txt')", sep = "")
        eval(parse(text = stmt))
        stmt <- paste("file.copy('results", dom, ".txt','results_dom", 
            dom, "_iter", i, ".txt')", sep = "")
        eval(parse(text = stmt))
        stmt <- paste("file.copy('plotdom", dom, ".pdf','plot_dom", 
            dom, "_iter", i, ".pdf')", sep = "")
        eval(parse(text = stmt))
		outstrata <- solution$aggr_strata
        newstrata <- updateStrata(strata, solution)
        framenew <- updateFrame(frame, newstrata)
        simula$nstrati[i] <- nrow(outstrata)
        simula$cost[i] <- sum(ceiling(outstrata$SOLUZ))
        for (j in (1:nsampl)) {
            samp <- selectSample(framenew, outstrata, verbatim=FALSE)
            for (k in 1:numY) {
                stmt <- paste("estim[i,j,k] <- sum(samp$Y", k, 
                  " * samp$WEIGHTS)", sep = "")
                eval(parse(text = stmt))
                differ$nopt[(i - 1) * nsampl + j] <- i
                differ$nsamp[(i - 1) * nsampl + j] <- j
                stmt <- paste("differ$diff", k, "[(i-1)*nsampl+j]<- estim[i,j,k] - Y", 
                  k, " ", sep = "")
                eval(parse(text = stmt))
            }
        }
        for (k in 1:numY) {
            stmt <- paste("simula$CV", k, "[i] <- sd(estim[i,,k]) / mean(estim[i,,k])", 
                sep = "")
            eval(parse(text = stmt))
        }
        sink("Iteration.txt", append = TRUE)
        cat("\nIteration :", i, "\n")
        simula[i, ]
        sink()
    }
    ind <- which(simula$cost == min(simula$cost))
    cat("--------------------")
    cat("Best optimization: run # ", ind)
    cat("Required sample cost: ", simula$cost[ind])
    cat("--------------------")
    cat("Values of parameters")
    cat("--------------------")
    cat("initialStrata: ", initialStrata[ind])
    cat("addStrataFactor: ", addStrataFactor[ind])
    cat("minnumstr: ", minnumstr[ind])
    cat("iter: ", iter[ind])
    cat("pops: ", pops[ind])
    cat("nmut_chance: ", mut_chance[ind])
    cat("nelitism_rate: ", elitism_rate[ind])
    cat("n--------------------")
    res <- paste("results_", dom, ".csv", sep = "")
    write.table(simula, res, row.names = FALSE, col.names = TRUE, 
        quote = FALSE, sep = ",")
    dif <- paste("details_", dom, ".csv", sep = "")
    write.table(differ, dif, row.names = FALSE, col.names = TRUE, 
        quote = FALSE, sep = ",")
		
    stmt <- paste("pdf('strata_costs_", dom, ".pdf', width=14, height=10)", 
        sep = "")
    eval(parse(text = stmt))
	split.screen(c(1, 2))
    screen(1)
    plot(simula$nsimul, simula$nstrati, xlim = c(1, nrow(simula)), 
        , xlab = "Runs", ylab = "Number of strata", 
        bty = "n")
    title("Number of optimal strata / optimization run",cex.main=1.0,font.main=1)
    lines(simula$nsimul, simula$nstrati)
    screen(2)
    plot(simula$nsimul, simula$cost, xlim = c(1, nrow(simula)), , xlab = "Runs", 
        ylab = "Sample cost", bty = "n")
    title("Solution cost / optimization run",cex.main=1.0,font.main=1)
    lines(simula$nsimul, simula$cost)
    close.screen(all.screens = TRUE)
    dev.off()
	
    stmt <- paste("pdf('estimates_", dom, ".pdf', width=14, height=10)", 
        sep = "")
    eval(parse(text = stmt))
	nopt <- differ$nopt
    k <- ceiling(numY/4)
    for (j in 1:k) {
        split.screen(c(2, 2))
        for (i in 1:4) {
            if (i + 4 * (j - 1) <= numY) {
                stmt <- paste("screen(", i, ")", sep = "")
                eval(parse(text = stmt))
                stmt <- paste("boxplot(diff", i + 4 * (j - 1), 
                  "~nopt,data=differ,xlab = 'Runs',col='orange')", 
                  sep = "")
                eval(parse(text = stmt))
                stmt <- paste("mtext(expression(Y", i + 
                  4 * (j - 1), "), side=3, adj=0, cex=1.0, line=1)", 
                  sep = "")
                eval(parse(text = stmt))
				title("Distribution of differences between true value and estimate",
					cex.main= 1.0,font.main=1)
            }
        }
        close.screen(all.screens = TRUE)
		dev.off()
    }
}
