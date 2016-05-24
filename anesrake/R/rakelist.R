rakelist <-
function(inputter, dataframe, caseid, 
    weightvec = NULL, cap = 999999, verbose = FALSE, maxit = 1000, 
    convcrit = 0.01) {
    mat <- dataframe
    if (is.null(weightvec)) {
        weightvec <- rep(1, length(caseid))
    }
    prevec <- weightvec
    if (sum(is.na(weightvec)) > 0) {
        stop("seed weights cannot have missing values, use filter to eliminate missing values or substitute 1 for missing cases")
    }
    if (length(weightvec) != length(caseid)) {
        stop("weight vector does not contain the same number of cases as data frame")
    }
    if (cap <= 1) {
        stop("cap may not be less than or equal to 1")
    }
    if (cap < 1.5) {
        print("cap is very low, the model may take a long time to run")
    }
    diferr <- 9999999
    diferrold <- 99999999999
    g <- 0
    pctstill <- 1 - convcrit
    pop <- 0
    while (diferr < pctstill * diferrold) {
        g <- g + 1
        wvold <- weightvec
        if (verbose == TRUE) {
            print(paste("Raking...Iteration", g))
        }
        for (i in 1:length(inputter)) {
            weightvec <- rakeonvar(eval(parse(text = paste("mat", 
                "$", names(inputter[i]), sep = ""))), inputter[[i]], 
                weightvec)
        }
        q <- 0
        while (range(weightvec)[2] > cap + 1e-04) {
            q <- q + 1
            if (verbose == TRUE) {
                print(paste("Capping...Iteration ", g, ".", q, 
                  sep = ""))
            }
            weightvec <- sapply(weightvec, function(x) if (x > 
                cap) {
                x <- cap
            }
            else {
                x <- x
            }, simplify = TRUE)
            weightvec <- weightvec/mean(weightvec)
        }
        if (g %in% seq(100, 10000, 50)) {
            print(paste(g, "iterations have occurred, convergence may not be possible...still working"))
        }
        diferrold <- diferr
        diferr <- sum(abs(weightvec - wvold))
        if (verbose == TRUE) {
            print(paste("Current iteration changed total weights by", 
                diferr))
        }
        if (g > maxit) {
            print(paste("convergence did not occur in", maxit, 
                "iterations"))
            print("output may not be accurate")
            warning("Raking Algorithm Did Not Converge, Results May Be Highly Inconsistent")
            diferrold <- 0
            pop <- 2
            converge <- paste("No convergence in", maxit, "iterations")
        }
    }
    if (diferr > 0.001) {
        print("raking achieved only partial convergence, please check the results to ensure that sufficient convergence was achieved.")
        print(paste("no improvement was apparent after", g, "iterations"))
        print(paste("current total change in the iteration is:", 
            diferr, "average change per weight is:", diferr/sum(weightvec)))
        warning(paste("Raking algorithm achieved only partial convergence, please check the results to ensure that sufficient convergence was achieved.  Average change in weight per case is", 
            diferr/sum(weightvec)))
        warning("Results are stable, but do not perfectly match population marginals")
        diferrx <- diferr
        pop <- 1
        converge <- "Results are stable, but do not perfectly match population marginals"
    }
    if (pop == 0) {
        print(paste("Raking converged in", g, "iterations"))
        diferrx <- diferr
        converge <- "Complete convergence was achieved"
    }
    out <- list(weightvec = weightvec, caseid = caseid, iterations = g, 
        nonconvergence = diferr, converge = converge, varsused = names(inputter), 
        targets = inputter, dataframe = dataframe, prevec=prevec)
    class(out) <- "anesrakelist"
    out
}

