deltaPlot<-function (data, type = "response", group, focal.name, thr = "norm", 
    purify = FALSE, purType = "IPP1", maxIter = 10, alpha = 0.05, 
    extreme = "constraint", const.range = c(0.001, 0.999), nrAdd = 1, 
    save.output = FALSE, output = c("out", "default")) 
{
    test <- switch(type, response = 1, prop = 2, delta = 3)
    if (is.null(test)) 
        stop("'type' must be either 'response','prop' or 'delta'", 
            call. = FALSE)
    test2 <- switch(extreme, constraint = 1, add = 2)
    if (is.null(test2)) 
        stop("'extreme' must be either 'constraint' or 'add'", 
            call. = FALSE)
    if (test > 1 & extreme == "add") 
        stop("'extreme' cannot be 'add' when 'type' is not 'response'", 
            call. = FALSE)
    if (test > 1 & ncol(data) > 2) 
        stop("'data' should not have more than two columns when 'type' is not 'response'", 
            call. = FALSE)
    internalDelta <- function() {
        if (test == 1) {
            if (is.character(group)) {
                DATA <- data[, colnames(data) != group]
                GROUP <- data[, colnames(data) == group]
            }
            else {
                DATA <- data[, (1:ncol(data)) != group]
                GROUP <- data[, group]
            }
            props <- matrix(NA, ncol(DATA), 2)
            for (i in 1:ncol(DATA)) {
                props[i, 1] <- mean(DATA[GROUP != focal.name, 
                  i], na.rm = TRUE)
                props[i, 2] <- mean(DATA[GROUP == focal.name, 
                  i], na.rm = TRUE)
            }
            adjProps <- adjustExtreme(data = DATA, group = GROUP, 
                focal.name = focal.name, prop = props, method = extreme, 
                const.range = const.range, nrAdd = nrAdd)$adj.prop
        }
        if (test == 2) {
            props <- data
            adjProps <- adjustExtreme(data = NULL, group = NULL, 
                focal.name = NULL, prop = props, method = extreme, 
                const.range = const.range, nrAdd = nrAdd)$adj.prop
        }
        if (test == 3) 
            Deltas <- data
        else Deltas <- 4 * qnorm(1 - adjProps) + 13
        SIG <- cov(Deltas)
        M <- colMeans(Deltas)
        b1 <- (SIG[2, 2] - SIG[1, 1] - sqrt((SIG[2, 2] - SIG[1, 
            1])^2 + 4 * SIG[1, 2]^2))/(2 * SIG[1, 2])
        b2 <- (SIG[2, 2] - SIG[1, 1] + sqrt((SIG[2, 2] - SIG[1, 
            1])^2 + 4 * SIG[1, 2]^2))/(2 * SIG[1, 2])
        b <- max(c(b1, b2))
        a <- M[2] - b * M[1]
        C <- c(b, -1)/sqrt(b^2 + 1)
        epsilon <- a/sqrt(b^2 + 1)
        DIST <- Deltas %*% C + epsilon
if (sum(is.na(DIST))>0) 
stop("Perpendicular distances cannot be computed - one set of Delta scores is probably constant",call.=FALSE)

        if (thr == "norm") {
            mat <- t(C) %*% SIG %*% C
            Q <- qnorm(1 - alpha/2, 0, sqrt(mat))
            rule <- "norm"
        }
        else {
            Q <- abs(thr)
            rule <- "fixed"
        }
        if (max(abs(DIST)) <= Q) 
            DIFitems <- "no DIF item detected"
        else DIFitems <- (1:nrow(Deltas))[abs(DIST) > Q]
        if (test < 3) {
            stats <- props
            stats2 <- adjProps
        }
        else stats <- stats2 <- NA
        if (!purify) 
            RES <- list(Props = stats, adjProps = stats2, Deltas = Deltas, 
                Dist = DIST, axis.par = c(a, b), thr = Q, rule = rule, 
                DIFitems = DIFitems, adjust.extreme = extreme, 
                const.range = const.range, nrAdd = nrAdd, purify = purify, 
                alpha = alpha, save.output = save.output, output = output)
        else {
            if (is.character(DIFitems)) 
                RES <- list(Props = stats, adjProps = stats2, 
                  Deltas = Deltas, Dist = DIST, axis.par = c(a, 
                    b), nrIter = 1, maxIter = maxIter, convergence = TRUE, 
                  thr = Q, rule = rule, purType = purType, DIFitems = DIFitems, 
                  adjust.extreme = extreme, const.range = const.range, 
                  nrAdd = nrAdd, purify = purify, alpha = alpha, 
                  save.output = save.output, output = output)
            else {
                convergence <- FALSE
                iter <- 1
                dif <- DIFitems
                difPur <- rep(0, nrow(Deltas))
                difPur[DIFitems] <- 1
                Qseries <- Q
                pars <- c(a, b)
                allDist <- DIST
                repeat {
                  iter <- iter + 1
                  if (iter > maxIter) 
                    break
                  else {
                    nodif <- NULL
                    for (i in 1:nrow(Deltas)) {
                      if (sum(i == dif) == 0) 
                        nodif <- c(nodif, i)
                    }
                    DeltaProv <- cbind(Deltas[nodif, 1], Deltas[nodif, 
                      2])
                    SSIG <- cov(DeltaProv)
                    MM <- colMeans(DeltaProv)
                    bb1 <- (SSIG[2, 2] - SSIG[1, 1] - sqrt((SSIG[2, 
                      2] - SSIG[1, 1])^2 + 4 * SSIG[1, 2]^2))/(2 * 
                      SSIG[1, 2])
                    bb2 <- (SSIG[2, 2] - SSIG[1, 1] + sqrt((SSIG[2, 
                      2] - SSIG[1, 1])^2 + 4 * SSIG[1, 2]^2))/(2 * 
                      SSIG[1, 2])
                    bb <- max(c(bb1, bb2))
                    aa <- MM[2] - bb * MM[1]
                    C <- c(bb, -1)/sqrt(bb^2 + 1)
                    epsilon <- aa/sqrt(bb^2 + 1)
                    DIST2 <- Deltas %*% C + epsilon
                    if (thr == "norm") {
                      if (purType == "IPP1") {
                        Q <- Qseries[1]
                        rule <- "norm"
                      }
                      else {
                        if (purType == "IPP3") {
                          mat <- t(C) %*% SSIG %*% C
                          Q <- qnorm(1 - alpha/2, 0, sqrt(mat))
                          rule <- "norm"
                        }
                        else {
                          mat <- t(C) %*% SIG %*% C
                          Q <- qnorm(1 - alpha/2, 0, sqrt(mat))
                          rule <- "norm"
                        }
                      }
                    }
                    else {
                      Q <- abs(thr)
                      rule <- "fixed"
                      purType <- "IPP1"
                    }
                    if (max(abs(DIST2)) <= Q) 
                      dif2 <- "no DIF item detected"
                    else dif2 <- (1:nrow(Deltas))[abs(DIST2) > 
                      Q]
                    difPur <- rbind(difPur, rep(0, nrow(Deltas)))
                    if (!is.character(dif2)) 
                      difPur[iter, dif2] <- 1
                    Qseries <- c(Qseries, Q)
                    pars <- rbind(pars, c(aa, bb))
                    allDist <- cbind(allDist, DIST2)
                    if (length(dif) != length(dif2)) 
                      dif <- dif2
                    else {
                      if (sum(abs(difPur[iter, ] - difPur[iter - 
                        1, ])) == 0) {
                        convergence <- TRUE
                        dif <- dif2
                        break
                      }
                      else dif <- dif2
                    }
                  }
                }
                RES <- list(Props = stats, adjProps = stats2, 
                  Deltas = Deltas, Dist = allDist, axis.par = pars, 
                  nrIter = nrow(difPur), maxIter = maxIter, convergence = convergence, 
                  difPur = difPur, thr = Qseries, rule = rule, 
                  purType = purType, DIFitems = dif, adjust.extreme = extreme, 
                  const.range = const.range, nrAdd = nrAdd, purify = purify, 
                  alpha = alpha, save.output = save.output, output = output)
            }
        }
        class(RES) <- "deltaPlot"
        return(RES)
    }
    resToReturn <- internalDelta()
    if (save.output) {
        if (output[2] == "default") 
            wd <- file.path(getwd())
        else wd <- output[2]
nameF<-paste(output[1],".txt",sep="")
        fileName <- file.path(wd, nameF)
        capture.output(resToReturn, file = fileName)
    }
    return(resToReturn)
}





## PRINT FUNCTION

print.deltaPlot<-function (x, only.final = TRUE, ...) 
{
    res <- x
    cat("\n")
    cat("Detection of Differential Item Functioning using Angoff's Delta method", 
        "\n")
    if (res$purify) 
        cat("  with item purification", "\n", "\n")
    else cat("  without item purification", "\n", "\n")
    if (res$purify) {
        if (res$convergence) {
            if (res$nrIter == 1) 
                cat("Convergence reached after", res$nrIter, 
                  "iteration", "\n", "\n")
            else cat("Convergence reached after", res$nrIter, 
                "iterations", "\n", "\n")
        }
        else {
            cat("WARNING: convergence was not reached after", 
                res$maxIter, "iterations!", "\n", "\n")
        }
        if (res$nrIter > 1) {
            if (res$purType == "IPP1") {
                cat("Threshold kept fixed to", res$thr[1], "\n")
                if (res$rule == "fixed") 
                  cat(" (as fixed by the user [IPP1])", "\n", 
                    "\n")
                else cat(" (as computed from normal approximation [IPP1])", 
                  "\n", "\n")
            }
            else {
                cat("Threshold adjusted iteratively using normal approximation", 
                  "\n")
                cat(" and ", round(res$alpha * 100), "% significance level", 
                  "\n", sep = "")
                if (res$purType == "IPP2") 
                  cat(" (only slope parameter updated [IPP2])", 
                    "\n", "\n")
                else cat(" (full update of the threshold [IPP3])", 
                  "\n", "\n")
            }
        }
    }
    if (res$adjust.extreme == "constraint") 
        cat("Extreme proportions adjusted by constraining to [", 
            round(res$const.range[1], 3), "; ", res$const.range[2], 
            "]", "\n", "\n", sep = "")
    else {
        if (res$nrAdd == 1) 
            cat("Extreme proportions adjusted by adding one success and one failure", 
                "\n", "\n")
        else cat("Extreme proportions adjusted by adding ", res$nrAdd, 
            " successes and ", res$nrAdd, " failures", "\n", 
            "\n", sep = "")
    }
    if (res$purify) 
        cat("Statistics (after the first iteration):", "\n", 
            "\n")
    else cat("Statistics:", "\n", "\n")
    if (is.null(dim(res$Props))) {
        m1 <- round(cbind(res$Deltas, res$Dist[, 1]), 4)
        symb <- symnum(abs(as.numeric(res$Dist[, 1])), c(0, abs(res$thr[1]), 
            Inf), symbols = c("", "***"))
        m1 <- noquote(cbind(format(m1, justify = "right"), symb))
        colnames(m1) <- c("Delta.Ref", "Delta.Foc", "Dist.", 
            "")
    }
    else {
        m1 <- round(cbind(res$Props, res$Deltas, res$Dist[, 1]), 
            4)
        symb <- symnum(abs(as.numeric(res$Dist[, 1])), c(0, abs(res$thr[length(res$thr)]), 
            Inf), symbols = c("", "***"))
        m1 <- noquote(cbind(format(m1, justify = "right"), symb))
        colnames(m1) <- c("Prop.Ref", "Prop.Foc", "Delta.Ref", 
            "Delta.Foc", "Dist.", "")
    }
    rn <- NULL
    for (i in 1:nrow(m1)) rn[i] <- paste("Item", i, sep = "")
    rownames(m1) <- rn
    print(m1)
    cat("\n")
    cat("Code: '***' if item is flagged as DIF", "\n", "\n")
    if (res$purify) {
        cat("Statistics (after the last iteration):", "\n", "\n")
        if (is.null(dim(res$Props))) {
            m1 <- round(cbind(res$Deltas, res$Dist[, ncol(res$Dist)]), 
                4)
            symb <- symnum(abs(as.numeric(res$Dist[, ncol(res$Dist)])), 
                c(0, abs(res$thr[length(res$thr)]), Inf), symbols = c("", 
                  "***"))
            m1 <- noquote(cbind(format(m1, justify = "right"), 
                symb))
            colnames(m1) <- c("Delta.Ref", "Delta.Foc", "Dist.", 
                "")
            rn <- NULL
            for (i in 1:nrow(m1)) rn[i] <- paste("Item", i, sep = "")
            rownames(m1) <- rn
        }
        else {
            m1 <- round(cbind(res$Props, res$Deltas, res$Dist[, 
                ncol(res$Dist)]), 4)
            symb <- symnum(abs(as.numeric(res$Dist[, ncol(res$Dist)])), 
                c(0, abs(res$thr[length(res$thr)]), Inf), symbols = c("", 
                  "***"))
            m1 <- noquote(cbind(format(m1, justify = "right"), 
                symb))
            colnames(m1) <- c("Prop.Ref", "Prop.Foc", "Delta.Ref", 
                "Delta.Foc", "Dist.", "")
            rn <- NULL
            for (i in 1:nrow(m1)) rn[i] <- paste("Item", i, sep = "")
            rownames(m1) <- rn
        }
        print(m1)
        cat("\n")
        cat("Code: '***' if item is flagged as DIF", "\n", "\n")
    }
    if (!only.final) {
        cat("Perpendicular distances:", "\n", "\n")
        m1 <- round(res$Dist, 4)
        rc <- NULL
        for (t in 1:ncol(res$Dist)) rc[t] <- paste("Iter", t, 
            sep = "")
        colnames(m1) <- rc
        rn <- NULL
        for (i in 1:nrow(m1)) rn[i] <- paste("Item", i, sep = "")
        rownames(m1) <- rn
        print(m1)
        cat("\n")
    }
myBool<-ifelse(!res$purify,TRUE,ifelse(res$nrIter==1,TRUE,FALSE))
    if (myBool) {
        cat("Parameters of the major axis:", "\n", "\n")
        np <- round(rbind(res$axis.par), 4)
        rownames(np) <- ""
        colnames(np) <- c("a", "b")
        print(np)
        cat("\n")
        if (res$rule == "norm") 
            cat("Detection threshold: ", round(res$thr, 4), " (significance level: ", 
                round(res$alpha * 100, 0), "%)", sep = "", "\n", 
                "\n")
        else cat("Detection threshold: ", round(res$thr, 4), 
            sep = "", "\n", "\n")
    }
    else {
        if (only.final) {
            cat("Parameters of the major axis (first and last iterations only):", 
                "\n", "\n")
            if (is.null(dim(res$axis.par))) 
                np <- round(rbind(res$axis.par, res$axis.par), 
                  4)
            else np <- round(rbind(res$axis.par[c(1, nrow(res$axis.par)), 
                ]), 4)
            rownames(np) <- c("First", "Last")
            colnames(np) <- c("a", "b")
            print(np)
            cat("\n")
            if (res$rule == "norm") {
                cat("First and last detection thresholds: ", 
                  round(res$thr[1], 4), " and ", round(res$thr[length(res$thr)], 
                    4), sep = "", "\n")
                cat(" (significance level: ", round(res$alpha * 
                  100, 0), "%)", sep = "", "\n", "\n")
            }
            else cat("First and last detection thresholds: ", 
                round(res$thr[1], 4), " and ", round(res$thr[length(res$thr)], 
                  4), sep = "", "\n")
        }
        else {
            cat("Parameters of the major axis:", "\n", "\n")
            np <- round(rbind(res$axis.par), 4)
            npr <- NULL
            for (i in 1:nrow(res$axis.par)) npr[i] <- paste("Iter", 
                i, sep = "")
            rownames(np) <- npr
            colnames(np) <- c("a", "b")
            print(np)
            cat("\n")
            cat("Detection thresholds:", "\n", "\n")
            mm <- rbind(res$thr)
            rownames(mm) <- ""
            cn <- NULL
            for (i in 1:length(res$thr)) cn[i] <- paste("Iter", 
                i, sep = "")
            colnames(mm) <- cn
            print(mm)
            cat("\n")
            if (res$rule == "norm") 
                cat("(significance level: ", round(res$alpha * 
                  100, 0), "%)", sep = "", "\n", "\n")
            else cat("\n")
        }
    }
    if (is.character(res$DIFitems)) 
        cat("Items detected as DIF items:", res$DIFitems, "\n", 
            "\n")
    else {
        cat("Items detected as DIF items:", "\n")
        namedif <- NULL
        for (i in 1:length(res$DIFitems)) namedif[i] <- paste("Item", 
            res$DIFitems[i], sep = "")
        m3 <- cbind(namedif)
        rownames(m3) <- rep("", length(res$DIFitems))
        colnames(m3) <- ""
        print(m3, quote = FALSE)
        cat("\n")
    }
    if (!res$save.output) 
        cat("Output was not captured!", "\n")
    else {
        if (res$output[2] == "default") 
            wd <- file.path(getwd())
        else wd <- res$output[2]
nameF<-paste(res$output[1], ".txt", sep = "")
        fileName <- file.path(wd, nameF)
        cat("Output was captured and saved into file", "\n", 
            " '", fileName, "'", "\n", "\n", sep = "")
    }
}

