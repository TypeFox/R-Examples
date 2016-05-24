###################################################################
##                                                                ##
##                    Pst Function                                ##
##             Transition Probability P(s,t)                      ##
##                                                                ##
####################################################################

## Function for P(s,t)
## Takes msSurv object, 1st (s), & last time (t)

Pst <- function (object, s = 0, t = "last", deci = 4, covar = FALSE) {
    if (!(0 <= s & s < t))
        stop("'s' and 't' must be positive, and s < t")
    if (t < et(object)[1] | s >= et(object)[length(et(object))])
        stop("Either 's' or 't' is an invalid time")
    if (t == "last")
        t <- et(object)[length(et(object))]
    ## NOTE - Changed s <= to s <
    idx <- which(s < et(object) & et(object) <= t)
    l.idx <- length(idx)
    cum.prod <- diag(ns(object))
    rownames(cum.prod) <- nodes(tree(object))
    red.AJs <- array(dim = c(ns(object), ns(object), nrow(dNs(object))),
        dimnames = list(rows = nodes(tree(object)), cols = nodes(tree(object)),
            dim = rownames(dNs(object))))
    for (i in idx) {
        cum.prod <- cum.prod %*% I.dA(object)[, , i]
        red.AJs[, , i] <- cum.prod
    }
    if (covar == TRUE) {
        bl.Id <- diag(1, (ns(object))^2)
        cov.Pst <- array(0, dim = c(dim(I.dA(object)[, , idx])[c(1,
            2)]^2, nrow(dNs(object))))
        colnames(cov.Pst) <- rownames(cov.Pst) <- paste(rep(nodes(tree(object)),
            ns(object)), sort(rep(nodes(tree(object)), ns(object))))
        Id <- diag(1, ns(object))
        for (i in idx) {
            if (i == idx[1])
                cov.Pst[, , i] <- bl.Id %*% cov.dA(object)[,
                  , i] %*% bl.Id
            else cov.Pst[, , i] <- (t(I.dA(object)[, , i]) %x%
                Id) %*% cov.Pst[, , i - 1] %*% ((I.dA(object)[,
                , i]) %x% Id) + (Id %x% red.AJs[, , i - 1]) %*%
                cov.dA(object)[, , i] %*% (Id %x% t(red.AJs[,
                , i - 1]))
        }
    }
    cat(paste("Estimate of P(", s, ", ", t, ")\n", sep = ""))
    print(round(cum.prod, digits = deci))
    cat("\n")
    if (!is.null(cov.AJs(object)) & covar == TRUE) {
        cat(paste("Estimate of cov(P(", s, ", ", t, "))\n", sep = ""))
        print(round(cov.Pst[, , max(idx)], digits = deci))
    }
    if (!is.null(cov.AJs(object)) & covar == TRUE) {
        return(invisible(list(Pst = cum.prod, cov.Pst = cov.Pst[, , max(idx)])))
    }
    else {
        return(invisible(list(Pst = cum.prod)))
    }
}


####################################################################
##                                                                ##
##                       SOPt Function                            ##
##                                                                ##
##                State Occ for Specific Time t                   ##
##                                                                ##
####################################################################


SOPt <- function(object, t="last", deci=4, covar=FALSE) {

    if (t=="last") t <- et(object)[length(et(object))]
    t.loc <- length(et(object)[et(object)<= t])

    cat(paste("The state occupation probabilities at time ", t, " are:\n", sep = ""))
    for (node in nodes(tree(object))) {
        p.node <- paste("p", node)
        cat(paste("State ", node, ": ", round(ps(object)[t.loc, p.node], deci), "\n", sep = ""))
    }
    cat("\n")

    if (!is.null(cov.AJs(object)) & covar == TRUE) {
        cat(paste("Covariance estimates for state occupation probabilities: \n", sep = ""))

        for (node in nodes(tree(object))) {
            p.node <- paste("Var p", node)
            cat(paste("State ", node, ": ", round(var.sop(object)[t.loc, p.node], deci), "\n", sep = ""))
        }
    }

    if (!is.null(var.sop(object)) & covar == TRUE) {
        return(invisible(list(SOPt=ps(object)[t.loc, ], var.sop=var.sop(object)[t.loc,])))
    } else {
        return(invisible(list(SOPt=ps(object)[t.loc, ])))
    }
}


####################################################################
##                                                                ##
##                      EntryExit Function                        ##
##                                                                ##
##                State Entry/Exit Time Distribution              ##
##           Calculates est of Fs & Gs at specified time t        ##
##              Also gives BS Variance Estimate & CIs             ##
##                                                                ##
####################################################################


## NOTE - Changed covar argument to var ...
EntryExit <- function(object, t="last", deci=4, var=FALSE, norm=TRUE) {

    if (var==TRUE & is.null(Fnorm.var(object))) {
        stop(paste("msSurv object does not have variance estimates for entry/exit time distributions.\n
Please re-run the msSurv function with the argument 'bs=TRUE'. \n", sep=""))
    }

    entry.st <- sapply(strsplit(colnames(Fsub(object)), " "), function(x) x[2])
    exit.st <- sapply(strsplit(colnames(Gsub(object)), " "), function(x) x[2])
    no.entry <- nodes(tree(object))[!nodes(tree(object)) %in% entry.st]
    no.exit <- nodes(tree(object))[!nodes(tree(object)) %in% exit.st]

    ## entry.st <- names(which(!(sapply(inEdges(tree(object)), function(x) length(x) == 0))))
    ## initial <- setdiff(nodes(tree(object)), entry.st) ## initial states, no Fs
    ## exit.st <- names(which((sapply(edges(tree(object)), function(x) length(x) > 0))))
    ## terminal <- setdiff(nodes(tree(object)), exit.st)

    if (t=="last") t <- et(object)[length(et(object))]
    t.loc <- length(et(object)[et(object)<= t])

    if (norm) {
	cat("The normalized state entry distributions at time ", t, " are:\n")
	for (state in entry.st) {
            f.state <- paste("F", state)
            cat("State ", state, ": ", round(Fnorm(object)[t.loc, f.state], deci), "\n")
	}
        cat("Normalized state entry distributions not estimated for state(s)", no.entry)
	cat("\n\n")

        if (var==TRUE) {
            cat("Variance estimates for normalized state entry distributions: \n")
            for (state in entry.st) {
                f.state <- paste("F", state)
                cat("State ", state, ": ", round(Fnorm.var(object)[t.loc, f.state], deci), "\n")
            }
        }

        cat("\n")

	cat("The normalized state exit distributions at time ", t, " are:\n")
	for (state in exit.st) {
            g.state <- paste("G", state)
            cat("State ", state, ": ", round(Gnorm(object)[t.loc, g.state], deci), "\n")
	}
	cat("Normalized state exit distributions not estimated for state(s)", no.exit)
	cat("\n\n")


        if (var==TRUE) {
            cat("Variance estimates for normalized state exit distributions: \n")
            for (state in exit.st) {
                g.state <- paste("G", state)
                cat("State ", state, ": ", round(Gnorm.var(object)[t.loc, g.state], deci), "\n")
            }
            cat("\n")
        }

        if (var==TRUE) {
            return(invisible(list(entry.norm = Fnorm(object)[t.loc, ],
                                  exit.norm = Gnorm(object)[t.loc, ],
                                  entry.var.norm = Fnorm.var(object)[t.loc, ],
                                  exit.var.norm = Gnorm.var(object)[t.loc, ])))
        } else {
            return(invisible(list(entry.norm = Fnorm(object)[t.loc, ],
                                  exit.norm = Gnorm(object)[t.loc, ])))
        }

    } ## End of Normalized Distrubutions


    if (norm==FALSE) {
	cat("The state entry subdistributions at time ", t, " are:\n")
	for (state in entry.st) {
            f.state <- paste("F", state)
            cat("State ", state, ": ", round(Fsub(object)[t.loc, f.state], deci), "\n")
	}
        cat("State entry subdistributions not estimated for state(s)", no.entry)
	cat("\n\n")

        if (var==TRUE) {
            cat("Variance estimates for normalized state entry distributions: \n")
            for (state in entry.st) {
                f.state <- paste("F", state)
                cat("State ", state, ": ", round(Fsub.var(object)[t.loc, f.state], deci), "\n")
            }
        }

        cat("\n")

	cat("The state exit subdistributions at time ", t, " are:\n")
	for (state in exit.st) {
            g.state <- paste("G", state)
            cat("State ", state, ": ", round(Gsub(object)[t.loc, g.state], deci), "\n")
	}
	cat("State exit subdistributions not estimated for state(s)", no.exit)
	cat("\n\n")


        if (var==TRUE) {
            cat("Variance estimates for normalized state exit distributions: \n")
            for (state in exit.st) {
                g.state <- paste("G", state)
                cat("State ", state, ": ", round(Gsub.var(object)[t.loc, g.state], deci), "\n")
            }
            cat("\n")
        }

        if (var==TRUE) {
            return(invisible(list(entry.sub = Fsub(object)[t.loc, ],
                                  exit.sub = Gsub(object)[t.loc, ],
                                  entry.var.sub = Fsub.var(object)[t.loc, ],
                                  exit.var.sub = Gsub.var(object)[t.loc, ])))
        } else {
            return(invisible(list(entry.sub = Fsub(object)[t.loc, ],
                                  exit.sub = Gsub(object)[t.loc, ])))
        }

    } ## End of Subdistribution Functions

} ## end of 'EntryExit'


####################################################################
##                                                                ##
##              MAIN FUNCTION - msSurv                            ##
##                                                                ##
####################################################################

## NOTE - Censoring state ASSUMED to be state 0
msSurv <- function(Data, tree, cens.type="ind", LT=FALSE, bs=FALSE, B=200) {

    ## cens.type = "ind" or "dep" - used for D-S calculation
    cens.type <- match.arg(cens.type, c("ind", "dep"))

    if (any(!(c("id", "stop", "start.stage", "end.stage")%in%colnames(Data))))
        stop("Incorrect column names for 'Data'.  Column names should have 'id', 'stop', 'start.stage', and 'end.stage'.")

    if (!("start" %in% colnames(Data)) & LT==TRUE)
	stop("The 'start' times must be specified for left truncated data.")

    if (bs==TRUE & length(unique(Data$id)) < 10)
        stop("At least 10 subjects required for the bootstrap")

    ## 02/20/15: Check if LT present and give warning if LT=TRUE is not specified
    ## browser()
    if("start" %in% colnames(Data)) {
        init.start <- with(Data, tapply(start, id, min))
        if (max(init.start) > 0 & LT == FALSE) {
            warning("It appears left truncation is present.  If so please specify the 'LT=TRUE'\n  argument to 'msSurv'")
        }
    }

    if (!("start" %in% colnames(Data)) & LT==FALSE)  Data <- Add.start(Data)

    ## Check to ensure that start times >= stop times
    if (any(Data$start >= Data$stop)) stop("'start' times must be < 'stop' times.")


    ## starting probabilities based on initial states for all individuals on study BEFORE
    ## 1st obs transition time
    idx <- which(Data$start < min(Data$stop[Data$end.stage!=0]))
    start.cnts  <- table(factor(Data$start.stage[idx], levels=nodes(tree), labels=nodes(tree)))
    start.probs <- prop.table(start.cnts)

    n <- length(unique(Data$id)) ## number of individuals in sample
    ns <- length(nodes(tree)) ## number of states

    ## adding dummy states for censoring and (if needed) left truncation
    Cens <- Add.States(tree, LT)
    if (LT) {
	Data <- LT.Data(Data)
	cp <- CP(tree, Cens$treeLT, Data, Cens$nt.states.LT)
    } else {
        cp <- CP(tree, Cens$tree0, Data, Cens$nt.states)
    }


    if (cens.type == "ind") {
        ## NOTE - Need condition for LT as above with CP function call??
        ds.est <- DS.ind(Cens$nt.states, cp$dNs, cp$sum_dNs, cp$Ys)
    } else {
        if (LT) {
            ds.est <- DS.dep(Data, Cens$treeLT, Cens$nt.states.LT, cp$dNs, cp$sum_dNs, cp$Ys, LT)
        } else {
            ds.est <- DS.dep(Data, Cens$tree0, Cens$nt.states, cp$dNs, cp$sum_dNs, cp$Ys, LT)
        }
    }


    ## Reduces CP to event times and removes LT and Cens states
    cp.red <- Red(tree, cp$dNs, cp$Ys, cp$sum_dNs, ds.est$dNs.K, ds.est$Ys.K, ds.est$sum_dNs.K)
    et <- as.numeric(rownames(cp.red$dNs))

    ## made a list of all possible transitions
    res.ci2 <- strsplit(colnames(cp.red$dNs), " ") ## string splits names
    a <- sapply(res.ci2, function(x) x[2])         ## pull 1st number
    b <- sapply(res.ci2, function(x) x[3])         ## pull 2nd number
    pos.trans <- paste(a, b, sep=" ")
    stay <- paste(Cens$nt.states, Cens$nt.states, sep=" ")
    pos.trans <- sort(c(stay, pos.trans)) ## sorting all possible "transitions" from P(s, t)

    ## AJ estimator based on the DS counting process which will automatically handle
    ##   dependent censoring (see ds.est above)
    AJest <- AJ.estimator(ns, tree, cp.red$dNs.K, cp.red$Ys.K, start.probs)

    ## Entry/Exit Distributions
    entry.exit <- Dist(AJest$ps, ns, tree)


    ## variances of AJs and SOPs when cens.type = "ind"
    ## if cens.type = "dep" then need 'bs=TRUE' to calculate variance
    if (cens.type=="ind") {
        variances <- var.fn(tree, ns, Cens$nt.states, cp.red$dNs, cp.red$Ys, cp.red$sum_dNs,
                            AJest$AJs, AJest$I.dA, AJest$ps)
    }
    if (bs == TRUE) {
        bs.var <- BS.var(Data, tree, ns, et, cens.type, B, LT,
                         colnames(entry.exit$Fsub), colnames(entry.exit$Gsub))
        cat("\nUsing bootstrap to calculate variances ... \n\n")
    }

    ## Assign variances below based on various conditions
    no.start.st <-  length(which(start.probs>0))
    if (cens.type=="ind" & no.start.st==1) {
        cov.AJs <- variances$cov.AJs; var.sop <- variances$var.sop; cov.dA <- variances$cov.dA
        if (bs==TRUE) {
            Fsub.var <- bs.var$Fsub.var; Gsub.var <- bs.var$Gsub.var ## subdistribution vars
            Fnorm.var <- bs.var$Fnorm.var; Gnorm.var <- bs.var$Gnorm.var ## normalized
        } else {
            Fsub.var <- Fnorm.var <- Gsub.var <- Gnorm.var <- NULL
        }
    } else if (cens.type=="ind" & no.start.st > 1) {
        ## Need Bootstrap for SOPs when no. starting states > 1
        cov.AJs <- variances$cov.AJs; cov.dA <- variances$cov.dA
        if (bs==TRUE) {
            var.sop <- bs.var$var.sop
            Fsub.var <- bs.var$Fsub.var; Gsub.var <- bs.var$Gsub.var ## subdistribution vars
            Fnorm.var <- bs.var$Fnorm.var; Gnorm.var <- bs.var$Gnorm.var ## normalized
        } else {
            var.sop <- Fsub.var <- Fnorm.var <- Gsub.var <- Gnorm.var <- NULL
        }
    } else if (cens.type=="dep") {
        ## Dependent Censoring - Requires Bootstrap for Variance of SOPs AND trans probs
        if (bs==TRUE) {
            cov.AJs <- bs.var$cov.AJs; cov.dA <- bs.var$cov.dA; var.sop <- bs.var$var.sop
            Fsub.var <- bs.var$Fsub.var; Gsub.var <- bs.var$Gsub.var ## subdistribution vars
            Fnorm.var <- bs.var$Fnorm.var; Gnorm.var <- bs.var$Gnorm.var ## normalized
        } else {
            cov.AJs <- cov.dA <- var.sop <- Fsub.var <- Fnorm.var <- Gsub.var <- Gnorm.var <- NULL
        }
    }

        ## Result returned below
        res <- new("msSurv",
                   ## State information
                   tree=tree, ns=ns, et=et, pos.trans=pos.trans, nt.states=Cens$nt.states,
                   ## Counting process information
                   dNs=cp.red$dNs, Ys=cp.red$Ys, sum_dNs=cp.red$sum_dNs,
                   dNs.K=cp.red$dNs.K, Ys.K=cp.red$Ys.K, sum_dNs.K=cp.red$sum_dNs.K,
                   ## State occupation probs and A-J estimators
                   ps=AJest$ps, AJs=AJest$AJs, I.dA=AJest$I.dA,
                   cov.AJs=cov.AJs, var.sop=var.sop, cov.dA=cov.dA,
                   ## State entry and exit distributions
                   ## NOTE - these might be NULL (eg, if have recursion)
                   Fnorm=entry.exit$Fnorm, Fsub=entry.exit$Fsub,
                   Gnorm=entry.exit$Gnorm, Gsub=entry.exit$Gsub,
                   Fnorm.var=Fnorm.var,  Fsub.var=Fsub.var,
                   Gnorm.var=Gnorm.var,  Gsub.var=Gsub.var)

    return(res)
}





