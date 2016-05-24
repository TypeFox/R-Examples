############################################################
##                 Adding Start Times                     ##
############################################################

## Adds starting time to Data
Add.start <- function(Data) {

    Data$start <- 0
    idx <- which(table(Data$id)>1)

    for (i in names(idx)) {
        ab <- Data[which(Data$id==i), ]
        ab <- with(ab, ab[order(ab$stop), ])
        ab2 <- which(Data$id==i) ## row numbers in Data
        start2 <- vector(length=length(ab2))
        start2[1] <- 0
        start2[2:length(ab2)] <- ab$stop[1:length(ab2)-1]
        Data$start[ab2] <- start2
    } ## end of for loop

    new.data <- data.frame(id=Data$id, start=Data$start, stop=Data$stop,
                           start.stage=Data$start.stage, end.stage=Data$end.stage)
    res <- new.data
}


####################################################################
##   Adding 'Dummy' States for Censoring and Left Truncation
####################################################################

Add.States <- function(tree, LT) {

    ## Adding censoring state to Nodes & Edges
    Nodes <- c(nodes(tree), "0")
    Edges <- edges(tree)
    Edges[["0"]] <- character(0)

    nt.states <- names(which(sapply(Edges, function(x) length(x)>0))) ## nonterminal states

    for (stage in nt.states) {
        Edges[[stage]] <- c("0", Edges[[stage]])
    }

    ## tree for censored data
    tree0 <- new("graphNEL", nodes=Nodes, edgeL=Edges, edgemode="directed")

    ## Adding "Left Truncated" State
    if (LT) {
        Nodes <- c(nodes(tree0), "LT")
        Edges[["LT"]] <- nt.states
        nt.states.LT <- names(which(sapply(Edges, function(x) length(x)>0))) ## nonterminal states
        treeLT <- new("graphNEL", nodes=Nodes, edgeL=Edges, edgemode="directed")
    }

    if (LT) {
        return(list(tree0=tree0, nt.states=nt.states, nt.states.LT=nt.states.LT, treeLT=treeLT))
    } else {
        return(list(tree0=tree0, nt.states=nt.states))
    }

}

############################################################
##          Adding Dummy "LT" obs to Data set             ##
############################################################

LT.Data <- function(Data) {

    Data <- Data[order(Data$id), ]  ## make sure id's line up below
    ids <- unique(Data$id)
    stop.time <- with(Data, tapply(start, id, min))
    enter.st <- by(Data, Data$id, function(x) x$start.stage[which.min(x$start)])
    ## dummy initial stage
    dummy <- data.frame(id = ids, start = -1, stop = stop.time, start.stage="LT",
                        end.stage=as.character(enter.st))

    Data <- rbind(Data, dummy)
    Data <- with(Data, Data[order(id, stop), ])

    return(Data=Data)

}



############################################################
##              Counting Process & At Risk                ##
############################################################

####################################################################
## NOTE - ASSUMES stage 0 is censoring - do we need to pass
##        argument to allow otherwise??
####################################################################

## Assuming stage 0 is censored and stage 1 is the initial state
## state 'LT' for left truncated data

CP <- function(tree, tree0, Data, nt.states) {
    ## tree0=tree for uncens, tree0=tree0 for cens, tree0=treeLT for LT
    ## NOTE - nt.states includes 'LT' for LT data

    ## unique stop times in entire data set
    ## Exclude stop times of zero (for LT data)
    times <- sort(unique(Data$stop[!Data$stop==0]))

    ## names for dNs transitions
    lng <- sapply(edges(tree0)[nodes(tree0)%in%nt.states], length)
    ds <- paste("dN", rep(nodes(tree0)[nodes(tree0)%in%nt.states], lng),
                unlist(edges(tree0)[nodes(tree0)%in%nt.states]))

    ## names for at-risk calcs
    nt.states2 <- nt.states[!nt.states=="LT"]
    ys <- paste("y", nt.states2)  ## used for Ys

    ##  matrix of # of transitions, initialize to zeros
    dNs <- matrix(0, nrow=length(times), ncol=length(ds))

    ##  matrix of total # of transitions from a state, initialize to zeros
    sum_dNs <- matrix(0, nrow=length(times), ncol=length(nt.states))

    ##  matrix of at-risk sets for each stage at each time
    Ys <- matrix(NA, nrow=length(times), ncol=length(ys))

    ## names of rows/columns for vectors/matrices
    rownames(dNs) <- rownames(sum_dNs) <- rownames(Ys) <- times
    colnames(dNs) <- ds
    colnames(Ys) <- ys
    colnames(sum_dNs) <- paste("dN", nt.states, ".")

    ## Calculations for transitions 'dNs'
    ## Outer loop = all nodes w/transitions into them
    nodes.into <- nodes(tree0)[sapply(inEdges(tree0), function(x) length(x) > 0)]
    for (i in nodes.into) {
        ## Inner loop = all nodes which transition into node i
        nodes.from <- inEdges(tree0)[[i]]
        for (j in nodes.from) {
            nam2 <- paste("dN", j, i)
            idx <- which(Data$end.stage==i & Data$start.stage==j)
            tmp.tab <- table(Data$stop[idx][!Data$stop[idx]==0])  ##  no. trans at each trans time
            dNs[names(tmp.tab), nam2] <- tmp.tab
        }
    }


    ## start counts at time == 0 for below
    start.stages <- Data$start.stage[Data$start==0]
    start.stages <- factor(start.stages, levels=nt.states2, labels=nt.states2)
    start.cnts <- table(start.stages)

    ## Calculations for at-risk sets 'Ys'
    ## only need for non-terminal nodes - use 'nt.states2' to exclude 'LT' state
    for (i in nt.states2) {

        n <- start.cnts[i]
        nam <-paste("y", i)

        if (length(inEdges(tree0)[[i]]) > 0)
            into.node <- paste("dN", inEdges(tree0)[[i]], i) else into.node <- NULL
        if (length(edges(tree0)[[i]])>0)
            from.node <- paste("dN", i, edges(tree0)[[i]]) else from.node <- NULL

        Ys[, nam] <- c(n, n + cumsum(rowSums(dNs[, into.node, drop=FALSE])) -
                      cumsum(rowSums(dNs[, from.node, drop=FALSE])))[-(nrow(Ys)+1)]

    } ## end of loop for Ys

    ## Counting transitions from different states (ie: state sums)
    sum_dNs <- matrix(nrow=nrow(dNs), ncol=length(nt.states))
    rownames(sum_dNs) <- rownames(dNs) ##
    colnames(sum_dNs) <- paste("dN", nt.states, ".")
    a <- strsplit(colnames(sum_dNs), " ")
    a2 <- strsplit(colnames(dNs), " ")
    uni <- unique(sapply(a, function(x) x[2]))##  gives the unique states exiting

    ## browser()
    ## bug fix 02/17/2015
    ## can't counting transitions into 'censored' as part of sum_dNs ...
    ## NOTE - Need to change to '0' to 'Cens' ...
    for (i in uni) { ## calculating the dNi.s
        b <- which(sapply(a, function(x) x[2]==i))
        b2 <- which(sapply(a2, function(x) x[2]==i & x[3]!=0))
        sum_dNs[, b] <- rowSums(dNs[, b2, drop=FALSE])
    } ## end of for loop for calculating dNi.s

    list(dNs=dNs, Ys=Ys, sum_dNs=sum_dNs)

} ## end of function


############################################################
##            Datta-Satten Estimation                     ##
############################################################

## for INDEPENDENT censoring
## NOTE - Dropped 'Cens' argument (assumed to be 0) since don't allow
##        this flexibility elsewhere
DS.ind <- function(nt.states, dNs, sum_dNs, Ys) {
    ## Calculating dNs, sum_dNs, and Y from D-S(2001) paper
    ## Dividing dNs*, sum_dNs*, & Y* by K to get dNs^, sum_dNs^, & Ys^
    ## Make sure nt.states is from the non-LT

    res <- strsplit(colnames(dNs), " ") ## string splits names
    res2 <- strsplit(colnames(Ys), " ")  ## string split names of Ys
    res3 <- strsplit(colnames(sum_dNs), " ") ## string splits names of dNs

    ##  Column indicator for transitions into censored states, needed for D-S est
    DS.col.idx <- which(sapply(res, function(x) x[3]==0))
    ##  Indicator for total at-risk in each transition state
    DS2.col.idx <- which(sapply(res2, function(x) x[2]%in%nt.states))
    ##  Indicator for total transitions out of each transition state
    DS3.col.idx <- which(sapply(res3, function(x) x[2]%in%nt.states))

    K <- vector(length=nrow(dNs))
    dN0 <- rowSums(dNs[, DS.col.idx, drop=FALSE])
    Y0 <- rowSums(Ys[, DS2.col.idx, drop=FALSE]) ## those at risk of being censored
    N.Y <- ifelse(dN0/Y0=="NaN", 0, dN0/Y0)
    colnames(N.Y) <- NULL
    H.t <- cumsum(N.Y) ## calculating the hazard
    k <- exp(-H.t)
    K <- c(1, k[-length(k)])

    dNs.K <- dNs/K ## D-S dNs
    Ys.K <- Ys/K ## D-S Ys
    sum_dNs.K <- sum_dNs/K

    res <- list(dNs.K=dNs.K, Ys.K=Ys.K, sum_dNs.K=sum_dNs.K)
    return(res)

}


## for DEPENDENT censoring
## NOTE - Dropped 'Cens' argument (assumed to be 0) since don't allow
##        this flexibility elsewhere
DS.dep <- function(Data, tree0, nt.states, dNs, sum_dNs, Ys, LT) {
    ## tree0=tree for uncens, tree0=tree0 for cens, tree0=treeLT for LT

    ####################################################################
    ## STEPS in Dependent Censoring Case:
    ## 1. Calculate state-dependent survival functions for censoring
    ## 2. Calculate K_i(T_ik-) - Weighting (K) value for each individual
    ##     i just prior to each of their observed transition times
    ## 3. Calculate K_i(t-) - Weighting (K) value for each individual
    ##     i at EVERY time prior to their observed censoring time (or last
    ##     transition time if to terminal node)
    ####################################################################


    ####################################################################
    ## 1. Calculate state-dependent survival functions for censoring
    ####################################################################
    res <- strsplit(colnames(dNs), " ") ## string splits names
    res2 <- strsplit(colnames(Ys), " ")  ## string split names of Ys
    res3 <- strsplit(colnames(sum_dNs), " ") ## string splits names of dNs

    ##  Column indicator for transitions into censored states, needed for D-S est
    DS.col.idx <- which(sapply(res, function(x) x[3]==0))
    ##  Indicator for total at-risk in each transition state
    DS2.col.idx <- which(sapply(res2, function(x) x[2]%in%nt.states))
    ##  Indicator for total transitions out of each transition state
    DS3.col.idx <- which(sapply(res3, function(x) x[2]%in%nt.states))

    dN0 <- dNs[, DS.col.idx, drop=FALSE]  ## Think need the drop=FALSE option?
    Y0 <- Ys[, DS2.col.idx, drop=FALSE] ## those at risk of being censored

    N.Y <- ifelse(dN0/Y0=="NaN", 0, dN0/Y0)
    colnames(N.Y) <- paste(colnames(dN0), "/", colnames(Y0))

    H.t <- apply(N.Y, 2, function(x) cumsum(x))
    Sj.cens <- exp(-H.t)   ## think need rbind(rep(1, XX), exp(-H.t))
    Sj.cens <- rbind(rep(1, ncol(Sj.cens)), Sj.cens) ## [,-nrow(Sj.cens)])
    rownames(Sj.cens)[1] <- 0
    colnames(Sj.cens) <- sapply(res2, function(x) x[[2]])
    times <- as.numeric(rownames(Sj.cens))


    ####################################################################
    ## 2. Calculate K_i(T_ik-) - Weighting (K) value for each individual
    ##     i just prior to each of their observed transition times
    ## 3. Calculate K_i(t-) - Weighting (K) value for each individual
    ##     i at EVERY time prior to their observed censoring time (or last
    ##     transition time if to terminal node)
    ####################################################################

    Data$dN.K <- rep(0, nrow(Data))
    ## Default to 0 = censoring value

    ## Need a matrix of states occupied by each subject at each time ...
    state.mat <- matrix(0, nrow = length(unique(Data$id)), ncol = nrow(Sj.cens))
    rownames(state.mat) <- unique(Data$id)
    colnames(state.mat) <- rownames(Sj.cens)
    ## Need a matrix of K values for each subject at each time ...
    Y.K.mat <- matrix(0, nrow = length(unique(Data$id)), ncol = nrow(Sj.cens))
    rownames(Y.K.mat) <- unique(Data$id)
    colnames(Y.K.mat) <- rownames(Sj.cens)


    for (i in 1:nrow(state.mat)) {  ## loop through subjects ...
        idx.subj <- which(Data$id %in% rownames(state.mat)[i])
        ## if (rownames(state.mat)[i] == 199) browser()

        ## Find the indexes corresponding to observed times
        if (!LT) {
            x <- Data[idx.subj, ]
            idx.time <- which(rownames(Sj.cens) %in% x$stop)
        } else {
            x <- Data[idx.subj[-1], ]
            idx.time <- which(rownames(Sj.cens) %in% x$stop)
            idx.start <- which(rownames(Sj.cens) == Data$stop[idx.subj][1])
        }
        ## NOTE - Presumes sorted by time so idx will be sorted also (should be the case)
        for (j in 1:length(idx.time)) {

            if (j == 1) {
                ## Take idx.time[j]-1 because we want time just PRIOR ...
                ## Here calculate for the dN.K's
                col.idx <- which(colnames(Sj.cens) == x$start.stage[j])
                x$dN.K[j] <- Sj.cens[(idx.time[j]-1), col.idx]
                ## Here for calculating the Y.K's
                ## NOTE - below method for storing 'state.mat' as character may not
                ##        be most efficient and maybe can change to POSITION of 'nt.states'
                ##        Need to be careful of order though
                ## NOTE - Here modify to account for LT (replace '1' with something else for LT)
                if (!LT) {
                    state.mat[i, 1:(idx.time[j]-1)] <- as.character(x$start.stage[j])
                    Y.K.mat[i, 1:(idx.time[j]-1)] <- Sj.cens[1:(idx.time[j]-1), col.idx]
                } else {
                    state.mat[i, idx.start:(idx.time[j]-1)] <- as.character(x$start.stage[j])
                    Y.K.mat[i, idx.start:(idx.time[j]-1)] <- Sj.cens[idx.start:(idx.time[j]-1), col.idx]
                }

            } else {
                ## Here calculate for the dN.K's
                ## Note the inclusion of the 'correction' factor:
                ##    S_{stage,C}(idx.time[j]) / S_{stage,C}(idx.time[j-1])
                ## This accounts for past transition history of subject i
                col.idx <- which(colnames(Sj.cens) == x$start.stage[j])
                x$dN.K[j] <- x$dN.K[(j-1)] * (Sj.cens[(idx.time[j]-1), col.idx]) /
                    (Sj.cens[(idx.time[j-1]-1), col.idx])
                ## Here for calculating the Y.K's
                state.mat[i, idx.time[j-1]:(idx.time[j]-1)] <- as.character(x$start.stage[j])
                Y.K.mat[i, idx.time[j-1]:(idx.time[j]-1)] <- Sj.cens[idx.time[j-1]:(idx.time[j]-1), col.idx] *
                                                               x$dN.K[(j-1)] /
                                                                  (Sj.cens[(idx.time[j-1]-1), col.idx])
            }
        }  ## End of 'j' loop (time)
        if (!LT) {
            Data$dN.K[idx.subj] <- x$dN.K
        } else {
            Data$dN.K[idx.subj[-1]] <- x$dN.K
        }
    }  ## End of 'i' loop (subject)


    ## Calculations for Ys.K
    ## Need to sum up 1/K's
    Ys.K <- Ys
    for (i in nt.states) {
        K.idx <- which(sapply(strsplit(colnames(N.Y), " "), function(x) x[2]==i))
        Ys.idx <- which(sapply(res2, function(x) x[2]==i))
        ## dNs.K[, dN.idx] <- dNs[, dN.idx]/K[, K.idx]
        ## sum_dNs.K[, sum_dNs.idx] <- sum_dNs[, sum_dNs.idx]/K[, K.idx]
        ## THIS SHOULD WORK FOR Ys.K ...
        Ys.K[, Ys.idx] <- colSums((state.mat[,-ncol(state.mat)]==i) * 1/Y.K.mat[,-ncol(Y.K.mat)], na.rm=TRUE)
        ## Ys[, Ys.idx]/K[, K.idx]
    }

    ## Calculations for transitions 'dNs'
    dNs.K <- dNs
    nodes.into <- nodes(tree0)[sapply(inEdges(tree0), function(x) length(x) > 0)]
    ## Outer loop = all nodes w/transitions into them
    for (i in nodes.into) {
        ## Inner loop = all nodes which transition into node i
        nodes.from <- inEdges(tree0)[[i]]
        for (j in nodes.from) {
            nam2 <- paste("dN", j, i)
            idx <- which(Data$end.stage==i & Data$start.stage==j)
            tmp.dNK <- tapply(Data$dN.K[idx][!Data$stop[idx]==0],
                              Data$stop[idx][!Data$stop[idx]==0],
                              function(x) sum(1/x))
            dNs.K[names(tmp.dNK), nam2] <- tmp.dNK
        }
    }

    ## Now calculate for SUM of dNs.K (sum_dNs.K)
    ## Counting transitions from different states (ie: state sums)
    sum_dNs.K <- matrix(nrow=nrow(dNs.K), ncol=length(nt.states))
    rownames(sum_dNs.K) <- rownames(dNs.K) ##
    colnames(sum_dNs.K) <- paste("dN", nt.states, ".")
    a <- strsplit(colnames(sum_dNs.K), " ")
    a2 <- strsplit(colnames(dNs.K), " ")
    uni <- unique(sapply(a, function(x) x[2])) ## gives the unique states exiting
    for (i in uni) { ## calculating the dNi.s
        b <- which(sapply(a, function(x) x[2]==i))
        b2 <- which(sapply(a2, function(x) x[2]==i))
        sum_dNs.K[, b] <- rowSums(dNs.K[, b2, drop=FALSE])
    } ## end of for loop for calculating dNi.s

    res <- list(dNs.K=dNs.K, Ys.K=Ys.K, sum_dNs.K=sum_dNs.K)
    return(res)
}


############################################################
##           Reducing dNs & Ys to event times             ##
############################################################

Red <- function(tree, dNs, Ys, sum_dNs, dNs.K, Ys.K, sum_dNs.K) {

    ## tree is original tree currently inputted by user
    ## dNs, sum.dNS, & Ys come from CP function
    ## K comes from DS

    ## reducing dNs & Ys to just event times & noncens/nontruncated states
    res <- strsplit(colnames(dNs), " ") ## string splits names
    ##  looks at noncensored columns
    col.idx <- which(sapply(res, function(x) x[2]%in%nodes(tree) & x[3]%in%nodes(tree)))
    ## identifies times where transitions occur
    row.idx <- which(apply(dNs[, col.idx, drop=FALSE], 1, function(x) any(x>0)))
    dNs.et <- dNs[row.idx, col.idx, drop=FALSE] ## reduces dNs

    res2 <- strsplit(colnames(Ys), " ") ## string split names of Ys
    nt.states.f <- names(which(sapply(edges(tree), function(x) length(x)>0))) ## nonterminal states
    col2.idx <- which(sapply(res2, function(x) x[2]%in%nt.states.f)) ## ids nonterminal columns
    Ys.et <- Ys[row.idx, col2.idx, drop=FALSE] ## reduces Ys

    col3.idx <- which(sapply(strsplit(colnames(sum_dNs), " "), function(x) x[2]%in%nodes(tree)))
    sum_dNs.et <- sum_dNs[row.idx, col3.idx, drop=FALSE]

    dNs.K.et <- dNs.K[row.idx, col.idx, drop=FALSE]
    Ys.K.et <- Ys.K[row.idx, col2.idx, drop=FALSE]
    sum_dNs.K.et <- sum_dNs.K[row.idx, col3.idx, drop=FALSE]

    ans <- list(dNs=dNs.et, Ys=Ys.et, sum_dNs=sum_dNs.et, dNs.K=dNs.K.et, Ys.K=Ys.K.et, sum_dNs.K=sum_dNs.K.et)
    return(ans)

}


############################################################
##     AJ Estimates and State Occupation Probabilities    ##
############################################################

AJ.estimator <- function(ns, tree, dNs.et, Ys.et, start.probs) {

    ## currently ns is defined in main function
    ## tree needs to be uncensored tree
    cum.tm <- diag(ns)
    colnames(cum.tm) <- rownames(cum.tm) <- nodes(tree)

    ps <- matrix(NA, nrow=nrow(dNs.et), ncol=length(nodes(tree)))
    rownames(ps) <- rownames(dNs.et); colnames(ps) <- paste("p", nodes(tree))
    all.dA <- all.I.dA <- all.AJs <- array(dim=c(ns, ns, nrow(dNs.et)),
                                           dimnames=list(rows=nodes(tree),
                                           cols=nodes(tree), time=rownames(dNs.et)))

    for (i in 1:nrow(dNs.et)) { ##loop through times


        I.dA <- diag(ns) ## creates trans matrix for current time
        dA <- matrix(0, nrow=ns, ncol=ns)
        colnames(I.dA) <- rownames(I.dA) <- colnames(dA) <- rownames(dA) <- nodes(tree)

        idx <- which(dNs.et[i, , drop=FALSE]>0)  ## transition time i
        t.nam <- colnames(dNs.et)[idx] ## gets names of transitions (ie:  dN##)
        tmp <- strsplit(t.nam, " ")    ## splits title of dN##
        start <- sapply(tmp, function(x) x[2])
        end <- sapply(tmp, function(x) x[3])  ## pulls start & stop states as character strings
        idxs <- matrix(as.character(c(start, end)), ncol=2)
        idxs2 <- matrix(as.character(c(start, start)), ncol=2)

        dA[idxs] <- dNs.et[i, idx]/Ys.et[i, paste("y", start)]
        if (length(idx)==1) {
            dA[start, start] <- -dNs.et[i, idx]/Ys.et[i, paste("y", start)]
        } else {
            dA[idxs2] <- -rowSums(dA[start, ])
        }

        I.dA <- I.dA + dA ## I+dA (transition) matrix

        all.dA[, , i] <- dA     ## stores all dA matrices
        all.I.dA[, , i] <- I.dA ## array for storing all tran matrices

        cum.tm <- cum.tm %*% I.dA  ## Multiply current.tm and cum.tm to get matrix for current time
        all.AJs[, , i] <- cum.tm   ## A-J estimates, stored in array

        ## multiply by start.probs to allow for starting states other than '1'
        ps[i, ] <- start.probs%*%all.AJs[, , i] ## state occupation probabilities

    } ## end of loop

    list(ps=ps, AJs=all.AJs, I.dA=all.I.dA)
} ## end of function


############################################################
##           State Entry/Exit Distributions               ##
############################################################

Dist <- function(ps, ns, tree) {
    ## ps from AJ.estimator function
    ## tree needs to be uncensored tree

    ## Recursive function to detect whether system is cyclic
    downstream <- function(nodes, tree) {
        all.edges <- unique(unlist(edges(tree)[nodes]))
        all.edges <- all.edges[!(all.edges %in% all.nodes)] ## remove nodes already visited
        if (length(all.edges) == 0) {
            return(NULL)
        } else {
            all.nodes <<- c(all.nodes, all.edges)
            return(c(all.edges, downstream(all.edges, tree)))
        }
    }


    later.nodes <- vector("list", ns)
    recur <- logical(ns)
    names(later.nodes) <- names(recur) <- nodes(tree)
    for (node in nodes(tree)) {
        all.nodes <- c()
        later.nodes[[node]] <- downstream(node, tree)
        recur[node] <- node %in% later.nodes[[node]]
    }
    ## NOTE - later.nodes truncated to length of those states HAVING downstream nodes

    ## separated = TRUE if boundary for later nodes of given node is unique
    ##        OR = TRUE if node is TERMINAL node
    separated <- logical(ns)
    names(separated) <- nodes(tree)
    for (i in 1:length(later.nodes)) {
        tmp <- boundary(later.nodes[[i]], tree)
        node <- names(later.nodes)[i]
        separated[node] <- (length(unique(unlist(tmp))) == 1) ## TRUE if separated
    }
    terminal <- names(separated)[!names(separated) %in% names(later.nodes)]
    separated[terminal] <- TRUE

    ## Estimate for states with recur == FALSE and separated = TRUE
    est.states <- !recur & separated
    if (sum(est.states) == 0) {
        cat("No states eligible for entry/exit distribution calculation. \n")
        return(list(Fnorm=NULL, Gsub=NULL, Fsub=NULL, Gnorm=NULL))
    }

    ## initial states have no Fs (entry dist)
    Fs.states <- which(!sapply(inEdges(tree), function(x) !length(x)>0) & est.states)

    ## terminal states have no Gs (exit dist)
    Gs.states <- which(!sapply(edges(tree), function(x) !length(x)>0) & est.states)
    ## NOTE - Fs.states and Gs.states give POSITION of node

    if (length(Fs.states) > 0) {
        Fnorm <- Fsub <- matrix(0, nrow=nrow(ps), ncol=length(Fs.states)) ## entry distn
        rownames(Fnorm) <- rownames(Fsub) <- rownames(ps)
        colnames(Fnorm) <- colnames(Fsub) <- paste("F", nodes(tree)[Fs.states])
    } else {
        cat("\nNo states eligible for entry distribution calculation.\n")
        Fsub <- Fnorm <- NULL
    }

    if (length(Gs.states) > 0) {
        Gsub <- Gnorm <- matrix(0, nrow=nrow(ps), ncol=length(Gs.states)) ## exit distn
        rownames(Gnorm) <- rownames(Gsub) <- rownames(ps)
        colnames(Gnorm) <- colnames(Gsub) <- paste("G", nodes(tree)[Gs.states])
    } else {
        cat("\nNo states eligible for exit distribution calculation.\n")
        Gsub <- Gnorm <- NULL
    }

    ## Fs (entry dist) calculation
    if (length(Fs.states) > 0) {
         cat("\nEntry distributions calculated for states", nodes(tree)[Fs.states], ".\n")
        for (i in 1:length(Fs.states)) {
            node <- nodes(tree)[Fs.states[i]]
            later.stages <- names(acc(tree, node)[[1]])
            stages <- c(node, later.stages)

            Fsub[, i] <- f.numer <- rowSums(ps[, paste("p", stages), drop=FALSE])
            Fnorm[, i] <- f.numer/f.numer[length(f.numer)]
        }
    }

    ## Gs (exit dist) calculation
    if (length(Gs.states) > 0) {
        cat("\nExit distributions calculated for states", nodes(tree)[Gs.states], ".\n")
        for (i in 1:length(Gs.states)) {
            node <- nodes(tree)[Gs.states[i]]
            later.stages <- names(acc(tree, node)[[1]])
            stages <- c(node, later.stages)
            f.numer <- rowSums(ps[, paste("p", stages), drop=FALSE])

            Gsub[, i] <- g.numer <- rowSums(ps[, paste("p", later.stages), drop=FALSE])
            Gnorm[, i] <- g.numer/f.numer[length(f.numer)]
        }
    }

    list(Fnorm=Fnorm, Gsub=Gsub, Fsub=Fsub, Gnorm=Gnorm)
} ## end of function

############################################################
###         Variance of AJ estimates                    ####
############################################################

var.fn <- function(tree, ns, nt.states, dNs.et, Ys.et, sum_dNs, AJs, I.dA, ps) {

    ## covariance matrices
    state.names <- nodes(tree)
    trans.names <- paste(rep(state.names, ns), rep(state.names, each=ns))
    cov.dA <- cov.AJs <- array(0, dim = c(ns^2, ns^2, nrow(dNs.et)))
    dimnames(cov.dA) <- dimnames(cov.AJs) <- list(trans.names, trans.names, rownames(dNs.et))
    res.array <- array(0, dim=c(ns^2, ns^2))
    colnames(res.array) <- rownames(res.array) <- trans.names


    ## Ident matrix for Kronecker products
    bl.Id <- diag(1, (ns)^2)
    Id <- diag(1, ns)

    ## matrix for var matrix of state occup prob
    var.sop <- matrix(0, nrow=nrow(dNs.et), ncol=ns)
    colnames(var.sop) <- paste("Var", "p", state.names)
    rownames(var.sop) <- rownames(ps)

    ## NOTE - see note below for use of v.p ...
    v.p <- matrix(0, ns, ns)

    for (i in 1:nrow(dNs.et)) { ## loop through times


        ## VARIANCE OF A-J ( TRANS PROB MATRIX P(0, t) )
        ## Equation 4.4.20 in Anderson 1993

        ## loop on the blocks (g) - only needed for non-terminal states
	for (outer in nt.states) {

            ## find positioning of outer in state.names
            idx.outer <- which(state.names==outer)
            ## tmp matrix for calculating variance of transitions in each block (g)
            tm <- matrix(0, nrow=ns, ncol=ns)
            colnames(tm) <- rownames(tm) <- state.names

            ## loop in the blocks
            for (j in 1:ns) {

                ## this just fills in upper diagonal matrix
                ## use symmetry to fill in rest
                for (k in j:ns) {

                    statej <- state.names[j]
                    statek <- state.names[k]
                    outer.Ys <- paste("y", outer)
                    outer.sum_dNs <- paste("dN", outer, ".")

                    if (Ys.et[i, outer.Ys]==0) {  ## if Y_g = 0 then covariance = 0
        	  	tm[j, k] <- 0
	        	next
                    }


                    if (statej == outer & statek == outer) {  ## 3rd formula
			tm[j, k] <- (Ys.et[i, outer.Ys] - sum_dNs[i, outer.sum_dNs])*sum_dNs[i, outer.sum_dNs] / Ys.et[i, outer.Ys]^3

                    }  else if (statej == outer & statek != outer) {  ## 2nd formula
                        name <- paste("dN", outer, statek)
			if (!name%in%colnames(dNs.et)) next
			tm[j, k] <- -(Ys.et[i, outer.Ys] - sum_dNs[i, outer.sum_dNs])*dNs.et[i, name] / Ys.et[i, outer.Ys]^3
                    } else if (statej != outer & statek == outer) {  ## 2nd formula pt 2, for recurrent
                        name <- paste("dN", outer, statej)
			if (!name%in%colnames(dNs.et)) next
			tm[j, k] <- -(Ys.et[i, outer.Ys] - sum_dNs[i, outer.sum_dNs])*dNs.et[i, name]/Ys.et[i, outer.Ys]^3
                    } else { ## 1st formula
			namek <- paste("dN", outer, statek)
			namej <- paste("dN", outer, statej)
			if (!(namej%in%colnames(dNs.et) & namek%in%colnames(dNs.et))) next
			tm[j, k] <- (ifelse(j==k, 1, 0)*Ys.et[i, outer.Ys] - dNs.et[i, namej])*dNs.et[i, namek]/Ys.et[i, outer.Ys]^3
                    } ## end of if/else statements
                } ## end of k loop
            } ## end of j loop

            tm[lower.tri(tm)] <- t(tm)[lower.tri(tm)]

            res.array[(seq(1, ns*(ns-1)+1, by=ns) + idx.outer - 1), (seq(1, ns*(ns-1)+1, by=ns) + idx.outer - 1)] <- tm

	}## end of outer loop

        ## array holding var-cov matrix for I+dA matrix at each time (differential of NA estim)
	cov.dA[, , i] <- res.array

	if (i==1) {
            cov.AJs[, , i] <- bl.Id%*% cov.dA[, , i] %*% bl.Id
        } else {
            cov.AJs[, , i] <- (t(I.dA[, , i]) %x% Id) %*% cov.AJs[, , i-1] %*%((I.dA[, , i]) %x% Id) +
            (Id %x% AJs[, , i-1]) %*% cov.dA[, , i]  %*% (Id%x% t(AJs[, , i-1]))
        }
        ## note:  AJs[, , i-1] corresponds to P(0, t-)
        ## cov.AJs is the var/cov est of P(0, t)

        ## calculating the variance of state occupation prob
	for (j in state.names) { ## loop through states

            name <- paste(state.names[1], j)
            part1 <- var.pkj0t <- cov.AJs[name, name, i]
            res <- strsplit(colnames(ps), " ")
            col.idx <- which(sapply(res, function(x) x[2]== j)) ## looks at state transitioned to
            b.t <- AJs[, col.idx, i] ## creating vector of col for current state from trans prob

            ## NOTE - This is ZERO when all indiv start from single state ...
            part2 <- t(b.t)%*%v.p%*%b.t ## should be 0 when P(0, t)
            ## right now forced to be 0 by the way v.p defined outside of time loop
            res.varp <- part1 + part2      ## calculating var/cov matrix for time i, state j
            var.sop[i, col.idx] <- res.varp ## storing var/cov calc for time i,n state j

	} ## closes states loop
    } ## end of time loop

    list(cov.AJs=cov.AJs, cov.dA=cov.dA, var.sop=var.sop)

}## end of function


#########################################################
##                BS Variance                          ##
#########################################################

## Needed For:
## 1. Dependent Censoring
## 2. Entry / Exit functions
## 3. SOPs when > 1 starting state

BS.var <- function(Data, tree, ns, et, cens.type, B, LT, entry.states, exit.states) {

    if(!is.null(entry.states))
        entry.states <- sapply(strsplit(entry.states, " "), function(x) x[2])

    if(!is.null(exit.states))
        exit.states <- sapply(strsplit(exit.states, " "), function(x) x[2])


    n <- length(unique(Data$id)) ## sample size
    ids <- unique(Data$id)

    ## storage for bootstrap estimates of transition probability matrices
    bs.est <- array(dim=c(length(nodes(tree)), length(nodes(tree)), length(et), B),
                    dimnames=list(rows=nodes(tree), cols=nodes(tree), time=et))
    bs.ps <- array(dim=c(length(et), ns, B))
    rownames(bs.ps) <- et
    colnames(bs.ps) <- paste("p", nodes(tree))

    ## For entry / exit distributions
    if (!is.null(entry.states)) {
        bs.Fnorm <- bs.Fsub <- array(dim=c(length(et), length(entry.states), B))
        colnames(bs.Fnorm) <- colnames(bs.Fsub) <- paste("F", entry.states)
    } else {
        bs.Fnorm <- bs.Fsub <- NULL
    }
    if (!is.null(exit.states)) {
        bs.Gnorm <- bs.Gsub <- array(dim=c(length(et), length(exit.states), B))
        colnames(bs.Gnorm) <- colnames(bs.Gsub) <- paste("G", exit.states)
    } else {
        bs.Gnorm <- bs.Gsub <- NULL
    }

    ## initial <- which(sapply(inEdges(tree), function(x) !length(x) > 0)) ## initial states, no Fs (entry)
    ## terminal <- which(sapply(edges(tree), function(x) !length(x) > 0)) ## terminal states, no Gs (exit)

    ## matrix for var matrix of state occup prob
    bs.var.sop <- matrix(0, nrow=length(et), ncol=ns)
    colnames(bs.var.sop) <- paste("Var", "p", nodes(tree))
    rownames(bs.var.sop) <- et

    res.array <- array(0, dim=c(ns^2, ns^2, length(et)),
                       dimnames=list(rows=paste(rep(nodes(tree), ns), sort(rep(nodes(tree), ns))),
                       cols=paste(rep(nodes(tree), ns), sort(rep(nodes(tree), ns))), time=et))

    for (b in 1:B) { ## randomly selects the indices

        ## Find the bootstrap sample, pull bs sample info from original data & put into data set
        bs <- sample(ids, n, replace=TRUE)
        bs <- factor(bs, levels=ids)
        bs.tab <- data.frame(table(bs)) ### table with the frequencies
        Data.bs <- merge(Data, bs.tab, by.x="id", by.y="bs") ## merging original data with bs freq table
        bs.id <- as.vector(unlist(apply(Data.bs[Data.bs$Freq>0, , drop=FALSE], 1, function(x) paste(x["id"], 1:x["Freq"], sep=".")))) ## creating bs id
        idx <- rep(1:nrow(Data.bs), Data.bs$Freq) ## indexing the bs sample
        Data.bs <- Data.bs[idx, ] ## creating a bs dataset
        Data.bs.originalID <- Data.bs$id
        Data.bs$id <- bs.id ## changing id column to bs.id to use functions
        Data.bs <- Data.bs[order(Data.bs$stop), ] ## ordered bs dataset

        ## Calling functions using bs dataset
        Cens <- Add.States(tree, LT)

        ## Here calculate start probs for BS data ...
        ## Data may be LT so need to account for that
        ## Put check here for whether have LT or not ...
        if (LT) {
            min.tran <- min(Data.bs$stop[Data.bs$end.stage!=0 & Data.bs$start.stage!="LT"])
            idx1 <- which(Data.bs$start.stage=="LT" & Data.bs$stop < min.tran)
            idx2 <- which(Data.bs$start.stage!="LT" & Data.bs$start < min.tran)
            start.stages <- by(Data[c(idx1, idx2), ], Data.bs$id[c(idx1, idx2)], function(x)
                               ifelse(min(x$start)<0, x$end.stage[which.min(x$start)],
                                      x$start.stage[which.min(x$start)]))
            start.stages <- factor(start.stages, levels=nodes(tree), labels=nodes(tree))
            start.cnts  <- table(start.stages)
            start.probs <- prop.table(start.cnts)
        } else {
                idx <- which(Data.bs$start < min(Data.bs$stop[Data.bs$end.stage!=0]))
                start.cnts  <- table(factor(Data.bs$start.stage[idx], levels=nodes(tree), labels=nodes(tree)))
                start.probs <- prop.table(start.cnts)
            }

        if (LT) {
            cp <- CP(tree, Cens$treeLT, Data.bs, Cens$nt.states.LT)
        }

        if (!LT) {
            cp <- CP(tree, Cens$tree0, Data.bs, Cens$nt.states)
        }

        if (cens.type == "ind") {
            ds.est <- DS.ind(Cens$nt.states, cp$dNs, cp$sum_dNs,
                             cp$Ys)
        } else {
            if (LT) {
                ds.est <- DS.dep(Data.bs, Cens$treeLT, Cens$nt.states.LT, cp$dNs,
                                 cp$sum_dNs, cp$Ys, LT)
            } else {
                ds.est <- DS.dep(Data.bs, Cens$tree0, Cens$nt.states, cp$dNs,
                                 cp$sum_dNs, cp$Ys, LT)
            }
        }
        cp.red <- Red(tree, cp$dNs, cp$Ys, cp$sum_dNs, ds.est$dNs.K,
                      ds.est$Ys.K, ds.est$sum_dNs.K)
        AJest <- AJ.estimator(ns, tree, cp.red$dNs.K, cp.red$Ys.K, start.probs)

        idx <- which(dimnames(bs.est)[[3]] %in% dimnames(AJest$I.dA)[[3]])
        idx2 <- which(!(dimnames(bs.est)[[3]] %in% dimnames(AJest$I.dA)[[3]]))
        bs.IA <- bs.est
        bs.IA[, , idx, b] <- AJest$I.dA
        bs.IA[, , idx2, b] <- diag(ns)

        bs.est[, , 1, b] <- bs.IA[, , 1, b]
        bs.ps[1, , b] <- start.probs%*%bs.est[, , 1, b]

        for (j in 2:length(et)) {
            bs.est[, , j, b] <- bs.est[, , j-1, b] %*% bs.IA[, , j, b]
            bs.ps[j, , b] <- start.probs%*%bs.est[, , j, b]
        } ## end of j for loop
        ## looks ok

        ## Entry / Exit variance as well
        ## Fs (entry dist) calculation
        if (!is.null(entry.states)) {
            for (i in 1:length(entry.states)) {
                node <- entry.states[i]
                later.stages <- names(acc(tree, node)[[1]])
                stages <- c(node, later.stages)
                bs.f.numer <- rowSums(bs.ps[, paste("p", stages), b, drop=FALSE])

                if (sum(bs.f.numer)==0)  bs.Fnorm[, i, b] <- bs.Fsub[, i, b] <- 0  else {
                    bs.Fsub[, i, b] <- bs.f.numer
                    bs.Fnorm[, i, b] <- bs.f.numer/bs.f.numer[length(bs.f.numer)]}
            }
        }
        ## Gs (exit dist) calculation
        if (!is.null(exit.states)) {
            for (i in 1:length(exit.states)) {
                node <- exit.states[i]
                later.stages <- names(acc(tree, node)[[1]])
                stages <- c(node, later.stages)
                bs.f.numer <- rowSums(bs.ps[, paste("p", stages), b, drop=FALSE])
                bs.g.numer <- rowSums(bs.ps[, paste("p", later.stages), b, drop=FALSE])

                if (sum(bs.g.numer)==0)  bs.Gnorm[, i, b] <- bs.Gsub[, i, b] <- 0 else {
                    bs.Gsub[, i, b] <- bs.g.numer
                    ## 02/25/15 - note bug here used bs.g.numer in denominator previously
                    bs.Gnorm[, i, b] <- bs.g.numer/bs.f.numer[length(bs.f.numer)]}
            } ## end of for loop
        }

    } ## end of bs loop

    ## Normalized entry / exit
    ## NOTE - dimension kept when using apply
    if (!is.null(bs.Fnorm)) {
        Fnorm.var <- apply(bs.Fnorm, c(1, 2), var)
    } else {
        Fnorm.var <- NULL
    }
    if (!is.null(bs.Gnorm)) {
        Gnorm.var <- apply(bs.Gnorm, c(1, 2), var)
    } else {
        Gnorm.var <- NULL
    }

    ## Subdistribution entry / exit
    if (!is.null(bs.Fsub)) {
        Fsub.var <- apply(bs.Fsub, c(1, 2), var)
    } else {
        Fsub.var <- NULL
    }
    if (!is.null(bs.Gsub)) {
        Gsub.var <- apply(bs.Gsub, c(1, 2), var)
    } else {
        Gsub.var <- NULL
    }

    ## SOP's
    bs.var.sop <- apply(bs.ps, c(1, 2), var)
    colnames(bs.var.sop) <- paste("Var", "p", nodes(tree))
    rownames(bs.var.sop) <- et

    state.names <- nodes(tree)
    trans.names <- paste(rep(state.names, ns), rep(state.names, each=ns))
    bs.cov <- array(dim=c(ns^2, ns^2, length(et)),
                    dimnames=list(rows=trans.names, cols=trans.names, time=et))
    bs.IdA.cov <- array(dim=c(ns^2, ns^2, length(et)),
                        dimnames=list(rows=trans.names, cols=trans.names, time=et))

    for (i in 1:length(et)) {
        bs.est.t <- matrix(bs.est[, , i, ], nrow=B, ncol=ns^2, byrow=TRUE)
        bs.IdA.t <- matrix(bs.IA[, , i, ], nrow=B, ncol=ns^2, byrow=TRUE)
        bs.cov[, , i] <- cov(bs.est.t)
        bs.IdA.cov[, , i] <- cov(bs.IdA.t)
    }  ## this for loop creates a B x (# of states)^2 x (# of event times)

    list(cov.AJs=bs.cov, var.sop=bs.var.sop, cov.dA=bs.IdA.cov,
         Fnorm.var=Fnorm.var, Gnorm.var=Gnorm.var, Gsub.var=Gsub.var, Fsub.var=Fsub.var)

} ## end of function




####################################################################
##          CONFIDENCE INTERVALS for p(t) & P(s, t)               ##
####################################################################

## x = msSurv object
MSM.CIs <- function(x, ci.level=0.95, ci.trans="linear", trans=TRUE, sop=TRUE) {

    if (ci.level < 0 || ci.level > 1)
        stop("confidence level must be between 0 and 1")

    z.alpha <- qnorm(ci.level + (1 - ci.level) / 2)

    ci.trans <- match.arg(ci.trans, c("linear", "log", "cloglog", "log-log"))

    ## CIs on state occup probability
    if (sop==TRUE) {
        CI.p <- array(0, dim=c(nrow(dNs(x)), 3, length(nodes(tree(x)))),
                      dimnames=list(rows=et(x),
                      cols=c("est", "lower limit", "upper limit"),
                      state=paste("p", nodes(tree(x)))))

        CI.p[ , 1, ] <- ps(x)
        switch(ci.trans[1],
               "linear" = {
                   CI.p[ , 2, ] <- ps(x) - z.alpha * sqrt(var.sop(x))
                   CI.p[ , 3, ] <- ps(x) + z.alpha * sqrt(var.sop(x))},
               "log" = {
                   CI.p[ , 2, ] <- exp(log(ps(x)) - z.alpha * sqrt(var.sop(x)) / ps(x))
                   CI.p[ , 3, ] <- exp(log(ps(x)) + z.alpha * sqrt(var.sop(x)) / ps(x))},
               "cloglog" = {
                   CI.p[ , 2, ] <- 1 - (1 - ps(x))^(exp(z.alpha * (sqrt(var.sop(x)) /
                                                                 ((1 - ps(x)) * log(1 - ps(x))))))
                   CI.p[ , 3, ] <- 1 - (1 - ps(x))^(exp(-z.alpha * (sqrt(var.sop(x)) /
                                                                  ((1 - ps(x)) * log(1 - ps(x))))))},
               "log-log" = {
                   CI.p[ , 2, ] <- ps(x)^(exp(-z.alpha * (sqrt(var.sop(x)) / (ps(x) * log(ps(x))))))
                   CI.p[ , 3, ] <- ps(x)^(exp(z.alpha * (sqrt(var.sop(x)) / (ps(x) * log(ps(x))))))})

        ## Need to loop through columns to apply pmax / pmin
        for (j in 1:length(nodes(tree(x)))) {
            CI.p[ , 2, j] <- pmax(CI.p[ , 2, j], 0)
            CI.p[ , 3, j] <- pmin(CI.p[ , 3, j], 1)
        } ## end states loop
    }

    ## CIs on transition probability matrices ##
    if (trans==TRUE) {
        CI.trans <- array(0, dim=c(nrow(dNs(x)), 4, length(pos.trans(x))),
                          dimnames = list(rows=et(x),
                          cols=c("est", "lower limit", "upper limit", "var.tp"),
                          trans=pos.trans(x)))

        for (j in 1:length(pos.trans(x))) { ## loop through possible transitions

            idx <- unlist(strsplit(pos.trans(x)[j], " "))
            CI.trans[ , 1, j] <- PE <- AJs(x)[idx[1], idx[2] , ]
            CI.trans[ , 4, j] <- var <- cov.AJs(x)[pos.trans(x)[j], pos.trans(x)[j], ]


            switch(ci.trans[1],
                   "linear" = {
                       CI.trans[ , 2, j] <- PE - z.alpha * sqrt(var)
                       CI.trans[ , 3, j] <- PE + z.alpha * sqrt(var)},
                   "log" = {
                       CI.trans[ , 2, j] <- exp(log(PE) - z.alpha * sqrt(var) / PE)
                       CI.trans[ , 3, j] <- exp(log(PE) + z.alpha * sqrt(var) / PE)},
                   "cloglog" = {
                       CI.trans[ , 2, j] <- 1 - (1 - PE)^(exp(z.alpha * (sqrt(var) /
                                                                       ((1 - PE) * log(1 - PE)))))
                       CI.trans[ , 3, j] <- 1 - (1 - PE)^(exp(-z.alpha * (sqrt(var) /
                                                                        ((1 - PE) * log(1 - PE)))))},
                   "log-log" = {
                       CI.trans[ , 2, j] <- PE^(exp(-z.alpha * (sqrt(var) / (PE * log(PE)))))
                       CI.trans[ , 3, j] <- PE^(exp(z.alpha * (sqrt(var) / (PE * log(PE)))))})

            CI.trans[ , 2, j] <- pmax(CI.trans[ , 2, j], 0)
            CI.trans[ , 3, j] <- pmin(CI.trans[ , 3, j], 1)

        } ## end j loop
    }

    if (trans==TRUE & sop==TRUE) {
        return(list(CI.p=CI.p, CI.trans=CI.trans))
    } else if (trans==TRUE) {
        return(list(CI.trans=CI.trans))
    } else if (sop==TRUE) {
        return(list(CI.p=CI.p))
    } else {
        cat("\nNothing to return \n\n")
    }

} ## end of function


####################################################################
##              CIs for Entry / Exit Functions                    ##
####################################################################


## NOTE - Currently Dist.CIs only called from plot method for msSurv object
##        x = msSurv object
##        type = plot.type = "entry.sub", "entry.norm", "exit.sub", or "exit.norm"
##        states = requested states for CI calculation

Dist.CIs <- function(x, ci.level=0.95, ci.trans="linear", type, states) {

    z.alpha <- qnorm(ci.level + (1 - ci.level) / 2)

    ci.trans <- match.arg(ci.trans, c("linear", "log", "cloglog", "log-log"))
    type <- match.arg(type, c("entry.sub", "entry.norm", "exit.sub", "exit.norm"))

    CI.Dist <- array(0, dim=c(length(et(x)), 3, length(states)),
                     dimnames=list(rows=et(x), cols=c("est", "lower limit", "upper limit"),
                         states=states))
    switch(type[1],
           "entry.sub" = {
               F.states <- paste("F", states)
               CI.Dist[,1,] <- Fsub(x)[,F.states]
               Dist.var <- Fsub.var(x)[,F.states]},
           "entry.norm" = {
               F.states <- paste("F", states)
               CI.Dist[,1,] <- Fnorm(x)[,F.states]
               Dist.var <- Fnorm.var(x)[,F.states]},
           "exit.sub" = {
               G.states <- paste("G", states)
               CI.Dist[,1,] <- Gsub(x)[,G.states]
               Dist.var <- Gsub.var(x)[,G.states]},
           "exit.norm" = {
               G.states <- paste("G", states)
               CI.Dist[,1,] <- Gnorm(x)[,G.states]
               Dist.var <- Gnorm.var(x)[,G.states]
           })


    switch(ci.trans[1],
           "linear" = {
               CI.Dist[ , 2, ] <- CI.Dist[,1,] - z.alpha * sqrt(Dist.var)
               CI.Dist[ , 3, ] <- CI.Dist[,1,] + z.alpha * sqrt(Dist.var)},
           "log" = {
               CI.Dist[ , 2, ] <- exp(log(CI.Dist[,1,]) - z.alpha * sqrt(Dist.var) / CI.Dist[,1,])
               CI.Dist[ , 3, ] <- exp(log(CI.Dist[,1,]) + z.alpha * sqrt(Dist.var) / CI.Dist[,1,])},
           "cloglog" = {
               CI.Dist[ , 2, ] <- 1 - (1 - CI.Dist[,1,])^(exp(z.alpha * (sqrt(Dist.var) / ((1 - CI.Dist[,1,]) * log(1 - CI.Dist[,1,])))))
               CI.Dist[ , 3, ] <- 1 - (1 - CI.Dist[,1,])^(exp(-z.alpha * (sqrt(Dist.var) / ((1 - CI.Dist[,1,]) * log(1 - CI.Dist[,1,])))))},
           "log-log" = {
               CI.Dist[ , 2, ] <- CI.Dist[,1,]^(exp(-z.alpha * (sqrt(Dist.var) / (CI.Dist[,1,] * log(CI.Dist[,1,])))))
               CI.Dist[ , 3, ] <- CI.Dist[,1,]^(exp(z.alpha * (sqrt(Dist.var) / (CI.Dist[,1,] * log(CI.Dist[,1,])))))
           })

    for (j in 1:length(states)) {
        CI.Dist[ , 2, j] <- pmax(CI.Dist[ , 2, j], 0)
        CI.Dist[ , 3, j] <- pmin(CI.Dist[ , 3, j], 1)
    }

    return(CI.Dist)

} ## end of function
