# for the current variable of interest, xname,
# create the list of variables to condition on:

create_cond_list <- function(cond, threshold, xname, input) {

   stopifnot(is.logical(cond))
   if (!cond) return(NULL)
   if (threshold > 0 & threshold < 1) {
           ctrl <- ctree_control(teststat = "quad", testtype = "Univariate", stump = TRUE)
           xnames <- names(input)
           xnames <- xnames[xnames != xname]
           ct <- ctree(as.formula(paste(xname, "~", paste(xnames, collapse = "+"), collapse = "")),
                       data = input, controls = ctrl)
           crit <- ct@tree$criterion[[2]]
           crit[which(is.na(crit))] <- 0
           return(xnames[crit > threshold])
       }
   stop()
}



## mincriterion = 0 so that complete tree is evaluated; 
## regulate size of considered tree here via, e.g., mincriterion = 0.95
## or when building the forest in the first place via cforest_control(mincriterion = 0.95)

varimp <- function (object, mincriterion = 0, conditional = FALSE, 
                    threshold = 0.2, nperm = 1, OOB = TRUE, pre1.0_0 = conditional)
{

    response <- object@responses
    if (length(response@variables) == 1 && 
        inherits(response@variables[[1]], "Surv"))
        return(varimpsurv(object, mincriterion, conditional, threshold, nperm, OOB, pre1.0_0))
    input <- object@data@get("input")
    xnames <- colnames(input)
    inp <- initVariableFrame(input, trafo = NULL)
    y <- object@responses@variables[[1]]
    if(length(response@variables) != 1)
        stop("cannot compute variable importance measure for multivariate response")

    if (conditional || pre1.0_0) {
        if(!all(complete.cases(inp@variables)))
            stop("cannot compute variable importance measure with missing values")
    }
    CLASS <- all(response@is_nominal)
    ORDERED <- all(response@is_ordinal)
    if (CLASS) {
        error <- function(x, oob) mean((levels(y)[sapply(x, which.max)] != 
            y)[oob])
    }
    else {
        if (ORDERED) {
            error <- function(x, oob) mean((sapply(x, which.max) != 
                y)[oob])
        }
        else {
            error <- function(x, oob) mean((unlist(x) - y)[oob]^2)
        }
    }

    w <- object@initweights
    if (max(abs(w - 1)) > sqrt(.Machine$double.eps))
        warning(sQuote("varimp"), " with non-unity weights might give misleading results")

    ## list for several permutations
    perror <- matrix(0, nrow = nperm*length(object@ensemble), ncol = length(xnames))
    ## this matrix is initialized with values 0 so that a tree that does not 
    ## contain the current variable adds importance 0 to its average importance
    colnames(perror) <- xnames
        for (b in 1:length(object@ensemble)){
            tree <- object@ensemble[[b]]


            ## if OOB == TRUE use only oob observations, otherwise use all observations in learning sample
            if(OOB){oob <- object@weights[[b]] == 0} else{ oob <- rep(TRUE, length(y))}
            p <- .Call("R_predict", tree, inp, mincriterion, -1L, PACKAGE = "party")
            eoob <- error(p, oob)

            ## for all variables (j = 1 ... number of variables) 
            for(j in unique(varIDs(tree))){
              for (per in 1:nperm){

                if (conditional || pre1.0_0) {
                    tmp <- inp
                    ccl <- create_cond_list(conditional, threshold, xnames[j], input)
                    if (is.null(ccl)) {
                        perm <- sample(which(oob))
                    } else {
                        perm <- conditional_perm(ccl, xnames, input, tree, oob)
                    }
                    tmp@variables[[j]][which(oob)] <- tmp@variables[[j]][perm]
                    p <- .Call("R_predict", tree, tmp, mincriterion, -1L,
                       PACKAGE = "party")
                } else {
                    p <- .Call("R_predict", tree, inp, mincriterion, as.integer(j),
                               PACKAGE = "party")
                }
                ## run through all rows of perror
                perror[(per+(b-1)*nperm), j] <- (error(p, oob) - eoob)

              } ## end of for (per in 1:nperm)
            } ## end of for(j in unique(varIDs(tree)))
        } ## end of for (b in 1:length(object@ensemble))

    perror <- as.data.frame(perror)
    #return(MeanDecreaseAccuracy = perror) ## return the whole matrix (= nperm*ntree values per variable)
    return(MeanDecreaseAccuracy = colMeans(perror)) ## return only averages over permutations and trees
}


varimpsurv <- function (object, mincriterion = 0, conditional = FALSE, 
                        threshold = 0.2, nperm = 1, OOB = TRUE, pre1.0_0 = conditional)
{

    cat("\n")
    cat("Variable importance for survival forests; this feature is _experimental_\n\n")
    response <- object@responses
    input <- object@data@get("input")
    xnames <- colnames(input)
    inp <- initVariableFrame(input, trafo = NULL)
    y <- object@responses@variables[[1]]
    weights <- object@initweights
    stopifnot(inherits(y, "Surv"))

    if (conditional || pre1.0_0) {
        if(!all(complete.cases(inp@variables)))
            stop("cannot compute variable importance measure with missing values")
    }
    stopifnot(requireNamespace("ipred", quietly = TRUE))
    error <- function(x, oob) ipred::sbrier(y[oob,,drop = FALSE], x[oob])

    pred <- function(tree, newinp, j = -1L) {

        where <- R_get_nodeID(tree, inp, mincriterion)
        wh <- .Call("R_get_nodeID", tree, newinp, mincriterion, as.integer(j), PACKAGE = "party")
        swh <- sort(unique(wh))
        RET <- vector(mode = "list", length = length(wh))
        for (i in 1:length(swh)) {
            w <- weights * (where == swh[i])
            RET[wh == swh[i]] <- list(mysurvfit(y, weights = w))
        }
        return(RET)
    }

    w <- object@initweights
    if (max(abs(w - 1)) > sqrt(.Machine$double.eps))
        warning(sQuote("varimp"), " with non-unity weights might give misleading results")

    ## list for several permutations
    perror <- matrix(0, nrow = nperm*length(object@ensemble), ncol = length(xnames))
    ## this matrix is initialized with values 0 so that a tree that does not 
    ## contain the current variable adds importance 0 to its average importance
    colnames(perror) <- xnames
        for (b in 1:length(object@ensemble)){
            tree <- object@ensemble[[b]]


            ## if OOB == TRUE use only oob observations, otherwise use all observations in learning sample
            if(OOB){oob <- object@weights[[b]] == 0} else{ oob <- rep(TRUE, length(y))}
            p <- pred(tree, inp)
            eoob <- error(p, oob)

            ## for all variables (j = 1 ... number of variables) 
            for(j in unique(varIDs(tree))){
              for (per in 1:nperm){
                 if (conditional || pre1.0_0) {
                    tmp <- inp
                    ccl <- create_cond_list(conditional, threshold, xnames[j], input)
                    if (is.null(ccl)) {
                        perm <- sample(which(oob))
                    } else {
                        perm <- conditional_perm(ccl, xnames, input, tree, oob)
                    }
                    tmp@variables[[j]][which(oob)] <- tmp@variables[[j]][perm]
                    p <- pred(tree, tmp, -1L)
                } else {
                    p <- pred(tree, inp, as.integer(j))
                }

                ## run through all rows of perror
                perror[(per+(b-1)*nperm), j] <- (error(p, oob) - eoob)

              } ## end of for (per in 1:nperm)
            } ## end of for(j in unique(varIDs(tree)))
        } ## end of for (b in 1:length(object@ensemble))

    perror <- as.data.frame(perror)
    #return(MeanDecreaseAccuracy = perror) ## return the whole matrix (= nperm*ntree values per variable)
    return(MeanDecreaseAccuracy = colMeans(perror)) ## return only averages over permutations and trees
}




# cutpoints_list() returns:
# - vector of cutpoints (length=number of cutpoints) 
#   if variable is continuous
# - vector of indicators (length=number of categories x number of cutpoints)
#   if variable is categorical (nominal or ordered)
cutpoints_list <- function(tree, variableID) {

    cutp <- function(node) {
       if (node[[4]]) return(NULL)
       cp <- NULL
       if (node[[5]][[1]] == variableID)
           cp <- node[[5]][[3]]
       nl <- cutp(node[[8]])
       nr <- cutp(node[[9]])
       return(c(cp, nl, nr))
    }
    return(cutp(tree))
}


conditional_perm <- function(cond, xnames, input, tree, oob){

    ## get cutpoints of all conditioning variables of the current variable of interest 
    ## and generate design matrix for permutation from factors in help
    blocks <- vector(mode = "list", length = length(cond))
                    
    for (i in 1:length(cond)) {

        ## varID is variable index or column number of input (predictor matrix) 
        ## not variable name!
        varID <- which(xnames == cond[i])


        ## if conditioning variable is not used for splitting in current tree
        ## proceed with next conditioning variable
        cl <- cutpoints_list(tree, varID)
        if (is.null(cl)) next

        ## proceed cutpoints for different types of variables
        x <- input[, varID]
        xclass <- class(x)[1]
        if (xclass == "integer") xclass <- "numeric"

        block <- switch(xclass, "numeric" = cut(x, breaks = c(-Inf, sort(unique(cl)), Inf)),
                        "ordered" = cut(as.numeric(x), breaks =  c(-Inf, sort(unique(cl)), Inf)),
                        "factor" = {
                            CL <- matrix(as.logical(cl), nrow = nlevels(x))                            
                            rs <- rowSums(CL)
                            dlev <- (1:nrow(CL))[rs %in% rs[duplicated(rs)]]
                            fuse <- c()
                            for (ii in dlev) {
                                for (j in dlev[dlev > ii]) {
                                    if (all(CL[ii,] == CL[j,])) fuse <- rbind(fuse, c(ii, j))
                                }
                            }
                            xlev <- 1:nlevels(x)
                            newl <- nlevels(x) + 1
                            block <- as.integer(x)
                            for (l in xlev) {
                                if (NROW(fuse) == 0) break
                                if (any(fuse[, 1] == l)) {
                                    f <- c(l, fuse[fuse[, 1] == l, 2])
                                    fuse <- fuse[!fuse[,1] %in% f, , drop = FALSE]
                                    block[block %in% f] <- newl
                                    newl <- newl + 1
                                 }
                            }
                            as.factor(block)
                         })
         blocks[[i]] <- block
    }

    ## remove non-splitting variables
    names(blocks) <- cond
    blocks <- blocks[!sapply(blocks, is.null)]

    ## if none of the conditioning variables are used in the tree
    if (!length(blocks)>0){
        perm <- sample(which(oob))
        return(perm)
    } else {
        blocks <- as.data.frame(blocks)
        ## from factors blocks create design matrix
        f <- paste("~ - 1 + ", paste(colnames(blocks), collapse = ":", sep = ""))
        des <- model.matrix(as.formula(f), data = blocks)

        ## one conditional permutation
        perm <- 1:nrow(input)
        for (l in 1:ncol(des)) {
           index <- which(des[,l] > 0 & oob)
           if (length(index) > 1)
               perm[index] <- sample(index)
           }
        return(perm[oob])
    }
}

varimpAUC <- function(object, mincriterion = 0, conditional = FALSE, 
                      threshold = 0.2, nperm = 1, OOB = TRUE, pre1.0_0 = conditional)
{

    response <- object@responses
    input <- object@data@get("input")
    xnames <- colnames(input)
    inp <- initVariableFrame(input, trafo = NULL)
    y <- object@responses@variables[[1]]
    if(length(response@variables) != 1)
        stop("cannot compute variable importance measure for multivariate response")

    if (conditional || pre1.0_0) {
        if(!all(complete.cases(inp@variables)))
            stop("cannot compute variable importance measure with missing values")
    }
    CLASS <- all(response@is_nominal)
    ORDERED <- all(response@is_ordinal)
    if (CLASS) {      
          if (nlevels(y)>2) {
            warning("AUC=TRUE works only for binary y\n error rate is used instead of AUC")
            error <- function(x, oob) mean((levels(y)[sapply(x, which.max)] != y)[oob])
          }   
          else {
             error <- function(x, oob) {
               xoob <- sapply(x, function(x) x[1])[oob]
               yoob <- y[oob]
               which1 <- which(yoob==levels(y)[1])
               noob1 <- length(which1)
               noob <- length(yoob)
               if (noob1==0|noob1==noob) { return(NA) }       # AUC cannot be computed if all OOB-observations are from one class
               return(1-sum(kronecker(xoob[which1] , xoob[-which1],">"))/(noob1*(length(yoob)-noob1)))       # calculate AUC
            }
       }
       ###  stop
    }
    else {
        if (ORDERED) {
            error <- function(x, oob) mean((sapply(x, which.max) != 
                y)[oob])
        }
        else {
            error <- function(x, oob) mean((unlist(x) - y)[oob]^2)
        }
    }

    w <- object@initweights
    if (max(abs(w - 1)) > sqrt(.Machine$double.eps))
        warning(sQuote("varimp"), " with non-unity weights might give misleading results")

    perror <- matrix(0, nrow = nperm*length(object@ensemble), ncol = length(xnames))
    colnames(perror) <- xnames
        for (b in 1:length(object@ensemble)){
            tree <- object@ensemble[[b]]

            if(OOB){oob <- object@weights[[b]] == 0} else{ oob <- rep(TRUE, length(xnames))}
            p <- .Call("R_predict", tree, inp, mincriterion, -1L, PACKAGE = "party")
            eoob <- error(p, oob)

            for(j in unique(varIDs(tree))){
              for (per in 1:nperm){

                if (conditional || pre1.0_0) {
                    tmp <- inp
                    ccl <- create_cond_list(conditional, threshold, xnames[j], input)
                    if (is.null(ccl)) {
                        perm <- sample(which(oob))
                    } else {
                        perm <- conditional_perm(ccl, xnames, input, tree, oob)
                    }
                    tmp@variables[[j]][which(oob)] <- tmp@variables[[j]][perm]
                    p <- .Call("R_predict", tree, tmp, mincriterion, -1L,
                       PACKAGE = "party")
                } else {
                    p <- .Call("R_predict", tree, inp, mincriterion, as.integer(j),
                               PACKAGE = "party")
                }
                perror[(per+(b-1)*nperm), j] <- (error(p, oob) - eoob)

              } 
            } 
        } 

    perror <- as.data.frame(perror)
    return(MeanDecreaseAccuracy = colMeans(perror, na.rm = TRUE)) ## na.rm = TRUE because with AUC-perm. VIM NA values occur whenever a tree's OOB-observations are all from the same class
}

