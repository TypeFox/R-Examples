RVineStructureSelect <- function(data, familyset = NA, type = 0, selectioncrit = "AIC", indeptest = FALSE, 
                                 level = 0.05, trunclevel = NA, progress = FALSE,  weights = NA, rotations = TRUE) { 
    d <- n <- dim(data)[2]
    N <- dim(data)[1]
    
    ## sanity checks 
    if (type == 0) 
        type <- "RVine" else if (type == 1) 
            type <- "CVine"
    if (type != "RVine" & type != "CVine") 
        stop("Vine model not implemented.")
    if (N < 2) 
        stop("Number of observations has to be at least 2.")
    if (d < 3) 
        stop("Dimension has to be at least 3.")
    if (any(data > 1) || any(data < 0)) 
        stop("Data has to be in the interval [0,1].")
    if (!is.na(familyset[1])) {
        for (i in 1:length(familyset)) { 
            if (!(familyset[i] %in% c(0, 1:10, 13, 14, 16:20,
                                      23, 24, 26:30, 33, 34, 36:40,
                                      104, 114, 124, 134, 204, 214, 224, 234))) 
                stop("Copula family not implemented.")
        }
    }
    if (selectioncrit != "AIC" && selectioncrit != "BIC") 
        stop("Selection criterion not implemented.")
    if (level < 0 & level > 1) 
        stop("Significance level has to be between 0 and 1.")
    
    ## set variable names and trunclevel if not provided 
    if (is.null(colnames(data))) 
        colnames(data) <- paste("V", 1:n, sep = "")
    if (is.na(trunclevel)) 
        trunclevel <- d
    
    ## adjust familyset 
    if (trunclevel == 0) 
        familyset <- 0
    if (rotations) 
        familyset <- with_rotations(familyset)
    
    ## initialize object for results
    RVine <- list(Tree = NULL, Graph = NULL)
    
    ## estimation in first tree ----------------------------
    # find optimal tree
    g <- initializeFirstGraph(data, weights)
    MST <- findMaximumTauTree(g, mode = type) 
    # estimate pair-copulas
    VineTree <- fit.FirstTreeCopulas(MST, 
                                     data, 
                                     familyset,
                                     selectioncrit,
                                     indeptest, 
                                     level,
                                     weights = weights)
    # store results
    RVine$Tree[[1]] <- VineTree
    RVine$Graph[[1]] <- g
    oldVineGraph <- VineTree
    
    ## estimation in higher trees --------------------------
    for (i in 2:(n - 1)) {
        # only estimate pair-copulas if not truncated
        if (trunclevel == i - 1) 
            familyset <- 0
        # find optimal tree
        g <- buildNextGraph(VineTree, weights)
        MST <- findMaximumTauTree(g, mode = type) 
        # estimate pair-copulas
        VineTree <- fit.TreeCopulas(MST,
                                    VineTree, 
                                    familyset, 
                                    selectioncrit,
                                    indeptest,
                                    level, 
                                    progress, 
                                    weights = weights)
        # store results
        RVine$Tree[[i]] <- VineTree
        RVine$Graph[[i]] <- g
    }
    
    ## return results as 'RVineMatrix' object
    as.RVM(RVine)
}

## full graph with Kendall's tau as edge weights
initializeFirstGraph <- function(data.univ, weights) {
    
    # C = cor(data.univ,method='kendall')
    q <- dim(data.univ)[2]
#     C <- matrix(rep(1, q * q), ncol = q)
#     
#     for (i in 1:(q - 1)) {
#         for (j in (i + 1):q) {
#             tau <- fasttau(data.univ[, i], data.univ[, j], weights)
#             C[i, j] <- tau
#             C[j, i] <- tau
#         }
#     }
#     rownames(C) <- colnames(C) <- colnames(data.univ)
    C <- TauMatrix(data = data.univ, weights = weights)
    
    g <- graph_from_adjacency_matrix(C, mode = "lower", weighted = TRUE, diag = FALSE)
    E(g)$tau <- E(g)$weight
    E(g)$name <- paste(as_edgelist(g)[, 1], as_edgelist(g)[, 2], sep = ",")
    
    for (i in 1:gsize(g)) {
        E(g)$conditionedSet[[i]] <- ends(g, i, names = FALSE)
    }
    return(g)
}

## find maximum spanning tree/ first vine tree
findMaximumTauTree <- function(g, mode = "RVine") {
    
    if (mode == "RVine") {
        return(mst(g, weights = 1 - abs(E(g)$weight)))
    } else if (mode == "CVine") {
        M <- abs(as_adjacency_matrix(g, attr = "weight", sparse = 0))
        sumtaus <- rowSums(M)
        root <- which.max(sumtaus)
        
        Ecken <- ends(g, 1:gsize(g), names = FALSE)
        pos <- Ecken[, 2] == root | Ecken[, 1] == root
        
        MST <- delete_edges(g, E(g)[!pos])
        
        return(MST)
    }
}

# not required any longer? Use TauMatrix instead
fasttau <- function(x, y, weights = NA) {
    if (any(is.na(weights))) {
        m <- length(x)
        n <- length(y)
        if (m == 0 || n == 0) 
            stop("both 'x' and 'y' must be non-empty")
        if (m != n) 
            stop("'x' and 'y' must have the same length")
        out <- .C("ktau",
                  x = as.double(x),
                  y = as.double(y),
                  N = as.integer(n),
                  tau = as.double(0),
                  S = as.double(0),
                  D = as.double(0),
                  T = as.integer(0), 
                  U = as.integer(0),
                  V = as.integer(0), 
                  PACKAGE = "VineCopula")
        ktau <- out$tau
    } else {
        ktau <- TauMatrix(matrix(c(x, y), length(x), 2), weights)[2, 1]
    }
    return(ktau)
}

## fit pair-copulas for the first vine tree
fit.FirstTreeCopulas <- function(MST, data.univ, type, copulaSelectionBy, testForIndependence, testForIndependence.level, weights = NA) {
    
    d <- gsize(MST)
    
    parameterForACopula <- list()
    
    for (i in 1:d) {
        parameterForACopula[[i]] <- list()
        
        a <- ends(MST, i, names = FALSE)
        
        parameterForACopula[[i]]$zr1 <- data.univ[, a[1]]
        parameterForACopula[[i]]$zr2 <- data.univ[, a[2]]
        
        E(MST)[i]$Copula.Data.1 <- list(data.univ[, a[1]])
        E(MST)[i]$Copula.Data.2 <- list(data.univ[, a[2]])
        
        if (is.null(V(MST)[a[1]]$name)) {
            E(MST)[i]$Copula.CondName.1 <- a[1] 
        } else { 
            E(MST)[i]$Copula.CondName.1 <- V(MST)[a[1]]$name
        }
        
        if (is.null(V(MST)[a[2]]$name)) {
            E(MST)[i]$Copula.CondName.2 <- a[2] 
        } else {
            E(MST)[i]$Copula.CondName.2 <- V(MST)[a[2]]$name
        }
        
        if (is.null(V(MST)[a[1]]$name) || is.null(V(MST)[a[2]]$name)) {
            E(MST)[i]$Copula.Name <- paste(a[1], a[2], sep = " , ") 
        } else {
            E(MST)[i]$Copula.Name <- paste(V(MST)[a[1]]$name, 
                                           V(MST)[a[2]]$name,
                                           sep = " , ")
        }
    }
    
    outForACopula <- lapply(X = parameterForACopula,
                            FUN = wrapper_fit.ACopula,
                            type, copulaSelectionBy,
                            testForIndependence,
                            testForIndependence.level, 
                            weights)
    
    for (i in 1:d) {
        E(MST)$Copula.param[[i]] <- c(outForACopula[[i]]$par,
                                      outForACopula[[i]]$par2)
        E(MST)[i]$Copula.type <- outForACopula[[i]]$family
        E(MST)[i]$Copula.out <- list(outForACopula[[i]])
        
        E(MST)[i]$Copula.CondData.1 <- list(outForACopula[[i]]$CondOn.1)
        E(MST)[i]$Copula.CondData.2 <- list(outForACopula[[i]]$CondOn.2)
    }
    
    return(MST)
}

## fit pair-copulas for vine trees 2,...
fit.TreeCopulas <- function(MST, oldVineGraph, type, copulaSelectionBy, testForIndependence, testForIndependence.level, progress, weights = NA) {
    d <- gsize(MST)
    
    parameterForACopula <- list()
    
    for (i in 1:d) {
        parameterForACopula[[i]] <- list()
        
        con <- as.vector(ends(MST, i, names = FALSE))
        
        temp <- ends(oldVineGraph, con, names = FALSE)
        
        if ((temp[1, 1] == temp[2, 1]) || (temp[1, 2] == temp[2, 1])) {
            same <- temp[2, 1]
        } else {
            if ((temp[1, 1] == temp[2, 2]) || (temp[1, 2] == temp[2, 2])) {
                same <- temp[2, 2]
            }
        }
        
        # other1 <- temp[1, temp[1, ] != same]
        # other2 <- temp[2, temp[2, ] != same]
        
        if (temp[1, 1] == same) {
            zr1 <- E(oldVineGraph)[con[1]]$Copula.CondData.2
            n1 <- E(oldVineGraph)[con[1]]$Copula.CondName.2
        } else {
            zr1 <- E(oldVineGraph)[con[1]]$Copula.CondData.1
            n1 <- E(oldVineGraph)[con[1]]$Copula.CondName.1
        }
        
        if (temp[2, 1] == same) {
            zr2 <- E(oldVineGraph)[con[2]]$Copula.CondData.2
            n2 <- E(oldVineGraph)[con[2]]$Copula.CondName.2
        } else {
            zr2 <- E(oldVineGraph)[con[2]]$Copula.CondData.1
            n2 <- E(oldVineGraph)[con[2]]$Copula.CondName.1
        }
        
        if (is.list(zr1)) {
            zr1a <- as.vector(zr1[[1]])
            zr2a <- as.vector(zr2[[1]])
            n1a <- as.vector(n1[[1]])
            n2a <- as.vector(n2[[1]])
        } else {
            zr1a <- zr1
            zr2a <- zr2
            n1a <- n1
            n2a <- n2
        }
        
        if (progress == TRUE) 
            message(n1a, " + ", n2a, " --> ", E(MST)[i]$name)
        
        parameterForACopula[[i]]$zr1 <- zr1a
        parameterForACopula[[i]]$zr2 <- zr2a
        
        E(MST)[i]$Copula.Data.1 <- list(zr1a)
        E(MST)[i]$Copula.Data.2 <- list(zr2a)
        
        E(MST)[i]$Copula.CondName.1 <- n1a
        E(MST)[i]$Copula.CondName.2 <- n2a
    }
    
    outForACopula <- lapply(X = parameterForACopula,
                            FUN = wrapper_fit.ACopula,
                            type, copulaSelectionBy,
                            testForIndependence, 
                            testForIndependence.level, 
                            weights)
    
    for (i in 1:d) {
        E(MST)$Copula.param[[i]] <- c(outForACopula[[i]]$par, outForACopula[[i]]$par2)
        E(MST)[i]$Copula.type <- outForACopula[[i]]$family
        E(MST)[i]$Copula.out <- list(outForACopula[[i]])
        
        E(MST)[i]$Copula.CondData.1 <- list(outForACopula[[i]]$CondOn.1)
        E(MST)[i]$Copula.CondData.2 <- list(outForACopula[[i]]$CondOn.2)
    }
    
    return(MST)
}

## initialize graph for next vine tree (possible edges)
buildNextGraph <- function(oldVineGraph, weights = NA) {
    
    # EL <- as_edgelist(oldVineGraph)
    d <- gsize(oldVineGraph)
    
    
    g <- make_full_graph(d)
    V(g)$name <- E(oldVineGraph)$name
    V(g)$conditionedSet <- E(oldVineGraph)$conditionedSet
    
    if (!is.null(E(oldVineGraph)$conditioningSet)) {
        V(g)$conditioningSet <- E(oldVineGraph)$conditioningSet
    }
    
    for (i in 1:gsize(g)) {
        
        con <- as.vector(ends(g, i, names = FALSE))
        
        temp <- ends(oldVineGraph, con, names = FALSE)
        
        ## check for proximity condition
        ok <- FALSE
        if ((temp[1, 1] == temp[2, 1]) || (temp[1, 2] == temp[2, 1])) {
            ok <- TRUE
            same <- temp[2, 1]
        } else {
            if ((temp[1, 1] == temp[2, 2]) || (temp[1, 2] == temp[2, 2])) {
                ok <- TRUE
                same <- temp[2, 2]
            }
        }
        
        ## if proximity condition is fulfilled
        if (ok) {
            # other1 <- temp[1, temp[1, ] != same]
            # other2 <- temp[2, temp[2, ] != same]
            
            if (temp[1, 1] == same) {
                zr1 <- E(oldVineGraph)[con[1]]$Copula.CondData.2
            } else {
                zr1 <- E(oldVineGraph)[con[1]]$Copula.CondData.1
            }
            
            if (temp[2, 1] == same) {
                zr2 <- E(oldVineGraph)[con[2]]$Copula.CondData.2
            } else {
                zr2 <- E(oldVineGraph)[con[2]]$Copula.CondData.1
            }
            # print(is.list(zr1))
            if (is.list(zr1)) {
                zr1a <- as.vector(zr1[[1]])
                zr2a <- as.vector(zr2[[1]])
            } else {
                zr1a <- zr1
                zr2a <- zr2
            }
            keine_nas <- !(is.na(zr1a) | is.na(zr2a))
            # print(keine_nas) print(zr1a) E(g)[i]$weight = cor(x=zr1[keine_nas],y=zr2[keine_nas], method='kendall')
            E(g)[i]$weight <- TauMatrix(cbind(zr1a[keine_nas], zr2a[keine_nas]), weights)[2,1]
            
            name.node1 <- strsplit(V(g)[con[1]]$name, split = " *[,;] *")[[1]]
            name.node2 <- strsplit(V(g)[con[2]]$name, split = " *[,;] *")[[1]]
            
            ## conditioning set
            schnitt <- intersect(name.node1, name.node2)      
#             for (j in 1:length(name.node1)) {
#                 for (k in 1:length(name.node2)) {
#                     if (name.node1[j] == name.node2[k]) {
#                         schnitt <- c(schnitt, name.node1[j])
#                         name.node1[j] <- ""
#                         name.node2[k] <- ""
#                         break
#                     }
#                 }
#             }
            
            ## conditioned set
            differenz <- c(setdiff(name.node1, name.node2), setdiff(name.node2, name.node1))
#             for (j in 1:length(name.node1)) {
#                 if (name.node1[j] != "") {
#                     differenz <- c(differenz, name.node1[j])
#                 }
#             }
#             for (j in 1:length(name.node2)) {
#                 if (name.node2[j] != "") {
#                     differenz <- c(differenz, name.node2[j])
#                 }
#             }
            
            E(g)[i]$name <- paste(paste(differenz, collapse = ","),
                                  paste(schnitt, collapse = ","),
                                  sep = " ; ")
            
            if (is.list(V(g)[con[1]]$conditionedSet)) {
                l1 <- c(as.vector(V(g)[con[1]]$conditionedSet[[1]]),
                        as.vector(V(g)[con[1]]$conditioningSet[[1]]))
                l2 <- c(as.vector(V(g)[con[2]]$conditionedSet[[1]]),
                        as.vector(V(g)[con[2]]$conditioningSet[[1]]))
            } else {
                l1 <- c(V(g)[con[1]]$conditionedSet,
                        V(g)[con[1]]$conditioningSet)
                l2 <- c(V(g)[con[2]]$conditionedSet,
                        V(g)[con[2]]$conditioningSet)
            }
            out <- intern_SchnittDifferenz(l1, l2)
            
            suppressWarnings({
                E(g)$conditionedSet[i] <- list(out$differenz)
            })
            suppressWarnings({
                E(g)$conditioningSet[i] <- list(out$schnitt)
            })
        }
        
        E(g)[i]$todel <- !ok
    }
    
    E(g)$tau <- E(g)$weight
    
    g <- delete_edges(g, E(g)[E(g)$todel])
    
    return(g)
}

wrapper_fit.ACopula <- function(parameterForACopula, type, ...) {
    return(fit.ACopula(parameterForACopula$zr1, 
                       parameterForACopula$zr2,
                       type,
                       ...))
}

intern_SchnittDifferenz <- function(liste1, liste2) {
    out <- list()
    out$schnitt <- intersect(liste1, liste2)
    out$differenz <- c(setdiff(liste1, liste2), setdiff(liste2, liste1))
    
#     for (j in 1:length(liste1)) {
#         for (k in 1:length(liste2)) {
#             if (!is.na(liste2[k]) && liste1[j] == liste2[k]) {
#                 out$schnitt <- c(out$schnitt, liste1[j])
#                 liste1[j] <- NA
#                 liste2[k] <- NA
#                 break
#             }
#         }
#     }
#     
#     for (j in 1:length(liste1)) {
#         if (!is.na(liste1[j])) {
#             out$differenz <- c(out$differenz, liste1[j])
#         }
#     }
#     for (j in 1:length(liste2)) {
#         if (!is.na(liste2[j])) {
#             out$differenz <- c(out$differenz, liste2[j])
#         }
#     }
    
    return(out)
}

## bivariate copula selection
fit.ACopula <- function(u1, u2, familyset = NA, selectioncrit = "AIC", indeptest = FALSE, level = 0.05, weights = NA) {
    
    ## select family and estimate parameter(s) for the pair copula
    out <- BiCopSelect(u1, u2,
                       familyset,
                       selectioncrit,
                       indeptest,
                       level,
                       weights = weights,
                       rotations = FALSE)
    
    ## change rotation if family is not symmetric wrt the main diagonal
    if (out$family %in% c(23, 24, 26:30, 124, 224)) {
        out$family <- out$family + 10
    } else if (out$family %in% c(33, 34, 36:40, 134, 234)) {
        out$family <- out$family - 10
    }
    
    ## tawn copulas also have to change type
    if (out$family%/%100 == 1) {
        out$family <- out$family + 100
    } else if (out$family%/%200 == 1) {
        out$family <- out$family - 100
    }
    
    ## store pseudo-observations for estimation in next tree
    out$CondOn.1 <- .C("Hfunc1", 
                       as.integer(out$family),
                       as.integer(length(u1)),
                       as.double(u1),
                       as.double(u2),
                       as.double(out$par),
                       as.double(out$par2), 
                       as.double(rep(0, length(u1))), 
                       PACKAGE = "VineCopula")[[7]]
    out$CondOn.2 <- .C("Hfunc2",
                       as.integer(out$family),
                       as.integer(length(u1)),
                       as.double(u2),
                       as.double(u1), 
                       as.double(out$par), 
                       as.double(out$par2), 
                       as.double(rep(0, length(u1))), 
                       PACKAGE = "VineCopula")[[7]]
    
    ## return results
    out
}

## build R-Vine matrix object based on nested set of trees
as.RVM <- function(RVine) {
    
    ## initialize objects
    n <- length(RVine$Tree) + 1
    con <- list()
    nam <- V(RVine$Tree[[1]])$name
    conditionedSets <- NULL
    corresppondingParams <- list()
    corresppondingTypes <- list()
    
    ## get selected pairs, families and estimated parameters
    if (is.list(E(RVine$Tree[[n - 1]])$conditionedSet)) {
        conditionedSets[[n - 1]][[1]] <- (E(RVine$Tree[[n - 1]])$conditionedSet[[1]])
        for (k in 1:(n - 2)) {
            # conditionedSets[[k]] = E(RVine$Tree[[k]])$conditionedSet[[1]]
            conditionedSets[[k]] <- E(RVine$Tree[[k]])$conditionedSet
            corresppondingParams[[k]] <- as.list(E(RVine$Tree[[k]])$Copula.param)
            corresppondingTypes[[k]] <- as.list(E(RVine$Tree[[k]])$Copula.type)
        }
        
        corresppondingParams[[n - 1]] <- list()
        corresppondingParams[[n - 1]] <- as.list(E(RVine$Tree[[n - 1]])$Copula.param)
        corresppondingTypes[[n - 1]] <- as.list(E(RVine$Tree[[n - 1]])$Copula.type)
        # print(corresppondingParams)
    } else {
        conditionedSets[[n - 1]][[1]] <- (E(RVine$Tree[[n - 1]])$conditionedSet)
        for (k in 1:(n - 2)) {
            conditionedSets[[k]] <- E(RVine$Tree[[k]])$conditionedSet
            corresppondingParams[[k]] <- as.list(E(RVine$Tree[[k]])$Copula.param)
            corresppondingTypes[[k]] <- as.list(E(RVine$Tree[[k]])$Copula.type)
        }
        # print(conditionedSets)
        corresppondingParams[[n - 1]] <- list()
        corresppondingParams[[n - 1]] <- as.list(E(RVine$Tree[[n - 1]])$Copula.param)
        corresppondingTypes[[n - 1]] <- as.list(E(RVine$Tree[[n - 1]])$Copula.type)
    }
    
    ## initialize matrices for RVineMatrix object
    Param <- array(dim = c(n, n))
    Params2 <- array(0, dim = c(n, n))
    Type <- array(dim = c(n, n))
    M <- matrix(NA, n, n)
    
    ## store structure, families and parameters in matrices
    for (k in 1:(n - 1)) {
        w <- conditionedSets[[n - k]][[1]][1]
        
        M[k, k] <- w
        M[(k + 1), k] <- conditionedSets[[n - k]][[1]][2]
        
        Param[(k + 1), k] <- corresppondingParams[[n - k]][[1]][1]
        Params2[(k + 1), k] <- corresppondingParams[[n - k]][[1]][2]
        
        Type[(k + 1), k] <- corresppondingTypes[[n - k]][[1]]
        
        if (k == (n - 1)) {
            M[(k + 1), (k + 1)] <- conditionedSets[[n - k]][[1]][2]
        } else {
            for (i in (k + 2):n) {
                for (j in 1:length(conditionedSets[[n - i + 1]])) {
                    cs <- conditionedSets[[n - i + 1]][[j]]
                    cty <- corresppondingTypes[[n - i + 1]][[j]]
                    if (cs[1] == w) {
                        M[i, k] <- cs[2]
                        Type[i, k] <- cty
                        break
                    } else if (cs[2] == w) {
                        # correct family for rotations
                        if (cty %in% c(23, 24, 26:30, 124, 224)) {
                            cty <- cty + 10
                        } else if (cty %in% c(33, 34, 36:40, 134, 234)) {
                            cty <- cty - 10
                        }
                        # change type for Tawn
                        if (cty%/%100 == 1) {
                            cty <- cty + 100
                        } else if (cty%/%200 == 1) {
                            cty <- cty - 100
                        }
                        M[i, k] <- cs[1]
                        Type[i, k] <- cty
                        break
                    }
                }
                Param[i, k] <- corresppondingParams[[n - i + 1]][[j]][1]
                Params2[i, k] <- corresppondingParams[[n - i + 1]][[j]][2]
                conditionedSets[[n - i + 1]][[j]] <- NULL
                corresppondingParams[[n - i + 1]][[j]] <- NULL
                corresppondingTypes[[n - i + 1]][[j]] <- NULL
            }
        }
        
    }
    
    ## clean NAs
    M[is.na(M)] <- 0
    Type[is.na(Type)] <- 0
    
    ## return RVineMatrix object
    RVineMatrix(M, family = Type, par = Param, par2 = Params2, names = nam)
}
