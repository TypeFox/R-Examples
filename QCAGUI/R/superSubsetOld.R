`superSubsetOld` <-
function(mydata, outcome = "", conditions = c(""), relation = "", incl.cut = 1, cov.cut = 0,
         index = 1, neg.out = FALSE, use.tilde = FALSE) {
    
    if (!isNamespaceLoaded("QCA")) {
        requireNamespace("QCA", quietly = TRUE)
    }
    
    if (all(conditions == c(""))) {
        conditions <- names(mydata)[-which(names(mydata) == outcome)]
    }
    
    verify.tt(mydata, outcome, conditions)
    
    if (relation == "" | is.na(match(relation, c("necessity", "sufficiency")))) {
        relation <- "necessity"
    }
    else if (length(relation) > 1) {
        relation <- relation[1]
    }
    
    colnames(mydata) <- toupper(colnames(mydata))
    conditions <- toupper(conditions)
    outcome <- toupper(outcome)
    
    mydata <- mydata[, c(conditions, outcome)]
    nofconditions <- length(conditions)
    
    uplow <- !use.tilde
    
    fuzzy.cc <- apply(mydata[, conditions], 2, function(x) any(x %% 1 > 0))
    if (mv.data <- any(mydata[, conditions] > 1)) {
        uplow <- use.tilde <- FALSE
    }
            
    check.equal <- function(x, y) {
        check.vector <- as.logical(unlist(lapply(x, all.equal, y)))
        check.vector[is.na(check.vector)] <- FALSE
        return(check.vector)
    }
    
    findClosest <- function(start.line, noflevels, mbase) {
        closest <- vector(length=sum(noflevels))
        clindex <- 1
        for (i in seq(length(noflevels))) {
            for (j in seq(noflevels[i])) {
                closest[clindex] <- j*mbase[i] + 1
                clindex <- clindex + 1
            }
        }
        return(closest + start.line)
    }
    
    getInclusion <- function(line.number, relation, connector, fuzzy.cc) {
        minscores <- maxscores <- vector(length=nrow(mydata))
        for (i in seq(nrow(mydata))) {
            values <- as.numeric(mydata[i, conditions])
            x <- getRow(noflevels + 1, line.number)
            
            if (any(onex3k <- x[fuzzy.cc] == 1)) {
                values[fuzzy.cc][onex3k] <- 1 - values[fuzzy.cc][onex3k]
            }
            
            copy.values <- values[!fuzzy.cc]
            if (length(copy.values) > 0) {
                values[!fuzzy.cc][x[!fuzzy.cc] != copy.values + 1] <- 0
                values[!fuzzy.cc][x[!fuzzy.cc] == copy.values + 1] <- 1
            }
            minscores[i] <- min(values[x != 0])
            maxscores[i] <- max(values[x != 0])
        }
        
        if (relation == "sufficiency") {
                inclusion <- sum(pmin(minscores, mydata[, outcome]))/sum(minscores)
                coverage  <- sum(pmin(minscores, mydata[, outcome]))/sum(mydata[, outcome])
        }
        else if (relation == "necessity") {
                inclusion <- sum(pmin(minscores, mydata[, outcome]))/sum(mydata[, outcome])
                coverage  <- sum(pmin(minscores, mydata[, outcome]))/sum(minscores)
            if (connector == "or") {
                inclusion <- sum(pmin(maxscores, mydata[, outcome]))/sum(mydata[, outcome])
                coverage  <- sum(pmin(maxscores, mydata[, outcome]))/sum(maxscores)
            }
        }
        
        inclusion[is.na(inclusion)] <- NA
        coverage[is.na(coverage)] <- NA
        return(c(inclusion, coverage))
    }
    
    extract <- function(closest, complex.list, noflevels.local, mbase.local, conditions.local, relation, connector, fuzzy.cc, uplow, use.tilde) {
        for (i in seq(length(closest))) {
            complex.list[[i]] <- vector("list", 2)
            row3k <- getRow(noflevels + 1, closest[i])
            colnames(row3k) <- conditions
            
            incov <- getInclusion(closest[i], relation, connector, fuzzy.cc)
            term <- QCA::writePrimeimp(row3k, uplow=uplow, use.tilde=use.tilde)
            if (relation == "necessity" & connector == "or") {
                term <- gsub("\\*", "+", term)
            }
            
            complex.list[[i]][[1]] <- c(row3k, tline=closest[i], term, inclusion=incov[1], coverage=incov[2])
            if (!is.na(incov[1])) {
                
                relation.inclusion <- FALSE
                if (relation == "sufficiency" | (relation == "necessity" & connector == "or")) {
                    relation.inclusion <- (incov[1] < incl.cut) & (length(noflevels.local) > 1)
                }
                else if (relation == "necessity") {
                    relation.inclusion <- (incov[1] > incl.cut | check.equal(incov[1], incl.cut)) & length(noflevels.local) > 1
                }
                
                if (relation.inclusion) {
                    noflevels.index <- which(cumsum(noflevels.local) >= i)[1]
                    next.closest <- findClosest(closest[i] - 1, noflevels.local[-noflevels.index], mbase.local[-noflevels.index])
                    match.vector <- match(next.closest, term.vector)
                    next.closest <- next.closest[!is.na(match.vector)]
                    match.vector <- match.vector[!is.na(match.vector)]
                    if (any(match.vector)) {
                        term.vector <<- term.vector[-match.vector]
                        complex.list[[i]][[2]] <- vector("list", length(next.closest))
                        complex.list[[i]][[2]] <- Recall(next.closest,
                                                         complex.list[[i]][[2]],
                                                         noflevels.local[-noflevels.index],
                                                         mbase.local[-noflevels.index],
                                                         conditions.local[-noflevels.index],
                                                         relation,
                                                         connector,
                                                         fuzzy.cc,
                                                         uplow,
                                                         use.tilde)
                    }
                }
                else {
                    complex.list[[i]] <- complex.list[[i]][-2]
                }
            }
            else {
                complex.list[[i]] <- complex.list[[i]][-2]
            }
        }
        return(complex.list)
    }
    
    if (neg.out) {
        mydata[, outcome] <- 1 - mydata[, outcome]
    }
    
    noflevels <- rep(2, nofconditions)
    mbase <- rev(c(1, cumprod(rev(noflevels + 1))))[-1]
        
    mostmin <- findClosest(0, noflevels, mbase)
    term.vector <- seq_len(prod(noflevels + 1))
    connector <- "and"
    result <- as.data.frame(matrix(unlist(extract(mostmin,
                                                  vector("list", length(mostmin)),
                                                  noflevels,
                                                  mbase,
                                                  conditions,
                                                  relation,
                                                  connector,
                                                  fuzzy.cc,
                                                  uplow,
                                                  use.tilde)), ncol=length(conditions) + 4, byrow=T))
    result[, ncol(result) - 1] <- as.numeric(as.character(result[, ncol(result) - 1]))
    result <- result[!is.na(result[, ncol(result) - 1]), ]
    result <- result[result[, ncol(result) - 1] > incl.cut | check.equal(result[, ncol(result) - 1], incl.cut), ]
    
    if (nrow(result) == 0) {
        connector <- "or"
        result <- as.data.frame(matrix(unlist(extract(mostmin,
                                                      vector("list", length(mostmin)),
                                                      noflevels,
                                                      mbase,
                                                      conditions,
                                                      relation,
                                                      connector,
                                                      fuzzy.cc,
                                                      uplow,
                                                      use.tilde)), ncol=length(conditions) + 4, byrow=T))
        result <- result[!is.na(result[, ncol(result) - 1]), ]
        result[, ncol(result) - 1] <- as.numeric(as.character(result[, ncol(result) - 1]))
        result <- result[result[, ncol(result) - 1] > incl.cut | check.equal(result[, ncol(result) - 1], incl.cut), ]
    }
    
    for (i in seq(ncol(result))) {
        result[, i] <- as.character(result[, i])
    }
    
    names(result) <- c(conditions, "tline", "grouping", "incl", "cov")
    result$tline <- as.numeric(result$tline)
    result$incl <- as.numeric(result$incl)
    result$cov <- as.numeric(result$cov)
    
    index <- 1
    while (index <= nrow(result)) {
        subs <- findSupersets(noflevels + 1, result$tline[index])
        subs.match <- match(subs, result$tline[-index])
        subs.match <- subs.match[!is.na(subs.match)]
        if (any(subs.match)) {
            result <- result[-index, ]
        }
        else {
            index <- index + 1
        }
    }
    
    result <- result[result$cov > cov.cut, ]
    
    if (nrow(result) > 0) {
        rownames(result) <- seq(nrow(result))
        return(structure(list(result=result[, seq(ncol(result) - 2, ncol(result))]), class="ss"))
    }
    else {
        cat("\nThere are no groupings that match given criteria.\n\n")
    }
}

