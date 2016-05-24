bundle.census <-
function (m, loops = FALSE) 
{
    if (isTRUE(is.array(m)) == FALSE) 
        stop("'m' sholud be an array.")
    ifelse(isTRUE(is.null(dimnames(m)[1]) == TRUE | is.null(dimnames(m)[1][[1]]) == 
        TRUE) == TRUE, LBS <- 1:nrow(m), LBS <- dimnames(m)[[1]])
    lbs <- seq(LBS)
    if (isTRUE(is.na(dim(m)[3])) == FALSE) {
        if (isTRUE(is.null(dimnames(m)[[3]]) == TRUE) | isTRUE(any(duplicated(dimnames(m)[[3]]))) == 
            TRUE) 
            dimnames(m)[[3]] <- 1:dim(m)[3]
    }
    if (is.na(dim(m)[3]) == FALSE) {
        x <- transf(dichot(m)[, , 1], "matlist", labels = lbs, 
            prsep = ", ", lb2lb = TRUE)
        if (isTRUE(dim(m)[3] > 1L) == TRUE) {
            for (k in 2:dim(m)[3]) x <- append(x, transf(dichot(m)[, 
                , k], "matlist", labels = lbs, prsep = ", ", 
                lb2lb = TRUE))
            rm(k)
        }
    }
    else {
        x <- transf(dichot(m), "matlist", labels = lbs, prsep = ", ", 
            lb2lb = TRUE)
    }
    dfl <- data.frame(matrix(ncol = 2L, nrow = 0L))
    for (i in 1:length(x)) {
        dfl[i, 1] <- strsplit(x[i], ", ")[[1]][1]
        dfl[i, 2] <- strsplit(x[i], ", ")[[1]][2]
    }
    rm(i)
    DF <- data.frame(matrix(ncol = 2L, nrow = 0L))
    k <- 1
    for (i in 1:nrow(dfl)) {
        if (isTRUE(dfl[i, 1] != dfl[i, 2]) == TRUE) 
            DF[k, ] <- dfl[i, ]
        k <- k + 1L
    }
    rm(i)
    rm(k)
    out <- list()
    inn <- list()
    All <- list()
    for (i in 1:length(lbs)) {
        out[[i]] <- as.numeric(DF[which(DF[, 1] == as.numeric(lbs[i])), 
            2])
        inn[[i]] <- as.numeric(DF[which(DF[, 2] == as.numeric(lbs[i])), 
            1])
        All[[i]] <- c(out[[i]], inn[[i]])
    }
    rm(i)
    rm(DF)
    finn <- list()
    fout <- list()
    for (i in seq(lbs)) {
        finn[[i]] <- which(tabulate(inn[[i]]) == dim(m)[3])
        fout[[i]] <- which(tabulate(out[[i]]) == dim(m)[3])
    }
    rm(i)
    full <- list()
    for (i in seq(lbs)) {
        full[[i]] <- intersect(fout[[i]], finn[[i]])
    }
    rm(i)
    rm(finn, fout)
    asym <- list()
    for (i in seq(lbs)) {
        asym[[i]] <- which(tabulate(All[[i]]) == 1L)
    }
    rm(i)
    dobl <- list()
    dout <- list()
    for (i in seq(lbs)) {
        dobl[[i]] <- which(tabulate(All[[i]]) == 2L)
        dout[[i]] <- which(tabulate(out[[i]]) == 2L)
    }
    rm(i)
    rete <- list()
    for (i in 1:length(dobl)) {
        tmprte <- vector()
        for (j in 1:length(dobl[[i]])) {
            if (isTRUE(dobl[[i]][j] %in% which(tabulate(inn[[i]]) == 
                1L)) == TRUE && isTRUE(dobl[[i]][j] %in% which(tabulate(out[[i]]) == 
                1L)) == TRUE) 
                tmprte[length(tmprte) + 1L] <- dobl[[i]][j]
        }
        rm(j)
        rete[[i]] <- tmprte
    }
    rm(i)
    rm(tmprte)
    if (isTRUE(is.na(dim(m)[3])) == FALSE) {
        tmp <- list()
        tt <- vector()
        if (isTRUE(is.null(dimnames(m)[[3]])) == FALSE) {
            for (i in 1:length(dimnames(m)[[3]])) tmp[i] <- dimnames(m)[[3]][i]
            rm(i)
            for (i in 1:length(tmp)) tt[i] <- (strsplit(tmp[[i]], 
                "")[[1]][1])
            rm(i)
        }
        else {
            for (i in 1:dim(m)[3]) tt[i] <- tmp[i] <- i
            rm(i)
        }
        for (k in 1:length(tt)) {
            allr <- paste("all", tt[k], sep = "_")
            assign(allr, transf(m[, , k], "matlist", labels = lbs, 
                prsep = ", ", lb2lb = TRUE))
            tmp <- transf(m[, , k], "matlist", labels = lbs, 
                prsep = ", ", lb2lb = TRUE)
            tDF <- data.frame(matrix(ncol = 2L, nrow = 0L))
            for (i in 1:length(tmp)) {
                tDF[i, 1] <- strsplit(tmp[i], ", ")[[1]][1]
                tDF[i, 2] <- strsplit(tmp[i], ", ")[[1]][2]
            }
            rm(i)
            rm(tmp)
            oud <- list()
            ind <- list()
            ald <- list()
            for (i in 1:length(lbs)) {
                oud[[i]] <- as.numeric(tDF[which(tDF[, 1] == 
                  as.numeric(lbs[i])), 2])
                ind[[i]] <- as.numeric(tDF[which(tDF[, 2] == 
                  as.numeric(lbs[i])), 1])
                ald[[i]] <- c(oud[[i]], ind[[i]])
            }
            rm(i)
            assign(allr, ald)
            rm(oud, ind, ald, allr)
            rm(tDF)
        }
        rm(k)
    }
    else {
        tt <- "R"
        tmp <- x
    }
    if (isTRUE(is.na(dim(m)[3])) == FALSE) {
        for (k in 1:length(tt)) {
            tmpxchg <- list()
            xchr <- paste("xch", tt[k], sep = "_")
            for (i in 1:length(rete)) {
                allr <- paste("all", tt[k], sep = "_")
                tmpxchr <- vector()
                for (j in 1:length(rete[[i]])) {
                  if (isTRUE(rete[[i]][j] %in% which(tabulate(eval(as.name(allr))[[i]]) == 
                    1L)) == TRUE) 
                    tmpxchr[length(tmpxchr) + 1L] <- rete[[i]][j]
                }
                rm(j)
                tmpxchg[[i]] <- tmpxchr
            }
            rm(i)
            assign(xchr, tmpxchg)
        }
        rm(k)
        rm(tmpxchr)
        rm(tmpxchg)
        xchg <- list()
        length(xchg) <- length(lbs)
        for (k in 1:length(tt)) {
            xchr <- paste("xch", tt[k], sep = "_")
            for (i in 1:length(rete)) {
                vecr <- paste("vec", i, sep = "")
                tmpxch <- vector()
                tmpxch <- eval(as.name(xchr))[[i]]
                if (sum(tmpxch) > 0L) {
                  assign(vecr, tmpxch)
                  xchg[[i]] <- eval(as.name(vecr))
                  rm(vecr)
                }
                else {
                  rm(vecr)
                }
            }
            rm(i)
        }
        rm(k)
        rm(tmpxch)
        rm(xchr)
    }
    else {
        xchg <- logical(0)
    }
    recp <- list()
    if (isTRUE(is.na(dim(m)[3])) == FALSE) {
        for (i in 1:length(rete)) {
            recp[[i]] <- rete[[i]][which(!(rete[[i]] %in% xchg[[i]]))]
        }
        rm(i)
        attr(recp, "names") <- lbs
    }
    else {
        recp <- rete
    }
    if (isTRUE(is.na(dim(m)[3])) == FALSE) {
        Eout <- list()
        length(Eout) <- length(lbs)
        for (i in 1:length(dobl)) {
            tmpout <- vector()
            for (j in 1:length(dobl[[i]])) {
                if (isTRUE(dobl[[i]][j] %in% dout[[i]]) == TRUE) 
                  tmpout[length(tmpout) + 1L] <- dobl[[i]][j]
            }
            rm(j)
            Eout[[i]] <- tmpout
        }
        rm(tmpout)
        trpr <- list()
        tinn <- list()
        tout <- list()
        for (i in seq(lbs)) {
            trpr[[i]] <- which(tabulate(All[[i]]) > 2L)
            tinn[[i]] <- which(tabulate(inn[[i]]) > 2L)
            tout[[i]] <- which(tabulate(out[[i]]) > 2L)
        }
        rm(i)
        teinn <- list()
        teout <- list()
        length(teinn) <- length(lbs)
        length(teout) <- length(lbs)
        for (i in 1:length(trpr)) {
            tmpinn <- vector()
            tmpout <- vector()
            for (j in 1:length(trpr[[i]])) {
                if (isTRUE(trpr[[i]][j] %in% tinn[[i]]) == TRUE) 
                  tmpinn[length(tmpinn) + 1L] <- trpr[[i]][j]
                if (isTRUE(trpr[[i]][j] %in% tout[[i]]) == TRUE) 
                  tmpout[length(tmpout) + 1L] <- trpr[[i]][j]
            }
            rm(j)
            teinn[[i]] <- tmpinn
            teout[[i]] <- tmpout
        }
        rm(tmpinn, tmpout)
        TEinn <- list()
        TEout <- list()
        for (i in 1:length(trpr)) {
            tmpinn <- vector()
            tmpout <- vector()
            for (j in 1:length(trpr[[i]])) {
                if (isTRUE(!(teinn[[i]][j] %in% out[[i]])) == 
                  TRUE) 
                  tmpinn[length(tmpinn) + 1L] <- teinn[[i]][j]
                if (isTRUE(!(teout[[i]][j] %in% inn[[i]])) == 
                  TRUE) 
                  tmpout[length(tmpout) + 1L] <- teout[[i]][j]
            }
            rm(j)
            TEinn[[i]] <- tmpinn
            TEout[[i]] <- tmpout
        }
        rm(i)
        rm(tmpinn, tmpout)
    }
    else {
        TEinn <- TEout <- Eout <- logical(0)
    }
    mix <- list()
    if (isTRUE(is.na(dim(m)[3])) == FALSE) {
        for (i in 1:length(trpr)) {
            mix[[i]] <- trpr[[i]][which(!(trpr[[i]] %in% TEinn[[i]] | 
                trpr[[i]] %in% TEout[[i]]))]
        }
        rm(i)
        mixe <- list()
        for (i in 1:length(mix)) {
            mixe[[i]] <- mix[[i]][which(!(mix[[i]] %in% full[[i]]))]
        }
        rm(i)
    }
    else {
        mixe <- logical(0)
    }
    if (loops) {
        if (isTRUE(is.na(dim(m)[3])) == FALSE) {
            vec <- vector()
            for (i in 1:dim(m)[3]) {
                vec <- append(vec, sum(dichot(diag(m[, , i]), 
                  c = 1L)))
            }
            rm(i)
            lop <- sum(vec)
        }
        else {
            lop <- sum(dichot(diag(m), c = 1L))
        }
    }
    ifelse(isTRUE(loops == FALSE) == TRUE, {
        bc <- cbind(abs(choose(nrow(m), 2) - (choose(nrow(m), 
            2)) - (length(unlist(full))/2L + length(unlist(asym))/2L + 
            length(unlist(recp))/2L + length(unlist(xchg))/2L + 
            (length(unlist(Eout)) + length(stats::na.exclude(unlist(TEout)))) + 
            length(unlist(mixe))/2L)), (choose(nrow(m), 2)) - 
            (length(unlist(asym))/2L + length(unlist(recp))/2L + 
                length(unlist(xchg))/2L + (length(unlist(Eout)) + 
                length(stats::na.exclude(unlist(TEout)))) + length(unlist(mixe))/2L + 
                length(unlist(full))/2L), length(unlist(asym))/2L, 
            length(unlist(recp))/2L, (length(unlist(Eout)) + 
                length(stats::na.exclude(unlist(TEout)))), length(unlist(xchg))/2L, 
            length(unlist(mixe))/2L, length(unlist(full))/2L)
        colnames(bc) <- c("BUNDLES", "NULL", "ASYMM", "RECIP", 
            "T.ENTR", "T.EXCH", "MIXED", "FULL")
        rownames(bc) <- "TOTAL"
    }, {
        bc <- cbind(abs(choose(nrow(m), 2) - (choose(nrow(m), 
            2)) - (length(unlist(full))/2L + length(unlist(asym))/2L + 
            length(unlist(recp))/2L + length(unlist(xchg))/2L + 
            (length(unlist(Eout)) + length(stats::na.exclude(unlist(TEout)))) + 
            length(unlist(mixe))/2L)), (choose(nrow(m), 2)) - 
            (length(unlist(asym))/2L + length(unlist(recp))/2L + 
                length(unlist(xchg))/2L + (length(unlist(Eout)) + 
                length(stats::na.exclude(unlist(TEout)))) + length(unlist(mixe))/2L + 
                length(unlist(full))/2L), length(unlist(asym))/2L, 
            length(unlist(recp))/2L, (length(unlist(Eout)) + 
                length(stats::na.exclude(unlist(TEout)))), length(unlist(xchg))/2L, 
            length(unlist(mixe))/2L, length(unlist(full))/2L, 
            lop)
        colnames(bc) <- c("BUNDLES", "NULL", "ASYMM", "RECIP", 
            "T.ENTR", "T.EXCH", "MIXED", "FULL", "LOOP")
        rownames(bc) <- "TOTAL"
    })
    return(bc)
}
