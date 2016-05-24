"booteval.relimp" <-
function (bootrun, bty = "perc", level = 0.95, sort = FALSE, norank = FALSE, 
    nodiff = FALSE, typesel = c("lmg", "pmvd", "last", "first", 
        "betasq", "pratt", "genizi", "car")) 
{

    # Author and copyright holder: Ulrike Groemping

    # This routine is distributed under GPL version 2 or newer.
    # The text of this license can be found at http://www.gnu.org/copyleft/gpl.html.

    # function for simulating percentage contribution and rank with CIs
    # result (ergebnis) contains estimated percentages, estimated ranks, CIs for percentages, CIs for ranks

    # bootrun is the result of bootstrap runs from boot.relimp
    # bty is the type of bootstrap intervals, default BCa, 
    #       eventually intended that vector can be given like in package boot
    #       currently bty only takes one single bty
    #      (rank intervals always with bty percentile, 
    #       since BCa does not work properly with ranks and normal not reasonable)
    # level is the confidence level, scalar or vector can be given like in package boot (option conf for boot.ci)
    # sort (if TRUE) requests that output is sorted by size of relative importances
    # norank (if TRUE) requests suppression of results for ranks (although ranks have been bootstrapped)
    # nodiff (if TRUE) requests suppression of results for differences (although differences have been bootstrapped)
    # typesel is the character string, character vector or list of metric types for which evaluation is requested

    # error control
    if (!(is(bootrun, "relimplmboot"))) 
        stop("bootrun must be output from routine boot.relimp")
    if ("type" %in% names(sys.call(1))) 
        stop("type is not a valid option for booteval.relimp. \n You may have intended to use bty or typesel.")
    if (!bty %in% c("perc", "bca", "norm", "basic")) 
        stop("bty type MUST be one of ", "perc ", "bca ", "norm ", 
            "basic")
    if (!(min(level) >= 0.5 && max(level) < 1)) 
        stop("invalid confidence levels requested: ", paste(level, 
            collapse = " "))
    type <- bootrun@type
    rank <- bootrun@rank
    diff <- bootrun@diff
    rela <- bootrun@rela
    always <- bootrun@always
    wt <- bootrun@wt
    if (is.null(wt)) wt <- rep(1, nrow(bootrun@boot$data))

    nlev <- length(level)
    if (length(bootrun@groupdocu)>0) {
    groups <- bootrun@groupdocu[[2]]
    groupnames <- bootrun@groupdocu[[1]][which(sapply(groups, length)>1)]
    groups <- groups[which(sapply(groups, length)>1)]
    }
    else 
    {
    groups <- NULL
    groupnames <- NULL
    }
    
    # prepare output object by first providing estimates themselves

    ## relimplm-Objekt generated
    ## values will subsequently be corrected, since they are wrong in case of interactions
    ## because ngroups information is not available here
    ausgabe <- calc.relimp(bootrun@boot$data, weights=wt, type = type, diff = diff, 
        rank = rank, rela = rela, always=always, groups=groups, 
        groupnames=groupnames)
    g <- length(ausgabe@namen)-1
    if (length(ausgabe@groupdocu)>0) g <- length(ausgabe@groupdocu[[2]])

    # extend output object
    ausgabe <- new("relimplmbooteval", ausgabe)  ## change UG for 1.3
    ausgabe@level <- level
    ausgabe@nboot <- bootrun@nboot
    ausgabe@bty <- bty
    ausgabe@rank <- rank
    ausgabe@diff <- diff
    ausgabe@rela <- rela
    ausgabe@fixed <- bootrun@fixed
    #provide names for identifying statistics (ungrouped case)
    names <- bootrun@namen
    if (!is.null(always)) names <- setdiff(names, bootrun@alwaysnam)
    if (!is.null(groupnames)) names <- c(names[1],as.character(bootrun@groupdocu[[1]]))
    diffnam <- paste(names[2:(g+1)][nchoosek(g, 2)[1, ]], names[2:(g+1)][nchoosek(g,2)[2, ]], sep = "-")
    ausgabe@var.y.boot <- bootrun@boot$t[, 1]
    ausgabe@R2.boot <- bootrun@boot$t[, 2]
    ausgabe@R2.decomp.boot <- bootrun@boot$t[, 3]

    #assign names to elements from bootrun in order to be able to refer to them later
    #columns of matrices can be referred to by their colnames (in quotes in square brackets instead of index)
    #elements of vectors analogously
    zaehl <- 4
    bootnames <- c("var.y", "R2", "R2.decomp")
    typname <- ""
    for (a in c("lmg", "pmvd", "last", "first", "betasq", "pratt", "genizi", "car")) {
        if (a %in% type) {
            bootnames <- c(bootnames, paste(names[2:(g+1)], ".", a, sep = ""))
            if (a %in% typesel) {
                slot(ausgabe, a) <- bootrun@boot$t0[zaehl:(zaehl + g - 1)]
                slot(ausgabe, paste(a, "boot", sep = ".")) <- bootrun@boot$t[, 
                  zaehl:(zaehl + g - 1)]}
            zaehl <- zaehl + g
            if (rank) {
                bootnames <- c(bootnames, paste(names[2:(g+1)], ".", a, 
                  "rank", sep = ""))
                if (a %in% typesel) 
                  slot(ausgabe, paste(a, "rank", sep = ".")) <- bootrun@boot$t0[zaehl:(zaehl + g - 1)]
                  slot(ausgabe, paste(a, "rank", "boot", sep = ".")) <- bootrun@boot$t[, 
                    zaehl:(zaehl + g - 1)]
                zaehl <- zaehl + g
            }
            if (diff) {
                bootnames <- c(bootnames, paste(diffnam, ".", 
                  a, "diff", sep = ""))
                if (a %in% typesel) 
                  slot(ausgabe, paste(a, "diff", sep = ".")) <- bootrun@boot$t0[zaehl:(zaehl + g*(g - 1)/2-1)]
                  slot(ausgabe, paste(a, "diff", "boot", sep = ".")) <- matrix(bootrun@boot$t[, 
                    zaehl:(zaehl + g * (g - 1)/2 - 1)],1,g * (g - 1)/2)
                zaehl <- zaehl + g * (g - 1)/2
            }
            if (a %in% typesel && !typname[1] == "") 
                typname <- c(typname, a)
            if (a %in% typesel && typname[1] == "") 
                typname <- a
        }
    }

    ## omit always in case it was appended as numeric variable
    bootrun@boot$t<-bootrun@boot$t[,1:length(bootnames)]
    bootrun@boot$t0<-bootrun@boot$t0[1:length(bootnames)]

    colnames(bootrun@boot$t) <- bootnames
    names(bootrun@boot$t0) <- bootnames
    if (!is.null(bootrun@vcov)){
    colnames(bootrun@vcov) <- bootnames
    rownames(bootrun@vcov) <- bootnames}
    ntype <- length(typname)
    percentages <- bootrun@boot$t0[paste(names[2:(g+1)], ".", matrix(typname, 
        g, ntype, byrow = T), sep = "")]
    if (rank) 
        ranks <- bootrun@boot$t0[paste(names[2:(g+1)], ".", matrix(typname, 
            g, ntype, byrow = T), "rank", sep = "")]
    if (diff) 
        diffs <- bootrun@boot$t0[paste(diffnam, ".", matrix(typname, 
            g * (g - 1)/2, ntype, byrow = T), "diff", sep = "")]
    perclower <- matrix(0, nlev, ntype * g, dimnames = list(level, 
        names(percentages)))
    if (rank) 
        ranklower <- matrix(0, nlev, ntype * g, dimnames = list(level, 
            names(percentages)))
    if (diff) 
        difflower <- matrix(0, nlev, ntype * g * (g - 1)/2, dimnames = list(level, 
            names(diffs)))
    percupper <- matrix(0, nlev, ntype * g, dimnames = list(level, 
        names(percentages)))
    if (rank) 
        rankupper <- matrix(0, nlev, ntype * g, dimnames = list(level, 
            names(percentages)))
    if (diff) 
        diffupper <- matrix(0, nlev, ntype * g * (g - 1)/2, dimnames = list(level, 
            names(diffs)))

    #strategy: if all bootstrap samples have same value-> lower=upper=that value
    #otherwise: boot.ci
        #determine confidence intervals
        #percentages
     for (j in 1:length(percentages)) {
        var <- var(bootrun@boot$t[, names(percentages)[j]])
        perclower[, j] <- percentages[j]
        percupper[, j] <- percentages[j]
        if (var > 0) {
            temp <- boot.ci(bootrun@boot, index = names(percentages)[j], 
                type = bty, conf = level)
            if (bty %in% c("perc", "bca", "basic")) {
                eval(parse(text = paste("perclower[,j]<-temp$", 
                  bty, "[,4]", sep = ""), n = 1))
                eval(parse(text = paste("percupper[,j]<-temp$", 
                  bty, "[,5]", sep = ""), n = 1))
            }
            if (bty == "norm") {
                eval(parse(text = paste("perclower[,j]<-temp$", 
                  "normal", "[,2]", sep = ""), n = 1))
                eval(parse(text = paste("percupper[,j]<-temp$", 
                  "normal", "[,3]", sep = ""), n = 1))
            }
        }
    }
        #determine confidence intervals
        #ranks
     if (rank && !norank) {
        for (j in 1:length(ranks)) {
            var <- var(bootrun@boot$t[, names(ranks)[j]])
            ranklower[, j] <- ranks[j]
            rankupper[, j] <- ranks[j]
            if (var > 0) {
                temp <- boot.ci(bootrun@boot, index = names(ranks)[j], 
                  type = "perc", conf = level)
                ranklower[, j] <- temp$percent[, 4]
                rankupper[, j] <- temp$percent[, 5]
            }
        }
    }
        #determine confidence intervals
        #diffs
     if (diff) {
        for (j in 1:length(diffs)) {
            var <- var(bootrun@boot$t[, names(diffs)[j]])
            difflower[, j] <- diffs[j]
            diffupper[, j] <- diffs[j]
            if (var > 0) {
                temp <- boot.ci(bootrun@boot, index = names(diffs)[j], 
                  type = bty, conf = level)
                if (bty %in% c("perc", "bca", "student", "basic")) {
                  eval(parse(text = paste("difflower[,j]<-temp$", 
                    bty, "[,4]", sep = ""), n = 1))
                  eval(parse(text = paste("diffupper[,j]<-temp$", 
                    bty, "[,5]", sep = ""), n = 1))
                }
                if (bty == "norm") {
                  eval(parse(text = paste("difflower[,j]<-temp$", 
                    "normal", "[,2]", sep = ""), n = 1))
                  eval(parse(text = paste("diffupper[,j]<-temp$", 
                    "normal", "[,3]", sep = ""), n = 1))
                }
            }
        }
    }
    ausgabe@type <- typname

    #write confidence bounds to output object
    for (a in typname) {
        slot(ausgabe, paste(a, "lower", sep = ".")) <- matrix(perclower[, 
            paste(names[2:(g+1)], a, sep = ".")], nlev, g)
        if (rank && !norank) 
            slot(ausgabe, paste(a, "rank", "lower", sep = ".")) <- matrix(ranklower[, 
                paste(names[2:(g+1)], a, sep = ".")], nlev, g)
        if (diff && !nodiff) 
            slot(ausgabe, paste(a, "diff", "lower", sep = ".")) <- matrix(difflower[, 
                paste(diffnam, paste(a, "diff", sep = ""), sep = ".")], 
                nlev, g * (g - 1)/2)
        slot(ausgabe, paste(a, "upper", sep = ".")) <- matrix(percupper[, 
            paste(names[2:(g+1)], a, sep = ".")], nlev, g)
        if (rank && !norank) 
            slot(ausgabe, paste(a, "rank", "upper", sep = ".")) <- matrix(rankupper[, 
                paste(names[2:(g+1)], a, sep = ".")], nlev, g)
        if (diff && !nodiff) 
            slot(ausgabe, paste(a, "diff", "upper", sep = ".")) <- matrix(diffupper[, 
                paste(diffnam, paste(a, "diff", sep = ""), sep = ".")], 
                nlev, g * (g - 1)/2)
    }
    #show confidence intervals with rank marks
    #only possible, if ranks are bootstrapped
       #initialize (character) matrix for showing (sorted) results with confidence info
    if (rank && !norank) 
        mark <- matrix(rep("", (g * ntype + ntype - 1) * (3 * 
            nlev + 1)), g * ntype + ntype - 1, 3 * nlev + 1, 
            dimnames = list(rep("", g * ntype + ntype - 1), c("percentage", 
                rep(level, 3))))
    else mark <- matrix(rep(0, (g * ntype + ntype - 1) * (2 * 
        nlev + 1)), g * ntype + ntype - 1, 2 * nlev + 1, dimnames = list(rep("", 
        g * ntype + ntype - 1), c("percentage", rep(level, 2))))
    if (sort) 
        marksort <- mark
    for (aa in 1:ntype) {
        a <- typname[aa]
        percent <- slot(ausgabe, a)
        names(percent) <- names[2:(g+1)]
        sortiert <- sort(percent, decreasing = T, index = T)
        percsort <- sortiert$x
        names(percsort) <- names[sortiert$ix]
        cilower <- matrix(slot(ausgabe, paste(a, "lower", sep = ".")), 
            nlev, g)
        ciupper <- matrix(slot(ausgabe, paste(a, "upper", sep = ".")), 
            nlev, g)
        if (rank && !norank) {
            lower <- matrix(slot(ausgabe, paste(a, "rank", "lower", 
                sep = ".")), nlev, g)
            upper <- matrix(slot(ausgabe, paste(a, "rank", "upper", 
                sep = ".")), nlev, g)
            for (j in 1:g) {
                # j is the rank that might or might not be in the confidence interval 
                # for the k-th sorted X
             hilf <- matrix(rep("", g * nlev), g, nlev)
                for (k in 1:g) {
                    #k is the k-th X, i.e. the k-th entry in percent for this type
                   for (i in 1:nlev) {
                         #i is the confidence level index
                      if (j < lower[i, k] | j > upper[i, k]) 
                      hilf[k, i] <- "_"
                    else hilf[k, i] <- LETTERS[j]
                  } # loop i
                } # loop k
                # append latest letter to mark
                mark[((aa - 1) * (g + 1) + 1):(aa * (g + 1) - 
                  1), 2:(1 + nlev)] <- matrix(paste(mark[((aa - 
                  1) * (g + 1) + 1):(aa * (g + 1) - 1), 2:(1 + 
                  nlev)], matrix(hilf, g, nlev), sep = ""), g, 
                  nlev)
            } # loop j
            mark[((aa - 1) * (g + 1) + 1):(aa * (g + 1) - 1), 
                ] <- matrix(cbind(format(percent,scientific=F), mark[((aa - 1) * (g + 
                1) + 1):(aa * (g + 1) - 1), 2:(1 + nlev)], format(t(cilower),scientific=F), 
                format(t(ciupper),scientific=F)), g, 3 * nlev + 1)
        } #if rank and !norank

        if (!rank || norank) {
            mark[((aa - 1) * (g + 1) + 1):(aa * (g + 1) - 1), 
                ] <- matrix(cbind(percent, t(cilower), t(ciupper)), 
                g, 2 * nlev + 1)
        }
        rownames(mark)[((aa - 1) * (g + 1) + 1):(aa * (g + 1) - 
            1)] <- paste(names(percent), a, sep = ".")
        if (sort) {
            marksort[((aa - 1) * (g + 1) + 1):(aa * (g + 1) - 
                1), ] <- mark[((aa - 1) * (g + 1) + 1):(aa * 
                (g + 1) - 1), ][sortiert$ix, ]
            rownames(marksort)[((aa - 1) * (g + 1) + 1):(aa * 
                (g + 1) - 1)] <- rownames(mark[((aa - 1) * (g + 
                1) + 1):(aa * (g + 1) - 1), ])[sortiert$ix]
        }
    } # loop aa

    # reduce number of displayed digits in percentages
    if (rank && !norank) 
        mark[, c(1, (2 + nlev):(3 * nlev + 1))] <- format(round(as.numeric(mark[, 
            c(1, (2 + nlev):(3 * nlev + 1))]), digits=4), nsmall=4, scientific=FALSE)
    else mark <- format(round(mark, digits = 4), nsmall=4, scientific=FALSE) 
    if (sort && rank && !norank) 
        marksort[, c(1, (2 + nlev):(3 * nlev + 1))] <- format(round(as.numeric(marksort[, 
            c(1, (2 + nlev):(3 * nlev + 1))]), digits=4), nsmall=4, scientific=FALSE)
    if (sort && (!rank || norank)) 
        marksort <- format(round(marksort, digits = 4), nsmall=4, scientific=FALSE)
    ## eliminate unwanted zeros or NAs between metrics in case of more than one metric
    ## and output desired result
    if (sort) {
        if (ntype>1) {
           for (i in 1:(ntype-1)) 
             marksort[i*(g+1),]<-rep("",ncol(marksort))
           }
        ausgabe@mark <- marksort
        }
    else {
        if (ntype>1) {
        for (i in 1:(ntype-1)) 
           mark[i*(g+1),]<-rep("",ncol(mark))
           }
        ausgabe@mark <- mark
        }

    # differences
    if (diff && !nodiff) {
        mark <- matrix(rep("", (g * (g - 1) * ntype/2 + ntype - 
            1) * (3 * nlev + 1)), g * (g - 1) * ntype/2 + ntype - 
            1, 3 * nlev + 1, dimnames = list(rep("", g * (g - 
            1) * ntype/2 + ntype - 1), c("difference", rep(level, 
            3))))
        if (sort) 
            marksort <- matrix(rep("", (g * (g - 1) * ntype/2 + 
                ntype - 1) * (3 * nlev + 1)), g * (g - 1) * ntype/2 + 
                ntype - 1, 3 * nlev + 1, dimnames = list(rep("", 
                g * (g - 1) * ntype/2 + ntype - 1), c("difference", 
                rep(level, 3))))
        for (aa in 1:ntype) {
            a <- typname[aa]
            differ <- slot(ausgabe, paste(a, "diff", sep = "."))
            sortiert <- sort(abs(differ), decreasing = T, index = T)
            names(differ) <- paste(diffnam, a, sep = ".")
            difflower <- matrix(slot(ausgabe, paste(a, "diff", 
                "lower", sep = ".")), nlev, g * (g - 1)/2)
            diffupper <- matrix(slot(ausgabe, paste(a, "diff", 
                "upper", sep = ".")), nlev, g * (g - 1)/2)
            hilf <- matrix(rep("", g * (g - 1) * nlev/2), g * 
                (g - 1)/2, nlev)
            for (k in 1:(g * (g - 1)/2)) {
                   # k refers to the k-th difference
                 for (i in 1:nlev) {
                       #i is the confidence level index
                   if (0 < difflower[i, k] | 0 > diffupper[i, 
                    k]) 
                    hilf[k, i] <- " * "
                  else hilf[k, i] <- " "
                } # loop i
            } # loop k

           #append latest letter to mark 
            mark[((aa - 1) * (g * (g - 1)/2 + 1) + 1):(aa * (g * 
                (g - 1)/2 + 1) - 1), 2:(1 + nlev)] <- matrix(paste(mark[((aa - 
                1) * (g * (g - 1)/2 + 1) + 1):(aa * (g * (g - 
                1)/2 + 1) - 1), 2:(1 + nlev)], matrix(hilf, g * 
                (g - 1)/2, nlev), sep = ""), g * (g - 1)/2, nlev)
            mark[((aa - 1) * (g * (g - 1)/2 + 1) + 1):(aa * (g * 
                (g - 1)/2 + 1) - 1), ] <- matrix(cbind(differ, 
                mark[((aa - 1) * (g * (g - 1)/2 + 1) + 1):(aa * 
                  (g * (g - 1)/2 + 1) - 1), 2:(1 + nlev)], t(difflower), 
                t(diffupper)), g * (g - 1)/2, 3 * nlev + 1)
            rownames(mark)[((aa - 1) * (g * (g - 1)/2 + 1) + 
                1):(aa * (g * (g - 1)/2 + 1) - 1)] <- names(differ)
            if (sort) {
                marksort[((aa - 1) * (g * (g - 1)/2 + 1) + 1):(aa * 
                  (g * (g - 1)/2 + 1) - 1), ] <- mark[((aa - 
                  1) * (g * (g - 1)/2 + 1) + 1):(aa * (g * (g - 
                  1)/2 + 1) - 1), ][sortiert$ix, ]
                rownames(marksort)[((aa - 1) * (g * (g - 1)/2 + 
                  1) + 1):(aa * (g * (g - 1)/2 + 1) - 1)] <- names(differ)[sortiert$ix]
            }
        } # loop aa

        ##reduce number of displayed digits in percentage differences
        if (!sort) {
            mark[, c(1, (2 + nlev):(3 * nlev + 1))] <- format(round(as.numeric(mark[, 
                c(1, (2 + nlev):(3 * nlev + 1))]), digits=4), nsmall=4, scientific=FALSE)
            if (ntype>1) {
              for (i in 1:(ntype-1)) 
                 mark[i * (g * (g - 1) / 2 + 1),] <- rep("", ncol(mark))
               }
            ausgabe@markdiff <- mark
            }
        if (sort) {
            marksort[, c(1, (2 + nlev):(3 * nlev + 1))] <- format(round(as.numeric(marksort[, 
                c(1, (2 + nlev):(3 * nlev + 1))]), digits=4), nsmall=4, scientific=FALSE)
            if (ntype>1) {
              for (i in 1:(ntype-1)) 
                 marksort[i * (g * (g - 1) / 2 + 1),] <- rep("", ncol(marksort))
               }
            ausgabe@markdiff <- marksort
        }
    } # if diff and !nodiff

    # set correct options for printing in output object
    ausgabe@nobs <- bootrun@nobs
    ausgabe@diff <- diff && !nodiff
    ausgabe@rank <- rank && !norank
    ausgabe@sort <- sort
    ausgabe@est <- bootrun@boot$t0
    ausgabe@vcov <- bootrun@vcov
    return(ausgabe)
}
