"mianalyze.relimp" <-
function (implist, level = 0.95, sort = FALSE, ..., 
    b = 50, type = "lmg", diff = TRUE, no.CI = FALSE, 
    rela = FALSE, always = NULL, groups = NULL, groupnames = NULL,
    deslist = NULL, bootlist.out = FALSE, formula=NULL,weights=NULL, strata=NULL, ids=NULL)
{
    # Author and copyright holder: Ulrike Groemping

    # This routine is distributed under GPL version 2 or newer.
    # The text of this license can be found at http://www.gnu.org/copyleft/gpl.html.

    # function for simulating percentage contribution and rank with CIs
    # result (ergebnis) contains estimated percentages, estimated ranks, CIs for percentages, CIs for ranks

    # implist is a list of imputed data frames (or matrices) or an imputationList object
    # level is the confidence level, scalar or vector can be given like in package boot (option conf for boot.ci)
    # b is number of bootstrap replicates
    # diff (if TRUE) requests investigation of pairwise differences
    # rela (if TRUE) requests forcing to sum 100%
    # type is a character vector or list (or single character string) of relative importance 
    #    types for which bootstrapping is requested
    # always gives regressors to adjust for
    # groups covers grouping of regressors
    # groupnames allows user-defined group labels instead of default G1 etc.
    # deslist is a list of survey-designs corresponding to the data frame list implist
    #      if given, it needs to be of the same length as implist; 
    #      otherwise, a list of unstructured designs reflecting simple random sampling is constructed here
    # bootlist.out indicates whether bootstrap results for each imputed data set are to be output
    # formula is a model formula for the linear model
    # weights is a set of weights used in constructing a design with function svydesign
    # strata is an expression for strata definition as needed in function svydesign
    # ids is an expression for Cluster ID definition as needed in function svydesign
    
    rank <- FALSE  ### rank = TRUE nicht unterstuetzt
    fixed <- FALSE  ### fixed = TRUE nicht unterstuetzt
    b <- b
    
    # error control
    if (!(is.list(implist) || is.data.frame(implist) || is.matrix(implist) || "imputationList" %in% class(implist)))
        stop(paste("implist has to be a list of imputed dataframes or matrices","\n",
           "or an object of class imputationList.",sep=""))
    if ("imputationList" %in% class(implist)) implist <- implist[[1]]   ## extract list of data frames etc.
    if (length(implist)<2)
        stop("You need AT LEAST two imputed datasets - and most likely more would be reasonable.")
    m <- length(implist)
    
    hilfdat <- sapply(implist,FUN="is.data.frame")
    hilfmat <- sapply(implist,FUN="is.matrix")
    if (!(all(hilfdat) || all(hilfmat))) stop("ERROR: implist must be a list of data frames or a list of matrices.")
    
    hilf <- sapply(implist,FUN="ncol")
    if (any(sapply(hilf[2:length(implist)],FUN=function(obj) !hilf[1]==obj)))
       stop ("ERROR: Imputed data sets in implist are of different length.")
    n <- nrow(implist[[1]])
    
    hilf <- sapply(implist,FUN="ncol")
    if (any(sapply(hilf[2:length(implist)],FUN=function(obj) hilf[1]!=obj)))
       stop ("ERROR: Imputed data sets in implist have a different number of columns.")
    hilf <- lapply(implist,FUN="colnames")
    if (any(sapply (hilf[2:length(implist)], function(obj) obj!=hilf[[1]])))
       stop ("ERROR: Imputed data sets in implist have different column names. This is not permitted.")

    if (!is.null(deslist) && (!length(implist)==length(deslist) || "survey.design" %in% class(deslist)))
        stop(paste("If you specify designs (deslist)","\n",
           "they must come as a list as long as the list in implist.",sep=""))
    if (!is.null(deslist) && !(is.null(weights) && is.null(strata) && is.null(ids)))
        stop("When specifying a deslist, weights, strata and ids MUST be NULL.")
    if ((!is.null(weights)) && !length(weights)==n) 
        stop("weights must be a numeric vector with an entry for each data row.")
    if ((!is.null(weights)) && !is.numeric(weights)) 
        stop("weights must be a numeric vector.")
    if (is.null(deslist)){
      if (is.null(weights)) weights <- rep(1/n,n)
      if (is.null(strata)) strata <- rep(1,n)
    }

    if (!(is.null(formula) || is(formula,"formula"))) stop ("formula must be NULL or a valid formula.")
    if (!(min(level) >= 0.5 && max(level) < 1))
        stop("invalid confidence levels requested: ", paste(level,
            collapse = " "))

## nur calc.relimp durchfuehren und ein Objekt des Typs relimplm zurueckgeben
if (no.CI){
  if (is.null(formula)) 
        if (!is.null(deslist)) 
        hilflist <- lapply(implist, function(obj){ calc.relimp(obj, design=deslist[[1]], 
            type = type, diff = diff, rank=rank, 
            rela = rela, always = always, groups = groups, groupnames = groupnames)})
        else 
        hilflist <- lapply(implist, function(obj){ calc.relimp(obj, weights=weights, 
                  type = type, diff = diff, rank=rank, 
                  rela = rela, always = always, groups = groups, groupnames = groupnames)})
  else 
        if (!is.null(deslist)) 
        hilflist <- lapply(implist, function(obj){ calc.relimp(formula, data=obj, design=deslist[[1]], 
            type = type, diff = diff, rank=rank, 
            rela = rela, always = always, groups = groups, groupnames = groupnames)})
        else {
        wt <- weights
        hilflist <- lapply(implist, function(obj){ calc.relimp(formula, data=obj, weights=wt, 
                  type = type, diff = diff, rank=rank, 
                  rela = rela, always = always, groups = groups, groupnames = groupnames)})
        }
   ausgabe <- hilflist[[1]]
   ausgabe@R2 <- mean(sapply(hilflist,function(obj){obj@R2}))
   ausgabe@R2.decomp <- mean(sapply(hilflist,function(obj){obj@R2}))
   ausgabe@var.y <- mean(sapply(hilflist,function(obj){obj@R2}))
   for (a in c("lmg", "pmvd", "last", "first", "betasq", "pratt")) {
      if (a %in% type) {
   slot(ausgabe,a) <- apply(sapply(hilflist,function(obj){slot(obj,a)}),1,"mean")
   if (diff) {
       if (length(slot(hilflist[[1]],paste(a,".diff",sep="")))==1)
       slot(ausgabe,paste(a,".diff",sep="")) <- 
            mean(sapply(hilflist,function(obj){slot(obj,paste(a,".diff",sep=""))})) 
       else
       slot(ausgabe,paste(a,".diff",sep="")) <- 
       apply(sapply(hilflist,function(obj){slot(obj,paste(a,".diff",sep=""))}),1,"mean") 
       }
   }
   }
   return(ausgabe)
   stop()
}

#### create bootstrap objects
bootlist <- list(m)
if (is.null(formula)) 
  for (i in 1:m) {
      if (!is.null(deslist)) 
      bootlist[[i]] <- boot.relimp(implist[[i]], design=deslist[[i]], 
          b = b, type = type, diff = diff, rank=rank, 
          rela = rela, always = always, groups = groups, groupnames = groupnames)
      else {
      if (is.null(ids)) 
      bootlist[[i]] <- boot.relimp(implist[[i]], design=svydesign(~1, data=implist[[i]], 
                                                        strata=~strata,weights=~weights), 
                b = b, type = type, diff = diff, rank=rank, 
                rela = rela, always = always, groups = groups, groupnames = groupnames)
      else
      bootlist[[i]] <- boot.relimp(implist[[i]], design=svydesign(ids=~ids, data=implist[[i]], 
                                                        strata=~strata,weights=~weights), 
                b = b, type = type, diff = diff, rank=rank, 
                rela = rela, always = always, groups = groups, groupnames = groupnames)
      }
      if (i>1 && !bootlist.out) bootlist[[i]]@boot$data <- NULL ## save space, first is needed further down
      }
else 
  for (i in 1:m) {
      if (!is.null(deslist)) 
      bootlist[[i]] <- boot.relimp(formula,implist[[i]], design=deslist[[i]], 
          b = b, type = type, diff = diff, rank=rank, 
          rela = rela, always = always, groups = groups, groupnames = groupnames)
      else {
      if (is.null(ids)) 
      bootlist[[i]] <- boot.relimp(formula,implist[[i]], design=svydesign(~1, data=implist[[i]], strata=~strata,weights=~weights), 
                b = b, type = type, diff = diff, rank=rank, 
                rela = rela, always = always, groups = groups, groupnames = groupnames)
      else
      bootlist[[i]] <- boot.relimp(formula,implist[[i]], design=svydesign(ids=~ids, data=implist[[i]], strata=~strata,weights=~weights), 
                b = b, type = type, diff = diff, rank=rank, 
                rela = rela, always = always, groups = groups, groupnames = groupnames)
      }
      if (i>1 && !bootlist.out) bootlist[[i]]@boot$data <- NULL ## save space, first is needed further down
      }
      if (!bootlist.out)
      bootlist[[i]]@boot$t <- bootlist[[i]]@boot$t0  ## save space

#### anschließend MIcombine fuer die Berechnungen der Kombinationen
#### Ausgabe-Verhalten aus diesem Programm (ohne den Kram fuer Raenge)

bootrun <- bootlist[[1]]  ## for checks on groups etc., since should be same for all list elements
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

    ## relimplm-Objekt generiert
    ## Werte werden anschließend unten korrigiert, da im Fall von WWen falsch
    ## wg. nichtvorhandener ngroups-Information hier
    ausgabe <- calc.relimp(bootrun@boot$data, design=deslist[[1]], type = type, diff = diff,
        rank = rank, rela = rela, always=always, groups=groups,
        groupnames=groupnames)
    g <- length(ausgabe@namen)-1
    if (length(ausgabe@groupdocu)>0) g <- length(ausgabe@groupdocu[[2]])

    # extend output object
    ausgabe <- new("relimplmbootMI", ausgabe)  ## change UG for 1.3
    ausgabe@level <- level
    ausgabe@nboot <- bootrun@nboot
    ausgabe@bty <- "MI"
    ausgabe@rank <- rank
    ausgabe@diff <- diff
    ausgabe@rela <- rela
    ausgabe@fixed <- bootrun@fixed
    #provide names for identifying statistics (ungrouped case)
    names <- bootrun@namen
    if (!is.null(always)) names <- setdiff(names, bootrun@alwaysnam)
    if (!is.null(groupnames)) names <- c(names[1],as.character(bootrun@groupdocu[[1]]))
    diffnam <- paste(names[2:(g+1)][nchoosek(g, 2)[1, ]], names[2:(g+1)][nchoosek(g,2)[2, ]], sep = "-")
#    ausgabe@var.y.boot <- bootrun@boot$t[, 1]
#    ausgabe@R2.boot <- bootrun@boot$t[, 2]
#    ausgabe@R2.decomp.boot <- bootrun@boot$t[, 3]

    if (bootlist.out) ausgabe@bootlist <- bootlist
        else ausgabe@bootlist <- NULL

    #define names for elements from estimate and vcov in order to be able to refer to them later
    zaehl <- 4
    bootnames <- c("var.y", "R2", "R2.decomp")
    typname <- ""
    for (a in c("lmg", "pmvd", "last", "first", "betasq", "pratt")) {
        if (a %in% type) {
            bootnames <- c(bootnames, paste(names[2:(g+1)], ".", a, sep = ""))
            zaehl <- zaehl + g
            if (diff) {
                bootnames <- c(bootnames, paste(diffnam, ".",
                  a, "diff", sep = ""))
                zaehl <- zaehl + g * (g - 1)/2
            }
            if (!typname[1] == "") typname <- c(typname, a)
            else typname <- a
        }
    }
    ## omit always in case it was appended as numeric variable
    ## don't know, whether really necessary
    for (i in 1:m) {
        bootlist[[i]]@boot$t0 <- bootlist[[i]]@boot$t0[1:length(bootnames)] 
        if (!bootlist.out)
        bootlist[[i]]@boot$t <- bootlist[[i]]@boot$t[1:length(bootnames)] 
        else 
        bootlist[[i]]@boot$t <- bootlist[[i]]@boot$t[,1:length(bootnames)] 
        bootlist[[i]]@vcov <- bootlist[[i]]@vcov[1:length(bootnames),1:length(bootnames)] 
    }
    mic <- MIcombine(lapply(bootlist,function(obj) obj@boot$t0), lapply(bootlist,function(obj) obj@vcov))
    names(mic$coefficients) <- bootnames
    names(mic$df) <- bootnames
    names(mic$missinfo) <- bootnames
    colnames(mic$variance) <- bootnames
    rownames(mic$variance) <- bootnames

    #columns of matrices can be referred to by their colnames 
    #(in quotes in square brackets instead of index)
    #elements of vectors analogously
    zaehl <- 4
    for (a in c("lmg", "pmvd", "last", "first", "betasq", "pratt")) {
        if (a %in% type) {
                slot(ausgabe, a) <- mic$coefficients[zaehl:(zaehl + g - 1)]
                zaehl <- zaehl + g
            if (diff) {
                  slot(ausgabe, paste(a, "diff", sep = ".")) <- mic$coefficients[zaehl:(zaehl + g*(g - 1)/2-1)]
                  zaehl <- zaehl + g * (g - 1)/2
            }
        }
    }
    ntype <- length(typname)
    tval <- outer(1-(1-level)/2, mic$df, FUN="qt")
    colnames(tval) <- names(mic$coefficients)
    percentages <- mic$coefficients[paste(names[2:(g+1)], ".", matrix(typname,
        g, ntype, byrow = T), sep = "")]
    if (diff)
        diffs <- mic$coefficients[paste(diffnam, ".", matrix(typname,
            g * (g - 1)/2, ntype, byrow = T), "diff", sep = "")]
    percupper <- matrix(rep(percentages, nlev),nlev,length(percentages),byrow=TRUE,dimnames = list(level,
        names(percentages)))
    perclower <- percupper - tval[,paste(names[2:(g+1)], ".", matrix(typname,
        g, ntype, byrow = T), sep = "")]*
           matrix(rep(sqrt(diag(mic$variance[paste(names[2:(g+1)], ".", matrix(typname,
           g, ntype, byrow = T), sep = ""),paste(names[2:(g+1)], ".", matrix(typname,
           g, ntype, byrow = T), sep = "")])),nlev),nlev,length(percentages),byrow=TRUE)
        dimnames(perclower) <- list(level, names(percentages))
    percupper <- percupper + tval[,paste(names[2:(g+1)], ".", matrix(typname,
        g, ntype, byrow = T), sep = "")]*
        matrix(rep(sqrt(diag(mic$variance[paste(names[2:(g+1)], ".", matrix(typname,
           g, ntype, byrow = T), sep = ""),paste(names[2:(g+1)], ".", matrix(typname,
           g, ntype, byrow = T), sep = "")])),nlev),nlev,length(percentages),byrow=TRUE)
        dimnames(percupper) <- list(level, names(percentages))
    if (diff){
        diffupper <- matrix(rep(diffs, nlev), nlev,length(diffs),byrow=TRUE,dimnames = list(level,names(diff)))
        difflower <- diffupper - tval[,paste(diffnam, ".", matrix(typname,
            g * (g - 1)/2, ntype, byrow = T), "diff", sep = "")]*
            matrix(rep(sqrt(diag(mic$variance)[paste(diffnam, ".", matrix(typname,
            g * (g - 1)/2, ntype, byrow = T), "diff", sep = "")]),nlev),nlev,length(diffs),byrow=TRUE)
        dimnames(difflower) <- list(level, names(diffs))
        diffupper <- diffupper + tval[,paste(diffnam, ".", matrix(typname,
            g * (g - 1)/2, ntype, byrow = T), "diff", sep = "")]*
        matrix(rep(sqrt(diag(mic$variance)[paste(diffnam, ".", matrix(typname,
            g * (g - 1)/2, ntype, byrow = T), "diff", sep = "")]),nlev),nlev,length(diffs),byrow=TRUE)
        dimnames(diffupper) <- list(level, names(diffs))
            }

    ausgabe@type <- typname

    #write confidence bounds to output object
    for (a in typname) {
        slot(ausgabe, paste(a, "lower", sep = ".")) <- matrix(perclower[,
            paste(names[2:(g+1)], a, sep = ".")], nlev, g)
        if (diff)
            slot(ausgabe, paste(a, "diff", "lower", sep = ".")) <- matrix(difflower[,
                paste(diffnam, paste(a, "diff", sep = ""), sep = ".")],
                nlev, g * (g - 1)/2)
        slot(ausgabe, paste(a, "upper", sep = ".")) <- matrix(percupper[,
            paste(names[2:(g+1)], a, sep = ".")], nlev, g)
        if (diff)
            slot(ausgabe, paste(a, "diff", "upper", sep = ".")) <- matrix(diffupper[,
                paste(diffnam, paste(a, "diff", sep = ""), sep = ".")],
                nlev, g * (g - 1)/2)
    }
       #initialize matrix for showing (sorted) results with confidence info
    mark <- matrix(rep(0, (g * ntype + ntype - 1) * (2 *
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

            mark[((aa - 1) * (g + 1) + 1):(aa * (g + 1) - 1),
                ] <- matrix(cbind(percent, t(cilower), t(ciupper)),
                g, 2 * nlev + 1)
        
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
    if (sort) {
        marksort <- format(round(marksort, digits = 4), nsmall=4, scientific=FALSE)
        ausgabe@mark <- marksort 
        }
    else {
        mark <- format(round(mark, digits = 4), nsmall=4, scientific=FALSE)
        ausgabe@mark <- mark
        }


    # differences
    if (diff) {
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
                    hilf[k, i] <- "*"
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
            ausgabe@markdiff <- mark
            }
        if (sort) {
            marksort[, c(1, (2 + nlev):(3 * nlev + 1))] <- format(round(as.numeric(marksort[,
                c(1, (2 + nlev):(3 * nlev + 1))]), digits=4), nsmall=4, scientific=FALSE)
            ausgabe@markdiff <- marksort
        }
    } # if diff 

    # set correct options for printing in output object
    ausgabe@nobs <- bootrun@nobs
    ausgabe@diff <- diff
    ausgabe@rank <- rank
    ausgabe@sort <- sort
    ausgabe@est <- mic$coefficients
    ausgabe@vcov <- mic$variance
    ausgabe@MIresult <- mic
    return(ausgabe)
}

