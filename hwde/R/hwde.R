hwde <-
function (data = hwde::IndianIrish, gp = "Population", termlist = NULL,
    refmodel = NULL, loci = paste("locus", 1:(dim(data)[2] -
        1), sep = ""), observed = "Observed", keep.models = FALSE,
    aovtable.print = TRUE, group.terms = TRUE, allele.chars = letters)
{
   tmaker <- function(trms = c("ma", "mb"), currbase = 1, group.terms = TRUE,
        pref = "") {
        if (length(trms) == 1) {
            refmodel <- currbase
            plus <- trms
        }
        else if (group.terms) {
            plus <- paste("(", paste(trms, collapse = "+"), ")",
                sep = "")
            refmodel <- currbase
        }
        else {
            plus <- paste(pref, c(trms, paste("(", paste(trms,
                collapse = "+"), ")", sep = "")), sep = "")
            refmodel <- rep(currbase, length(plus))
        }
        list(plus, refmodel)
    }
    nloci <- length(loci)
    colnames <- names(data)
    for (i in 1:length(loci)) if (!(loci[i] %in% colnames)) {
        nloci <- i - 1
        loci <- loci[1:nloci]
        break
    }
    obs <- data[, observed]
    nobs <- length(obs)
    gploc <- loci
    contrasts.info <- make.contrasts(data = data[, loci], allele.chars = allele.chars)
    contr.df <- contrasts.info$contrasts.df
    oset <- contrasts.info$oset
    list.columns <- contrasts.info$list.columns
    if (gp %in% colnames)
        data.df <- cbind(data.frame(obs = obs, gp = data[, gp]),
            contr.df, oset=oset)
    else data.df <- cbind(data.frame(obs = obs), contr.df,
                          oset=oset)
    gpslot <- as.numeric("gp" %in% names(data.df))
    if (is.null(termlist)) {
        addterms <- NULL
        newbase <- 1
        refmodel <- NULL
        for (colvec in list.columns) {
            trms <- tmaker(colvec, currbase = newbase, group.terms = group.terms)
            addterms <- c(addterms, trms[[1]])
            newbase <- newbase + length(trms[[1]])
            refmodel <- c(refmodel, trms[[2]])
        }
        if (gpslot) {
            addtermsg <- paste("gp:", addterms, sep = "")
            addterms <- c("gp", addterms, addtermsg)
            refmodel <- c(1, refmodel + 1, refmodel + newbase)
        }
        addterms <- paste("+", addterms, sep = "")
    }
    else {
        addterms <- termlist
    }
    form1 <- formula(paste("obs", "~", "1"))
    m1 <- glm(form1, family = poisson, offset = log(oset), data = data.df)
    n <- length(addterms)
    firstchar <- rep(" ", n)
    currmod <- m1
    mlist <- list(m1)
    for (i in 1:n) {
        modi <- paste("m", i + 1, sep = "")
        updi <- formula(paste(".~.", addterms[i], sep = ""))
        newmod <- update(currmod, updi)
        if ((i < n) & (refmodel[i + 1] > refmodel[i])) {
            if (!group.terms)
                firstchar[i] <- "r"
            else firstchar <- " "
            currmod <- newmod
        }
        assign(paste("m", i + 1, sep = ""), newmod)
        mlist <- c(mlist, list(newmod))
    }
    anovatab <- do.call("anova", mlist)
    nam <- c("m", paste(firstchar, addterms, sep = ""))
    if (!group.terms)
        nam[1] <- paste("r", nam[1], sep = "")
    anovatab[2:(n + 1), "Deviance"] <- anovatab[refmodel, "Resid. Dev"] -
        anovatab[2:(n + 1), "Resid. Dev"]
    anovatab[2:(n + 1), "Df"] <- anovatab[refmodel, "Resid. Df"] -
        anovatab[2:(n + 1), "Resid. Df"]
    aovtab.terms <- attributes(anovatab)$heading[2]
    row.names(anovatab) <- c(refmodel[1], addterms)
    attributes(anovatab)$heading <- NULL
    if (aovtable.print == TRUE) {
        print("Analysis of Deviance Table")
        print(anovatab, rowlab = nam)
    }
    if (!keep.models)
        invisible(list(anovatab = anovatab, data.df = data.df,
            aovtab.terms = aovtab.terms))
    else invisible(list(anovatab = anovatab, data.df = data.df,
        aovtab.terms = aovtab.terms, models = mlist))
}
