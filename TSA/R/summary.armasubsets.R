`summary.armasubsets` <-
function (object, all.best = TRUE, matrix = TRUE, matrix.logical = FALSE, 
    df = NULL, ...) 
{
# modified by KSC 
# date: 3/9/11
#
# add the calculation of AIC and AICc
    ll <- object
    triangle <- function(k) {
        j <- k - 1
        1 + j * (j + 1)/2
    }
    nmodl <- ll$nbest * ll$nvmax
    if (all.best) 
        nshow <- ll$nbest
    else nshow <- 1
    if (!is.null(df)) 
        n1 <- df
    else n1 <- ll$nn - ll$intercept
    outmat <- NULL
    rmat <- NULL
    rnames <- NULL
    outnames <- NULL
    rsqvec <- NULL
    cpvec <- NULL
    adjr2vec <- NULL
    bicvec <- NULL
    aicvec <- NULL
    aiccvec <- NULL
    rssvec <- NULL
    sigma2 <- ll$sserr/(n1 + ll$intercept - ll$last)
    for (i in ll$first:min(ll$last, ll$nvmax)) {
        if (!matrix) 
            outmat <- NULL
        for (j in 1:nshow) {
            if (ll$ress[i, j] >= 1e+35) 
                next
            if ((j > 1) & (all(ll$lopt[triangle(i):(triangle(i + 
                1) - 1), j - 1] == ll$lopt[triangle(i):(triangle(i + 
                1) - 1), j]))) 
                next
            rline <- rep(FALSE, ll$np)
            rline[ll$lopt[triangle(i):(triangle(i + 1) - 1), 
                j]] <- TRUE
            outnames <- c(outnames, paste(i - ll$intercept, " (", 
                j, ")"))
            rnames <- c(rnames, as.character(i - ll$intercept))
            rmat <- rbind(rmat, rline)
            vr <- ll$ress[i, j]/ll$nullrss
            rssvec <- c(rssvec, ll$ress[i, j])
            rsqvec <- c(rsqvec, 1 - vr)
            adjr2vec <- c(adjr2vec, 1 - vr * n1/(n1 + ll$intercept - 
                i))
            cpvec <- c(cpvec, ll$ress[i, j]/sigma2 - (n1 + ll$intercept - 
                2 * i))
            bicvec <- c(bicvec, (n1 + ll$intercept) * log(vr) + 
                i * log(n1 + ll$intercept))
            aicvec <- c(aicvec, (n1 + ll$intercept) * log(vr) + 
                i *2)
            aiccvec <-c(aiccvec,  (n1 + ll$intercept) * log(vr) + 
                i *2+2*(i+1)*(i+2)/(n1+ll$intercept-i-2))
        }
    }
    rownames(rmat) <- rnames
    cn <- ll$xnames
    colnames(rmat) <- cn
    reorder <- if (is.null(ll$reorder)) 
        1:NCOL(rmat)
    else c(1, 1 + ll$reorder)
    rmat <- rmat[, order(reorder), drop = FALSE]
    if (matrix) {
        if (!matrix.logical) 
            outmat <- ifelse(rmat, "*", " ")
        else outmat <- rmat
        rownames(outmat) <- outnames
        if (ll$intercept) 
            outmat <- outmat[, -1, drop = FALSE]
    }
    rval <- list(which = rmat, rsq = rsqvec, rss = rssvec, adjR2 = adjr2vec, 
        Cp = cpvec, BIC = bicvec,AIC=aicvec,AICc=aiccvec, outmat = outmat, obj = ll)
    class(rval) <- "summary.regsubsets"
    rval
}

