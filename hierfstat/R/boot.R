boot.vc<-function (levels = levels, loci = loci, diploid = TRUE, nboot = 1000, 
    quant = c(0.025, 0.5, 0.975)) 
{
    gf <- function(dat, num, den) {
        sum(dat[num])/sum(dat[den])
    }
    nloc <- dim(loci)[2]
    if (nloc < 5) {
        stop("Not enough loci to bootstrap. Exiting")
    }
    x <- varcomp.glob(levels = levels, loci = loci, diploid = diploid)
    x.loc <- data.frame(x$loc)
    rows<-complete.cases(x.loc)
    if(sum(rows)<5){
	stop("Not enough polymorphic loci to bootstrap. Exiting")
    }
    nloc<-sum(rows)
    x.loc<-x.loc[rows,]
    nlev <- dim(x.loc)[2]
    names(x.loc) <- names(x$overall)
    mat.boot <- data.frame(matrix(rep(0, nboot * nlev), ncol = nlev))
    for (i in 1:nboot) {
        mat.boot[i, ] <- apply(x.loc[sample(nloc, replace = TRUE), 
            ], 2, sum)
    }
    nam <- vector(length = nlev + 1)
    nam[-1] <- names(x.loc)
    nam[1] <- "Total"
    names(mat.boot) <- names(x.loc)
    mat.res <- data.frame(matrix(rep(0, nboot * (nlev * (nlev + 
        1)/2)), nrow = nboot))
    names.res <- vector(length = nlev * (nlev + 1)/2)
    mat.res[, 1] <- apply(mat.boot, 1, sum)
    acc <- 0
    for (i in 1:(nlev - 1)) {
        acc <- acc + 1
        mat.res[, acc] <- apply(mat.boot[, c(i:nlev)], 1, sum)/nloc
        names.res[acc] <- paste("H-", nam[i], sep = "")
        for (j in i:(nlev - 1)) {
            acc <- acc + 1
            mat.res[, acc] <- apply(mat.boot, 1, gf, num = i:j, 
                den = i:nlev)
            names.res[acc] <- paste("F-", nam[j + 1], "/", nam[i], 
                sep = "")
        }
    }
    acc <- acc + 1
    mat.res[, acc] <- mat.boot[, nlev]/nloc
    names.res[acc] <- "Hobs"
    names(mat.res) <- names.res
    return(list(boot = round(mat.boot,digits=4), res = round(mat.res,digits=4), ci = round(apply(mat.res, 2, quantile, quant),digits=4)))
}

