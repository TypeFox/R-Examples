## tests
#devtools::install_github("datacloning/dcmle")

library(dcmle)

## data type classes

as(new("gsFit"), "dcFit")
as(new("dcFit"), "gsFit")

## dcmle model fit classes

str(as.mcmc.list(regmod))
str(as(regmod, "mcmc.list"))
str(as(regmod, "MCMClist"))
str(as(regmod, "codaMCMC"))
str(as(regmod, "dcCodaMCMC"))
str(as(regmod, "dcmle"))

str(as(as(regmod, "MCMClist"), "MCMClist"))
str(as(as(regmod, "MCMClist"), "codaMCMC"))
str(as(as(regmod, "MCMClist"), "dcCodaMCMC"))
str(as(as(regmod, "MCMClist"), "dcmle"))

str(as(as(regmod, "codaMCMC"), "MCMClist"))
str(as(as(regmod, "codaMCMC"), "codaMCMC"))
str(as(as(regmod, "codaMCMC"), "dcCodaMCMC"))
str(as(as(regmod, "codaMCMC"), "dcmle"))

str(as(as(regmod, "dcCodaMCMC"), "MCMClist"))
str(as(as(regmod, "dcCodaMCMC"), "codaMCMC"))
str(as(as(regmod, "dcCodaMCMC"), "dcCodaMCMC"))
str(as(as(regmod, "dcCodaMCMC"), "dcmle"))

str(as(as(regmod, "dcmle"), "MCMClist"))
str(as(as(regmod, "dcmle"), "codaMCMC"))
str(as(as(regmod, "dcmle"), "dcCodaMCMC"))
str(as(as(regmod, "dcmle"), "dcmle"))

## testing plot methods

regmod_MCMClist <- as(regmod, "MCMClist")
regmod_codaMCMC <- as(regmod, "codaMCMC")
regmod_dcCodaMCMC <- as(regmod, "dcCodaMCMC")
regmod_dcmle <- as(regmod, "dcmle")

evalfun <- function(FUN) {
    evalfun_int <- function(x, FUN) {
        eval(parse(text=
          paste(FUN, "(", x, ")", sep="")
          ))
    }
    out <- rep(0L, 4)
    names(out) <- c("MCMClist", "codaMCMC", "dcCodaMCMC", "dcmle")

    cat("\n\n***", FUN, "***")

    cat("\n\n---", FUN, "--- MCMClist ---\n")
    ooo <- try(evalfun_int("regmod_MCMClist", FUN), silent=TRUE)
    if (inherits(ooo, "try-error")) {
        cat(as.character(ooo), "\n\n")
        out[1] <- 1L
    } else {
        cat("OK\n\n")
    }

    cat("\n\n---", FUN, "--- codaMCMC ---\n")
    ooo <- try(evalfun_int("regmod_codaMCMC", FUN), silent=TRUE)
    if (inherits(ooo, "try-error")) {
        cat(as.character(ooo), "\n\n")
        out[2] <- 1L
    } else {
        cat("OK\n\n")
    }

    cat("\n\n---", FUN, "--- dcCodaMCMC ---\n")
    ooo <- try(evalfun_int("regmod_dcCodaMCMC", FUN), silent=TRUE)
    if (inherits(ooo, "try-error")) {
        cat(as.character(ooo), "\n\n")
        out[3] <- 1L
    } else {
        cat("OK\n\n")
    }

    cat("\n\n---", FUN, "--- dcmle ---\n")
    ooo <- try(evalfun_int("regmod_dcmle", FUN), silent=TRUE)
    if (inherits(ooo, "try-error")) {
        cat(as.character(ooo), "\n\n")
        out[4] <- 1L
    } else {
        cat("OK\n\n")
    }
    out
}
#evalfun("str")

toEval <- c("plot",
  "traceplot",
  "densplot",
  "pairs",
  "densityplot",
  "qqmath",
  "xyplot",
  "acfplot",
  "crosscorr.plot",
  "dcdiag",
  "dctable",
  "nclones",
  "dcsd",
  "as.matrix",
  "as.array",
  "nvar",
  "varnames",
  "chanames",
  "nchain",
  "niter",
  "crosscorr",
  "mcpar",
  "thin",
  "coef",
  "vcov",
#  "confint",
  "quantile",
  "start",
  "end",
  "frequency",
  "time",
  "window",
  "stack",
  "str",
  "head",
  "tail",
  "autocorr.diag",
  "lambdamax.diag",
  "chisq.diag",
  "gelman.diag",
  "geweke.diag",
  "raftery.diag",
  "heidel.diag")

res <- list()
for (i in toEval) {
  res[[i]] <- evalfun(i)
}

evalfun("confint")
attr(regmod, "n.clones") <- 2
regmod_MCMClist <- as(regmod, "MCMClist")
regmod_codaMCMC <- as(regmod, "codaMCMC")
regmod_dcCodaMCMC <- as(regmod, "dcCodaMCMC")
regmod_dcmle <- as(regmod, "dcmle")
evalfun("confint")

## testing stats/coda methods

#tmp <- mcmc(cbind(z=c(1,1,1,1), a=c(0,0,0,0), q2=c(-1,-1,-1,-1)),
#    start=101, end=107, thin=2)
#x <- as.mcmc.list(list(tmp, tmp+0.5, tmp-0.5))
#rm(tmp)

## update
## show, summary
## [, [[

tmp <- do.call(rbind, res)

## --------- TOTAL -----------
colSums(tmp)
tmp[rowSums(tmp) != 0,]

