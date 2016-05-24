## test that should make sure that the results remain constant
require(robustlmm)

fit <- function(formula, data, methods =  c("DASvar", "DAStau"),
                rho.e = cPsi, rho.b = cPsi, ...) {
    fits <- list()
    ## compare with result of lmer if rho arguments are not given
    classic <- ! any(c("rho.e", "rho.b") %in% names(match.call())[-1])
    if (classic) fm <- lmer(formula, data)
    for (method in methods) {
        fits[[method]] <- list()
        if (classic) fits[[method]][["lmer"]] <- fm
        cat("\n########", method, "########\n")
        try({cat("Time elapsed:",
                 system.time(m <- rlmer(formula, data, method=method,
                                        rho.e = rho.e, rho.b = rho.b, ...)),
                 "\n")
             fits[[method]][["IRWLS"]] <- m
             print(summary(m))
             print(robustlmm:::u.rlmerMod(m), 4)
             if (classic) {
                 ## compare with lmer fit
                 cat("#### Checking equality with lmer... ####\n")
                 cat("Fixed effects: ", all.equal(fixef(fm), fixef(m), tolerance = 1e-4), "\n")
                 cat("Random effects:", all.equal(ranef(fm), ranef(m), tolerance = 1e-4,
                                                  check.attributes=FALSE), "\n")
                 cat("Theta:         ", all.equal(theta(fm), theta(m), tolerance = 1e-4), "\n")
                 cat("Sigma:         ", all.equal(sigma(fm), sigma(m), tolerance = 1e-4), "\n")
                 if (packageVersion("lme4") >= "0.99999911.0") {
                     tmp <- all.equal(fm@pp$unsc(), unname(m@pp$unsc()), tolerance = 1e-4)
                     if (!isTRUE(tmp))
                         cat("Unsc:          ", tmp , "\n")
                 }
             }
         })
        fits[[method]][["dnames"]] <- names(fits[[method]])
    }
    cat("\n################################################\n")
    cat("################################################\n")
    cat("################################################\n")
    for (method in methods) {
        cat("\n################ results for",method," ##############\n")
        cmp <- do.call(compare, fits[[method]])
        cmp <- cmp[grep("^rho", rownames(cmp), invert=TRUE),,drop=FALSE]
        print.default(cmp, quote="FALSE")
    }
}

if (FALSE) {
    ## classic (REML)
    fit(Yield ~ (1 | Batch), Dyestuff)
    fit(Yield ~ (1 | Batch), Dyestuff2)
    fit(diameter ~ (1|plate) + (1|sample), Penicillin)

    ## classic (no init)
    fit(Yield ~ (1 | Batch), Dyestuff, init = lmerNoFit)
    fit(Yield ~ (1 | Batch), Dyestuff2, init = lmerNoFit)
    fit(diameter ~ (1|plate) + (1|sample), Penicillin, init = lmerNoFit)

    ## smoothPsi, wExp = 1
    fit(Yield ~ (1 | Batch), Dyestuff,
        rho.e = smoothPsi, rho.b = smoothPsi,
        rho.sigma.b = smoothPsi, rho.sigma.e = smoothPsi)
    fit(Yield ~ (1 | Batch), Dyestuff2,
        rho.e = smoothPsi, rho.b = smoothPsi,
        rho.sigma.b = smoothPsi, rho.sigma.e = smoothPsi)
    fit(diameter ~ (1|plate) + (1|sample), Penicillin,
        rho.e = smoothPsi, rho.b = smoothPsi,
        rho.sigma.b = smoothPsi, rho.sigma.e = smoothPsi)

    ## smoothPsi Proposal II for estimating sigma only
    fit(Yield ~ (1 | Batch), Dyestuff,
        rho.e = smoothPsi, rho.b = smoothPsi,
        rho.sigma.b = smoothPsi)
    fit(Yield ~ (1 | Batch), Dyestuff2,
        rho.e = smoothPsi, rho.b = smoothPsi,
        rho.sigma.b = smoothPsi)
    fit(diameter ~ (1|plate) + (1|sample), Penicillin,
        rho.e = smoothPsi, rho.b = smoothPsi,
        rho.sigma.b = smoothPsi)
}

## smoothPsi Proposal II (default)
fit(Yield ~ (1 | Batch), Dyestuff,
    rho.e = smoothPsi, rho.b = smoothPsi)
fit(Yield ~ (1 | Batch), Dyestuff2,
    rho.e = smoothPsi, rho.b = smoothPsi)
fit(diameter ~ (1|plate) + (1|sample), Penicillin,
    rho.e = smoothPsi, rho.b = smoothPsi)

if (FALSE) {
    ## correlated random effects
    fit(Reaction ~ Days + (Days|Subject), sleepstudy,
        methods = c("DASvar", "DAStau"))
    fit(Reaction ~ Days + (Days|Subject), sleepstudy,
        methods = c("DASvar", "DAStau"), init = lmerNoFit)
    ## robust
    fit(Reaction ~ Days + (Days|Subject), sleepstudy,
        rho.e = smoothPsi, rho.b = smoothPsi,
        methods = c("DASvar")) ##, "DAStau"))
    fit(Reaction ~ Days + (Days|Subject), sleepstudy,
        rho.e = smoothPsi, rho.b = smoothPsi,
        methods = c("DASvar"), ##, "DAStau"),
        init = lmerNoFit)

    ## ## including a 0 variance compontent
    sleepstudy2 <- within(sleepstudy, Group <- letters[1:4])
    fit(Reaction ~ Days + (Days|Subject) + (1|Group), sleepstudy2,
        methods = c("DASvar", "DAStau"))
    fit(Reaction ~ Days + (Days|Subject) + (1|Group), sleepstudy2,
        methods = c("DASvar", "DAStau"), init = lmerNoFit)
    ## robust
    fit(Reaction ~ Days + (Days|Subject) + (1|Group), sleepstudy2,
        rho.e = smoothPsi, rho.b = smoothPsi,
        methods = c("DASvar")) ##, "DAStau"))
    fit(Reaction ~ Days + (Days|Subject) + (1|Group), sleepstudy2,
        rho.e = smoothPsi, rho.b = smoothPsi,
        methods = c("DASvar"), ##, "DAStau"),
        init = lmerNoFit)
}

cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''
