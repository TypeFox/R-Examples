modgam =function (rdata, rgrid, family = binomial, permute = 0, conditional = TRUE, m = "adjusted", sp = NULL, keep = FALSE, verbose = TRUE, ...)
{
    if (is.character(family)) 
        family = get(family, mode = "function", envir = parent.frame())
    if (is.function(family)) 
        family = family()
    if (is.null(family$family)) {
        print(family)
        stop("`family' not recognized")
    }
    if (!is.null(sp) && !conditional) {
	warning("User-specified span size ignored due to setting conditional = FALSE") 
	sp = NULL
    }
    output = "fit"
    x = length(rdata)
    if (x < 3) 
        stop("rdata must include at least 3 columns")
    if (tolower(m) == "adjusted") {
        m = "Adjusted"
        if (x == 3) 
            stop("cannot run adjusted model without additional covariates in rdata")
        addvars = names(rdata)[!names(rdata) %in% names(rgrid)][-1]
        if (length(addvars) > 0) {
            for (i in addvars) {
                if (is.factor(rdata[, i])) 
                  rgrid[i] = factor(names(which.max(table(rdata[, 
                    i]))), levels = levels(rdata[, i]))
                else rgrid[i] = quantile(rdata[, i], 0.5)
                cat(paste("GAM predictions will use ", i, " = ", 
                  rgrid[1, i], " at all grid points.", sep = ""), 
                  fill = T)
            }
        }
    }     else if (tolower(m) == "unadjusted" | tolower(m) == "crude") {
        x = 3
        m = "Unadjusted"
    }     else stop(paste("model type", m, "not recognized"))
    if (is.null(sp)) 
        sp = optspan(rdata, m, family = family, verbose=F, ...)
    if (x == 3) {
        fmla = as.formula(paste(names(rdata)[1], paste("lo(", 
            paste(names(rdata)[2:3], collapse = ","), ",span=", 
            sp, ")"), sep = "~"))
        fmla.0 = as.formula(paste(names(rdata)[1], paste("1"), 
            sep = "~"))
    }     else {
        fmla = as.formula(paste(names(rdata)[1], paste("lo(", 
            paste(names(rdata)[2:3], collapse = ","), ",span=", 
            sp, ")+", paste(names(rdata)[-(1:3)], collapse = "+")), 
            sep = "~"))
        fmla.0 = as.formula(paste(names(rdata)[1], paste(names(rdata)[-(1:3)], 
            collapse = "+"), sep = "~"))
    }
    if (verbose) {
        cat(paste("The ", tolower(m), " model is: ", sep = ""), 
            fill = T)
        print(fmla, showEnv = F)
        cat(paste(c("Family:", "Link:"), family[1:2]), fill = T)
    }
    model = gam(fmla, data = rdata, family = family, ...)
    model.0 = gam(fmla.0, data = rdata, family = family, ...)
    nullmod = ifelse(x == 3, mean(predict.gam(model.0)), mean(predict.gam(model.0, 
        rgrid)))
    origresults = predict.gam(model, rgrid)
    if (family[1] == "binomial" && family[2] == "logit") 
        output = "OR"
    results = list(grid = rgrid[, 1:2], m = m, span = sp, gamobj = model, 
        family = family, fit = as.vector(origresults))
    if (output == "OR") 
        results$OR = as.vector(exp(origresults - nullmod))
    if (permute > 0) {
        n = length(rgrid[, 1])
        nobs = length(rdata[, 1])
        if (keep) 
            permresults = matrix(NA, n, permute - 1)
        ptranks = rep(1, n)
        devstat = rep(NA, permute)
        ucspans = rep(NA, permute)		# vector of optimal spans
        ucspans[1] = sp
        devstat[1] = anova(model.0, model)$Deviance[2]
        coords = rdata[, 2:3]
        m.data = rdata
        for (i in 2:permute) {
            index = sample(1:nobs, replace = F)
            m.data[, 2:3] = coords[index, ]
            if (!conditional) {
                ucsp = optspan(m.data, m, family = family, verbose=F, ...)
                ucspans[i] = ucsp
                if (x == 3) {
                    fmla = as.formula(paste(names(rdata)[1], paste("lo(",
                    paste(names(rdata)[2:3], collapse = ","), ",span=",
                    ucsp, ")"), sep = "~"))
                    }
                else {
                    fmla = as.formula(paste(names(rdata)[1], paste("lo(",
                    paste(names(rdata)[2:3], collapse = ","), ",span=",
                    ucsp, ")+", paste(names(rdata)[-(1:3)], collapse = "+")),
                    sep = "~"))
                }
            }
            m.gam = gam(fmla, data = m.data, family = family, ...)
            devstat[i] = anova(model.0, m.gam)$Deviance[2]
            tempresults = predict.gam(m.gam, rgrid)
            ptranks = ptranks + (origresults > tempresults)
            if (keep)
                permresults[, i - 1] = tempresults
            if (verbose && i%%10 == 0)
                cat(paste("Permutation", i, "of", permute), fill = TRUE)
        }
        devglobp = (permute - rank(devstat)[1])/permute
        globprint = if (devglobp == 0) 
            paste("<", round(1/permute, 3), sep = "")
        else devglobp
        cat(paste("The global statistic for the ", tolower(m), 
            " model is ", globprint, sep = ""), fill = TRUE)
        results$global = devglobp
        results$deviance = devstat	# return the vector of model deviances
        results$pointwise = as.vector(ptranks/permute)
        if (keep) {
            if (output == "OR") 
                results$permutations = exp(permresults - nullmod)
            else results$permutations = permresults
            results$globaldevs = devstat
            if (!conditional) results$span = ucspans
        }
    }
    return(results)
}
