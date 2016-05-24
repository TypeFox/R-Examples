escalc.formula <-
function (measure, formula, weights, data, add = 1/2, to = "only0", 
    drop00 = FALSE, vtype = "LS", var.names = c("yi", "vi"), 
    digits = 4, ...) 
{
    if (!is.element(measure, c("RR", "OR", "PETO", "RD", "AS", 
        "PHI", "YUQ", "YUY", "RTET", "PBIT", "OR2D", "OR2DN", 
        "OR2DL", "IRR", "IRD", "IRSD", "MD", "SMD", "SMDH", "ROM", 
        "RPB", "RBIS", "D2OR", "D2ORN", "D2ORL", "COR", "UCOR", 
        "ZCOR", "PR", "PLN", "PLO", "PAS", "PFT", "IR", "IRLN", 
        "IRS", "IRFT", "MN", "MC", "SMCC", "SMCR", "SMCRH", "ROMC", 
        "ARAW", "AHW", "ABT"))) 
        stop("Unknown 'measure' specified.")
    if (is.element(measure, c("MC", "SMCC", "SMCR", "SMCRH", 
        "ROMC"))) {
        stop("Formula interface (currently) not implemented for these outcome measures.")
    }
    if (!requireNamespace("Formula", quietly = TRUE)) 
        stop("Please install the 'Formula' package to use the formula interface.")
    na.act <- getOption("na.action")
    options(na.action = "na.pass")
    on.exit(options(na.action = na.act))
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "weights"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    formula <- Formula::as.Formula(formula)
    if (length(formula)[2] < 2L) 
        stop("Right-hand side of formula must specify both a grouping and a study factor (i.e., ~ group | study).")
    mf$formula <- formula
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    weights <- model.weights(mf)
    lhs <- Formula::model.part(formula, data = mf, lhs = 1)
    rhs1 <- Formula::model.part(formula, data = mf, rhs = 1)
    study <- Formula::model.part(formula, data = mf, rhs = 2)
    if (length(study) != 1) 
        stop("A single study factor must be specified.")
    if (!is.factor(study[[1]])) 
        stop("Study variable must be a factor.")
    study <- study[[1]]
    if (anyNA(study)) 
        stop("Study factor must not contain NAs.")
    if (is.element(measure, c("RR", "OR", "RD", "AS", "PETO", 
        "PHI", "YUQ", "YUY", "RTET", "PBIT", "OR2D", "OR2DN", 
        "OR2DL"))) {
        if (is.null(weights)) 
            stop("Must specify the 'weights' argument.")
        if (length(lhs) != 1) 
            stop("Left-hand side of formula must be a single outcome factor.")
        outcome <- lhs[[1]]
        if (!is.factor(outcome)) 
            stop("Left-hand side of formula must be a factor.")
        if (nlevels(outcome) != 2) 
            stop("Outcome factor on left-hand side of formula should have two levels.")
        if (length(rhs1) != 1) 
            stop("A single grouping factor must be specified.")
        if (!is.factor(rhs1[[1]])) 
            stop("Grouping variable must be a factor.")
        group <- rhs1[[1]]
        if (nlevels(group) != 2) 
            stop("Grouping factor should have two levels.")
        if (anyNA(group) || anyNA(outcome)) 
            stop("Grouping and outcome factors must not contain NAs.")
        ai <- weights[group == levels(group)[1] & outcome == 
            levels(outcome)[1]]
        bi <- weights[group == levels(group)[1] & outcome == 
            levels(outcome)[2]]
        ci <- weights[group == levels(group)[2] & outcome == 
            levels(outcome)[1]]
        di <- weights[group == levels(group)[2] & outcome == 
            levels(outcome)[2]]
        names(ai) <- mf$study[group == levels(group)[1] & outcome == 
            levels(outcome)[1]]
        names(bi) <- mf$study[group == levels(group)[1] & outcome == 
            levels(outcome)[2]]
        names(ci) <- mf$study[group == levels(group)[2] & outcome == 
            levels(outcome)[1]]
        names(di) <- mf$study[group == levels(group)[2] & outcome == 
            levels(outcome)[2]]
        return(escalc(measure = measure, ai = ai, bi = bi, ci = ci, 
            di = di, add = add, to = to, drop00 = drop00, vtype = vtype, 
            var.names = var.names, append = "FALSE", digits = digits))
    }
    if (is.element(measure, c("IRR", "IRD", "IRSD"))) {
        if (length(lhs) != 2) 
            stop("Left-hand side of formula must specify the number of events and the total person-time at risk (i.e., events/times ~).")
        events <- lhs[, 1]
        times <- lhs[, 2]
        if (!is.vector(events) || !is.vector(times)) 
            stop("The events and person-time at risk variables should be vectors.")
        if (length(rhs1) != 1) 
            stop("A single grouping factor must be specified.")
        if (!is.factor(rhs1[[1]])) 
            stop("Grouping variable must be a factor.")
        group <- rhs1[[1]]
        if (nlevels(group) != 2) 
            stop("Grouping factor should have two levels.")
        if (anyNA(group)) 
            stop("Grouping factor must not contain NAs.")
        x1i <- events[group == levels(group)[1]]
        x2i <- events[group == levels(group)[2]]
        t1i <- times[group == levels(group)[1]]
        t2i <- times[group == levels(group)[2]]
        names(x1i) <- mf$study[group == levels(group)[1]]
        names(x2i) <- mf$study[group == levels(group)[2]]
        return(escalc(measure = measure, x1i = x1i, x2i = x2i, 
            t1i = t1i, t2i = t2i, add = add, to = to, drop00 = drop00, 
            vtype = vtype, var.names = var.names, append = "FALSE", 
            digits = digits))
    }
    if (is.element(measure, c("MD", "SMD", "SMDH", "ROM", "RPB", 
        "RBIS", "D2OR", "D2ORN", "D2ORL"))) {
        if (is.null(weights)) 
            stop("Must specify the 'weights' argument.")
        if (length(lhs) != 2) 
            stop("Left-hand side of formula must specify the means and standard devations (i.e., means/sds ~).")
        means <- lhs[, 1]
        sds <- lhs[, 2]
        if (!is.vector(means) || !is.vector(sds)) 
            stop("The mean and standard devation variables should be vectors.")
        if (length(rhs1) != 1) 
            stop("A single grouping factor must be specified.")
        if (!is.factor(rhs1[[1]])) 
            stop("Grouping variable must be a factor.")
        group <- rhs1[[1]]
        if (nlevels(group) != 2) 
            stop("Grouping factor should have two levels.")
        if (anyNA(group)) 
            stop("Grouping factor must not contain NAs.")
        m1i <- means[group == levels(group)[1]]
        m2i <- means[group == levels(group)[2]]
        sd1i <- sds[group == levels(group)[1]]
        sd2i <- sds[group == levels(group)[2]]
        n1i <- weights[group == levels(group)[1]]
        n2i <- weights[group == levels(group)[2]]
        names(m1i) <- mf$study[group == levels(group)[1]]
        names(m2i) <- mf$study[group == levels(group)[2]]
        return(escalc(measure = measure, m1i = m1i, m2i = m2i, 
            sd1i = sd1i, sd2i = sd2i, n1i = n1i, n2i = n2i, vtype = vtype, 
            var.names = var.names, append = "FALSE", digits = digits))
    }
    if (is.element(measure, c("COR", "UCOR", "ZCOR"))) {
        if (is.null(weights)) 
            stop("Must specify the 'weights' argument.")
        if (length(lhs) != 1) 
            stop("Left-hand side of formula must specify the correlations (i.e., cors ~).")
        ri <- lhs[[1]]
        if (!is.vector(ri)) 
            stop("The variable specifying the correlation should be a vector.")
        ni <- weights
        names(ri) <- mf$study
        return(escalc(measure = measure, ri = ri, ni = ni, vtype = vtype, 
            var.names = var.names, append = "FALSE", digits = digits))
    }
    if (is.element(measure, c("PR", "PLN", "PLO", "PAS", "PFT"))) {
        if (length(lhs) != 1) 
            stop("Left-hand side of formula must be a single outcome factor.")
        outcome <- lhs[[1]]
        if (!is.factor(outcome)) 
            stop("Left-hand side of formula must be a factor.")
        if (nlevels(outcome) != 2) 
            stop("Outcome factor on left-hand side of formula should have two levels.")
        if (anyNA(outcome)) 
            stop("Outcome factor must not contain NAs.")
        xi <- weights[outcome == levels(outcome)[1]]
        mi <- weights[outcome == levels(outcome)[2]]
        names(xi) <- mf$study[outcome == levels(outcome)[1]]
        names(mi) <- mf$study[outcome == levels(outcome)[2]]
        return(escalc(measure = measure, xi = xi, mi = mi, add = add, 
            to = to, vtype = vtype, var.names = var.names, append = "FALSE", 
            digits = digits))
    }
    if (is.element(measure, c("IR", "IRLN", "IRS", "IRFT"))) {
        if (length(lhs) != 2) 
            stop("Left-hand side of formula must specify the number of cases and the total person-time at risk (i.e., cases/times ~).")
        events <- lhs[, 1]
        times <- lhs[, 2]
        if (!is.vector(events) || !is.vector(times)) 
            stop("The events and person-time at risk variables should be vectors.")
        xi <- events
        ti <- times
        names(xi) <- mf$study
        return(escalc(measure = measure, xi = xi, ti = ti, add = add, 
            to = to, vtype = vtype, var.names = var.names, append = "FALSE", 
            digits = digits))
    }
    if (is.element(measure, c("MN"))) {
        if (is.null(weights)) 
            stop("Must specify the 'weights' argument.")
        if (length(lhs) != 2) 
            stop("Left-hand side of formula must specify the means and standard devations (i.e., means/sds ~).")
        means <- lhs[, 1]
        sds <- lhs[, 2]
        if (!is.vector(means) || !is.vector(sds)) 
            stop("The mean and standard devation variables should be vectors.")
        mi <- means
        sdi <- sds
        ni <- weights
        names(mi) <- mf$study
        return(escalc(measure = measure, mi = mi, sdi = sdi, 
            ni = ni, vtype = vtype, var.names = var.names, append = "FALSE", 
            digits = digits))
    }
    if (is.element(measure, c("ARAW", "AHW", "ABT"))) {
        if (is.null(weights)) 
            stop("Must specify the 'weights' argument.")
        if (length(lhs) != 2) 
            stop("Left-hand side of formula must specify the alpha values and number of items (i.e., alphas/items ~).")
        alphas <- lhs[, 1]
        items <- lhs[, 2]
        if (!is.vector(alphas) || !is.vector(items)) 
            stop("The alpha and item variables should be vectors.")
        ai <- alphas
        mi <- items
        ni <- weights
        names(ai) <- mf$study
        return(escalc(measure = measure, ai = ai, mi = mi, ni = ni, 
            vtype = vtype, var.names = var.names, append = "FALSE", 
            digits = digits))
    }
}
