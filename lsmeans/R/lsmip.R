# lsmip code - interaction plots

lsmip = function(object, formula, ...)
    UseMethod("lsmip")

# object - a model object supported by lsmeans
# formula - a formula of the form  x.factors ~ trace.factors | panel.factors
lsmip.default = function(object, formula, type,  
        pch=c(1,2,6,7,9,10,15:20), lty=1, col=NULL, plotit = TRUE, ...) {
    if (!requireNamespace("lattice"))
        stop("This function requires the 'lattice' package be installed.")
    if (length(formula) < 3)
        formula = reformulate(as.character(formula)[[2]], response = ".single.")
        ###stop("'formula' must be two-sided, e.g. trace.factor ~ x.factor")
        ### NEW: Allow lhs to be empty, so then we get a single trace
    
    # Glean the parts of ... to use in lsmeans call
    # arguments allowed to be passed
    lsa.allowed = c("at","trend","cov.reduce","fac.reduce")
    xargs = list(...)
    lsmopts = list(...)
    for (arg in names(xargs)) {
        idx = pmatch(arg, lsa.allowed)
        if (!is.na(idx)) {
            opt = lsa.allowed[idx]
            lsmopts[[opt]] = xargs[[arg]]
            xargs[[arg]] = NULL
        }
    }
    
    allvars = setdiff(.all.vars(formula), ".single.")
    lsmopts$object = object
    lsmopts$specs = reformulate(allvars)
    lsmo = do.call("lsmeans", lsmopts)
    if(missing(type)) {
        type = get.lsm.option("summary")$predict.type
        if (is.null(type))
            type = .get.predict.type(lsmo@misc)
    }
    type = .validate.type(type)

    lsm = predict(lsmo, type = type)
    lsms = cbind(lsmo@grid, lsmean = lsm)

    # Set up trace vars and key
    tvars = .all.vars(update(formula, . ~ 1))
    if (all(tvars == ".single.")) {
        lsms$.single. = 1
        my.key = function(tvars) list()
    }
    else {
        my.key = function(tvars) 
            list(space="right", 
                 title = paste(tvars, collapse=" * "), 
                 points = TRUE, 
                 lines=length(lty) > 1,
                 cex.title=1)
    }
    tv = do.call(paste, lsms[tvars])
    lsms$tvar = factor(tv, levels=unique(tv))
    
    # figure out 'x' and 'by' vars
    rhs = strsplit(as.character(formula[3]), "\\|")[[1]]
    xvars = .all.vars(reformulate(rhs[[1]]))
    xv = do.call(paste, lsms[xvars])
    lsms$xvar = factor(xv, levels = unique(xv))
    lsms = lsms[order(lsms$xvar), ]
    plotform = lsmean ~ xvar
    
    # see if we have any 'by' vars
    if (length(rhs) > 1) {
        byvars = .all.vars(reformulate(rhs[[2]]))
        plotform = as.formula(paste("lsmean ~ xvar |", paste(byvars, collapse="*")))
    }

    # The strips the way I want them
    my.strip = lattice::strip.custom(strip.names = c(TRUE,TRUE), strip.levels = c(TRUE,TRUE), sep = " = ")
    
    TP = TP.orig = lattice::trellis.par.get()
    TP$superpose.symbol$pch = pch
    TP$superpose.line$lty = lty
    if (!is.null(col)) TP$superpose.symbol$col = TP$superpose.line$col = col
    lattice::trellis.par.set(TP)
    
    xlab = ifelse(is.null(xargs$xlab),
        paste("Levels of", paste(xvars, collapse=" * ")), xargs$xlab)
    rspLbl = paste("Predicted", 
        ifelse(is.null(lsmo@misc$inv.lbl), "response", lsmo@misc$inv.lbl))
    ylab = ifelse(is.null(xargs$ylab),
        ifelse(type == "response", rspLbl, "Linear prediction"),
        xargs$ylab)
    
    # remove the unneeded stuff from xlabs
    xargs = xargs[setdiff(names(xargs), c("xlab","ylab"))]
    plotspecs = list(x = plotform, data = lsms, groups = ~ tvar, 
        xlab = xlab, ylab = ylab,
        strip = my.strip, auto.key = my.key(tvars), type=c("p","l"))
    grobj = do.call(lattice::xyplot, c(plotspecs, xargs))
    if (plotit)
        print(grobj)
    attr(lsms, "lattice") = grobj
    lattice::trellis.par.set(TP.orig)
    invisible(lsms)
}
