`make.par.file` <-
function(covariates, form, par.file="parameters.RData", root.dir=".", ...){
    default.vals <- list(
        add.norm = TRUE,
        add.par = 0,
        align.fcn = NA,
        align.method = "spline",
        base.dir = paste(root.dir, "/Baselines", sep=""),
        bhbysubj = FALSE,
        calc.all.peaks = FALSE,
        cluster.constant = 10,
        cluster.method = "ppm",
        cor.thresh = 0.8,
        FDR = 0.1,
        FTICRMS.version = "0.8",
        gengamma.quantiles = TRUE,
        halve.search = FALSE,
        isotope.dist = 7,
        lrg.dir = paste(root.dir, "/Large_Peaks", sep=""),
        lrg.file = "lrg_peaks.RData",
        lrg.only = TRUE,
        masses = NA,
        max.iter = 20,
        min.spect = 1,
        neg.div = NA,
        neg.norm.by = "baseline",
        norm.peaks = "common",
        norm.post.repl = FALSE,
        num.pts = 5,
        oneside.min = 1,
        overwrite = FALSE,
        par.file = "parameters.RData",
        peak.dir = paste(root.dir, "/All_Peaks", sep=""),
        peak.method = "parabola",
        peak.thresh = 3.798194,
        pre.align = FALSE,
        pval.fcn = "default",
        R2.thresh = 0.98,
        raw.dir = paste(root.dir, "/Raw_Data", sep=""),
        rel.conv.crit = TRUE,
        repl.method = "max",
        res.dir = paste(root.dir, "/Results", sep=""),
        res.file = "analyzed.RData",
        root.dir = ".",
        sm.div = NA,
        sm.norm.by = "baseline",
        sm.ord = 2,
        sm.par = 1e-11,
        subtract.base = FALSE,
        tol = 5e-8,
        trans.method = "shiftedlog",
        use.model = "lm",
        zero.rm = TRUE
    )
    new.vals <- c(list(...), root.dir=root.dir, par.file=par.file)
    for(i in names(default.vals)){
        if(i %in% names(new.vals)){
            assign(i, new.vals[[i]])
        } else {
            assign(i, default.vals[[i]])        
        }
    }
    tmp <- !(names(new.vals) %in% names(default.vals))
    if(sum(tmp)==1){
        warning(paste(names(new.vals)[tmp], "is not a valid parameter name"))
    } else if(sum(tmp)==2){
        warning(paste(paste(names(new.vals)[tmp], collapse=" and "), 
            "are not valid parameter names"))
    } else if(sum(tmp)>2){
        names(new.vals)[tmp][sum(tmp)] <- paste("and", names(new.vals)[tmp][sum(tmp)])
        warning(paste(paste(names(new.vals)[tmp], collapse=", "), 
            "are not valid parameter names"))
    }
    save(list=setdiff(ls(), c("default.vals", "new.vals", "tmp", "i")), 
        file=paste(root.dir, "/", par.file, sep=""))
}

