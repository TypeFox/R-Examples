`ensemble.plot` <- function(
    RASTER.species.name="Species001", RASTER.stack.name="base",
    plot.method="suitability", 
    dev.new.width=7, dev.new.height=7,
    main=paste(RASTER.species.name, " ", plot.method, " for ", RASTER.stack.name, sep=""),
    positive.filters=c("grd","_MEAN_"), negative.filters=c("xml"), 
    p=NULL, a=NULL,
    threshold=-1,
    threshold.method="spec_sens", threshold.sensitivity=0.9, threshold.PresenceAbsence=FALSE,
    abs.breaks=6, pres.breaks=6, ...
)
{
    .BiodiversityR <- new.env()
#    if (! require(dismo)) {stop("Please install the dismo package")}
    if (threshold < 0  && plot.method=="suitability") {
        if (is.null(p)==T || is.null(a)==T) {stop(paste("Please provide locations p and a to calculate thresholds", "\n", sep = ""))}
    }

# new methods for calculating thresholds
    threshold2 <- function(eval, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, 
            threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=Pres, Abs=Abs) {
        if (threshold.PresenceAbsence == T){        
            Pres2 <- cbind(rep(1, length(Pres)), Pres)
            Abs2 <- cbind(rep(0, length(Abs)), Abs)
            data1 <- rbind(Pres2, Abs2)
            data2 <- cbind(seq(1:nrow(data1)), data1)
            auc.value <- PresenceAbsence::auc(data2, st.dev=F)
            cat(paste("\n", "AUC from PresenceAbsence package (also used to calculate threshold): ", auc.value, "\n", sep = ""))
            if (threshold.method=="kappa") {threshold.method <- "MaxKappa"}
            if (threshold.method=="spec_sens") {threshold.method <- "MaxSens+Spec"}
            if (threshold.method=="prevalence") {threshold.method <- "ObsPrev"}
            if (threshold.method=="equal_sens_spec") {threshold.method <- "Sens=Spec"}
            if (threshold.method=="sensitivity") {threshold.method <- "ReqSens"}
            req.sens <- threshold.sensitivity
            if (threshold.method=="no_omission") {
                threshold.method <- "ReqSens"
                req.sens <- 1.0
            }
            result <- PresenceAbsence::optimal.thresholds(data2, threshold=seq(from=0, to=1, by=0.005), req.sens=req.sens)
            result2 <- as.numeric(result[, 2])
            names(result2) <- result[, 1]
            if (threshold.method == "threshold.min") {
                t1 <- result2[["MaxSens+Spec"]]
                t2 <- result2[["Sens=Spec"]]            
                t3 <- result2[["ObsPrev"]]
                thresholds <- as.numeric(c(t1, t2, t3))
                thresholds <- thresholds[thresholds > 0]
                return(min(thresholds))
            }
            if (threshold.method == "threshold.mean") {
                t1 <- result2[["MaxSens+Spec"]]
                t2 <- result2[["Sens=Spec"]]            
                t3 <- result2[["ObsPrev"]]
                thresholds <- as.numeric(c(t1, t2, t3))
                thresholds <- thresholds[thresholds > 0]
                return(mean(thresholds))
            }
            return(as.numeric(result2[[threshold.method]]))
        }else{
            result <- dismo::threshold(eval, sensitivity=threshold.sensitivity)        
            if (threshold.method == "threshold.min") {
                t1 <- result[["spec_sens"]]
                t2 <- result[["equal_sens_spec"]]            
                t3 <- result[["prevalence"]]
                thresholds <- as.numeric(c(t1, t2, t3))
                thresholds <- thresholds[thresholds > 0]
                return(min(thresholds))
            }
            if (threshold.method == "threshold.mean") {
                t1 <- result[["spec_sens"]]
                t2 <- result[["equal_sens_spec"]]            
                t3 <- result[["prevalence"]]
                thresholds <- as.numeric(c(t1, t2, t3))
                thresholds <- thresholds[thresholds > 0]
                return(mean(thresholds))
            }
            return(result[[threshold.method]])
        }
    }
#
# get raster files
# basic assumption is that different ensemble files are named as species_ENSEMBLE_1, species_ENSEMBLE_2, ... i.e. output from ensemble.batch

    species_focus <- RASTER.species.name
    if (gsub(".", "_", RASTER.species.name, fixed=T) != RASTER.species.name) {cat(paste("\n", "WARNING: species name (", RASTER.species.name, ") contains '.'", "\n\n", sep = ""))}
    if (plot.method == "suitability") {
        ensemble.files <- list.files(path=paste(getwd(), "//ensembles", sep=""), pattern=species_focus, full.names=TRUE)
        if (length(ensemble.files) < 1) {
            cat(paste("\n", "NOTE: not meaningful to provide means as there are no raster files for this species:", RASTER.species.name, "\n", sep = ""))
            return(NULL)
        }
    }
    if (plot.method == "presence") {
        ensemble.files <- list.files(path=paste(getwd(), "//ensembles//presence", sep=""), pattern=species_focus, full.names=TRUE)
        if (length(ensemble.files) < 1) {
            cat(paste("\n", "NOTE: not meaningful to provide means as there are no raster files for this species:", RASTER.species.name, "\n", sep = ""))
            return(NULL)
        }
    }
    if (plot.method == "count") {
        ensemble.files <- list.files(path=paste(getwd(), "//ensembles//count", sep=""), pattern=species_focus, full.names=TRUE)
        if (length(ensemble.files) < 1) {
            cat(paste("\n", "NOTE: not meaningful to provide means as there are no raster files for this species:", RASTER.species.name, "\n", sep = ""))
            return(NULL)
        }
    }
    RASTER.stack.name2 <- RASTER.stack.name
    if (gsub(".", "_", RASTER.stack.name, fixed=T) != RASTER.stack.name) {cat(paste("\n", "WARNING: title of stack (", RASTER.stack.name, ") contains '.'", "\n\n", sep = ""))}
    if (RASTER.stack.name != "") {
        ensemble.files <- ensemble.files[grepl(pattern=RASTER.stack.name, x=ensemble.files)]
        RASTER.stack.name2 <- paste("_", RASTER.stack.name, sep="")
        if (length(ensemble.files) < 1) {
            cat(paste("\n", "NOTE: not meaningful to provide means as there are no raster files for this stack:", RASTER.stack.name, "\n", sep = ""))
            return(NULL)
        }
    }
    for (i in 1:length(positive.filters)) {
        ensemble.files <- ensemble.files[grepl(pattern=positive.filters[i], x=ensemble.files)]
    }
    for (i in 1:length(negative.filters)) {
        ensemble.files <- ensemble.files[grepl(pattern=negative.filters[i], x=ensemble.files) == FALSE]
    }
    if (length(ensemble.files) < 1) {
        cat(paste("\n", "NOTE: not meaningful to plot as there are no raster files", "\n", sep = ""))
        return(NULL)
    }
    cat(paste("\n", "Files used to create plots", "\n\n", sep = ""))
    print(ensemble.files)

    subtitle <- NULL

    for (i in 1:length(ensemble.files)) {
        raster.focus <- raster::raster(ensemble.files[i])

        if (length(ensemble.files) > 1) {subtitle <- ensemble.files[i]}

        if (plot.method == "suitability") {
# thresholds apply to probabilities, also plot for probabilities
            raster.focus <- raster.focus / 1000
            raster.min <- raster::minValue(raster.focus)
            raster.max <- raster::maxValue(raster.focus)
            threshold.mean <- threshold
            if (threshold.mean < 0) {
                eval1 <- NULL
                if (ncol(p) == 3) {p <- p[p[,1]==species_focus, c(2:3)]}
                if (ncol(a) == 3) {a <- a[a[,1]==species_focus, c(2:3)]}
                cat(paste("\n", "Evaluation of suitability raster layer at locations p and a", "\n", sep = ""))
                cat(paste("Note that threshold is only meaningful for calibration stack suitabilities", "\n\n", sep = ""))
                pres_consensus <- raster::extract(raster.focus, p)
                abs_consensus <- raster::extract(raster.focus, a)
                eval1 <- dismo::evaluate(p=pres_consensus, a=abs_consensus)
                print(eval1)
                threshold.mean <- threshold2(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=pres_consensus, Abs=abs_consensus)
                cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
                print(as.numeric(threshold.mean))
            }
            seq1 <- round(seq(from=raster.min, to=threshold.mean, length.out=abs.breaks), 4)
            seq1 <- seq1[1:(abs.breaks-1)]
            seq1[-abs.breaks]
            seq2 <- round(seq(from = threshold.mean, to = raster.max, length.out=pres.breaks), 4)
            if (dev.new.width > 0 && dev.new.height > 0) {grDevices::dev.new(width=dev.new.width, height=dev.new.height)}
            raster::plot(raster.focus, breaks = c(seq1, seq2), col = c(grDevices::rainbow(n=(abs.breaks-1), start = 0, end = 1/6), grDevices::rainbow(n=(pres.breaks-1), start=3/6, end=4/6)), colNA = NA, 
                legend.shrink=0.8, cex.axis=0.8, main=main, sub=subtitle, ...)
        }

        if (plot.method == "presence") {
            if (dev.new.width > 0 && dev.new.height > 0) {grDevices::dev.new(width=dev.new.width, height=dev.new.height)}
            raster::plot(raster.focus, breaks=c(0, 0.5, 1), col = c("grey", "green"), colNA = NA, 
                legend.shrink=0.6, cex.axis=0.8, lab.breaks=c("", "abs", "pres"), main=main, sub=subtitle, ...)
        }

        if (plot.method == "count") {
            nmax <- raster::maxValue(raster.focus)
            if (dev.new.width > 0 && dev.new.height > 0) {grDevices::dev.new(width=dev.new.width, height=dev.new.height)}
            if (nmax > 3) {
                raster::plot(raster.focus, breaks=seq(from=-1, to=nmax, by=1), col=c("grey", grDevices::rainbow(n=(nmax-1), start=0, end=1/3), "blue"), 
                    legend.shrink=0.8, cex.axis=0.8, main=main, sub=subtitle, ...)
            }else{
                raster::plot(raster.focus, breaks=seq(from=-1, to=nmax, by=1), col=c("grey", grDevices::rainbow(n=nmax, start=0, end=1/3)), 
                    legend.shrink=0.8, cex.axis=0.8, main=main, sub=subtitle, ...)
            }
        }

    }
}
