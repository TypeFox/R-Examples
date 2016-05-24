enirg.predict <-
function(enirg.results, qtegv.maps = NULL, qlegv.maps = NULL, load.map = FALSE, 
    method = "normal", prediction.name = "predicted") {
    if (class(enirg.results) != "enirg") 
        stop("The function predict.enirg needs an object of class 'enirg'!")
    egv.map.names <- c(qlegv.maps, qtegv.maps)
    li.map.names <- paste(enirg.results$species, "_li_", colnames(enirg.results$co), sep="")
    pred.map.name <- paste(enirg.results$species, "_", prediction.name, "_pred", sep = "")
    hsm.map.name <- paste(enirg.results$species, "_", prediction.name, "_hsm", sep = "")
    HSM <- list()
    if (!is.null(egv.map.names)) {
        if (length(egv.map.names) != length(enirg.results$egvs))
            stop("Number of prediction maps does not match number of maps used to estimate the niche!")
        cat("\n\nCalculating li prediction maps from new set of EGVs ...\n\n")
        li.map.names <- paste(enirg.results$species, "_", prediction.name,"_li_Mar", sep="")
        for(i in 1:enirg.results$nf) {
            li.map.names <- c(li.map.names, paste(enirg.results$species, "_", prediction.name,"_li_Spec", i, sep=""))
        }
        li.calc <- c()
        for (j in 1:ncol(enirg.results$co)) {
            pre1 <- paste(egv.map.names, enirg.results$co[, j], sep = " * ")
            calc <- paste(li.map.names[j], "= 0", sep = "")
            for (i in 1:length(pre1)) calc <- paste(calc, pre1[i], sep = " + ")
            li.calc <- c(li.calc, calc)
        }
        if (method == "normal") {
            for(i in 1:length(li.calc)) execGRASS("r.mapcalc", expression = li.calc[i], flags = "overwrite")
        }
        if (method == "large") {
            for(i in 1:length(li.calc)) system(paste("r.mapcalc --overwrite", li.calc[i], sep = ""))
        }
    }
    HSM[["EGV_map_names"]] <- egv.map.names
    HSM[["li_map_names"]] <- li.map.names
    HSM[["HSM_map_name"]] <- hsm.map.name
    cat("\n\nCalculating Mahalanobis distances ...\n\n")
    tmp <- as.matrix(enirg.results$obs.li[, c(3:(enirg.results$nf + 3), 2)])
    mat.presences <- matrix(nrow = sum(tmp[, enirg.results$nf + 2]), ncol = ncol(enirg.results$co))
    colnames(mat.presences) <- colnames(enirg.results$co)
    samples <- 1:nrow(tmp)
    position <- 1
    for (i in samples) {
        variables <- c(tmp[i, 1:(enirg.results$nf + 1)])
        rep <- tmp[i, enirg.results$nf + 2]
        ocupa <- (position):(position + rep - 1)
        for (j in ocupa) {
            mat.presences[j, 1:(enirg.results$nf + 1)] <- variables
        }
        position <- position + rep
    }
    mean.centre <- apply(mat.presences, 2, median)
    cov <- t(as.matrix(mat.presences)) %*% as.matrix(mat.presences)/nrow(mat.presences)
    cova <- solve(cov)
    cat(paste("\n\nCalculating the Habitat Suitability Map ", enirg.results$species, "_", prediction.name,
        "_hsm ...\n\n", sep = ""))
    f3 <- function(i) enirg.results$co[, i]/sqrt(crossprod(enirg.results$co[, i])/length(enirg.results$co[, 
        i]))
    c1 <- matrix(unlist(lapply(1:(enirg.results$nf + 1), f3)), length(enirg.results$egvs))
    li_number <- length(li.map.names)
    c1 <- paste("(", li.map.names, "-", mean.centre, ")", sep = "")
    c <- (paste(c1, "*", cova, sep = ""))
    c3 <- NULL
    for (i in 1:(enirg.results$nf + 1)) {
        c3[i] <- ""
        for (j in 1:enirg.results$nf) c3[i] <- paste(c3[i], c[j + (i - 1) * (enirg.results$nf + 
            1)], "+", sep = "")
        j <- j + 1
        c3[i] <- paste(c3[i], c[j + (i - 1) * (enirg.results$nf + 1)], sep = "")
    }
    calc <- paste("(", c3[1], ")", "*", c1[1], sep = "")
    i <- 1
    if (length(c3) > 2) 
        for (i in 2:(length(c3) - 1)) calc <- paste(calc, "+", "(", c3[i], ")", 
            "*", c1[i], sep = "")
    calc <- paste(calc, "+", "(", c3[length(c3)], ")", "*", c1[i + 1], sep = "")
    calc.pred <- paste(pred.map.name, " = ", calc, sep = "")
    if (method == "normal") {
        execGRASS("r.mapcalc", expression = calc.pred, flags = "overwrite", legacyExec = TRUE)
        statistic <- execGRASS("r.univar", map = pred.map.name, flags = "g", intern = TRUE, legacyExec = TRUE)
        map.max <- as.numeric(strsplit(statistic[agrep(statistic, pattern = "max")], 
            "=")[[1]][[2]])
        calc.hsm <- paste(hsm.map.name, "=1-(", pred.map.name, "/", map.max, ")", sep = "")
        execGRASS("r.mapcalc", expression = calc.hsm, flags = "overwrite", legacyExec = TRUE)
        if (load.map == TRUE) {
            HSM[[hsm.map.name]] <- raster(readRAST(hsm.map.name))
        }
    }
    if (method == "large") {
        calc.pred <- paste("r.mapcalc --overwrite '", calc.pred, "'", sep = "")
        system(calc.pred)
        statistic <- execGRASS("r.univar", map = pred.map.name, flags = "g", intern = TRUE, legacyExec = TRUE)
        map.max <- as.numeric(strsplit(statistic[agrep(statistic, pattern = "max")], 
            "=")[[1]][[2]])
        calc.hsm <- paste("r.mapcalc --overwrite '", hsm.map.name, " = 1-(", pred.map.name, "/", map.max + 0.01, ")'", sep = "")
        system(calc.hsm)
    }
    cat(paste("'", hsm.map.name, "' HSM map was sucessfully created!\n\n", sep = ""))
    cat("You can find HSM in your mapset in GRASS. You can see it through QGIS or using 'raster' library instead.\n\n")
    presences_map <- SpatialPointsDataFrame(enirg.results$presences, data = data.frame(presences = enirg.results$presences[, 3]))
    pred.points <- paste(enirg.results$species, "_", prediction.name, sep = "")
    if(pred.points %in% list.maps()$vector) {
        cat(paste("'", pred.points, "' map is already present in your mapset!\n"))
        cat(paste("Removing '", pred.points, "' map ...\n\n"))
        execGRASS("g.remove", type = "vector",
            name = pred.points, flags = "f")
    }
    writeVECT(SDF = presences_map, vname = pred.points,
        v.in.ogr_flags = c("overwrite", "o"))
    execGRASS("v.db.addcolumn", map = pred.points, columns = "pred double precision")
    execGRASS("v.what.rast", map = pred.points, raster = hsm.map.name, column = "pred")
    vect <- readVECT(pred.points)
    vect <- na.exclude(data.frame(cbind(vect@data$presences, vect@data$pred)))
    names(vect) <- c("observed", "predicted")
    HSM[["predictions"]] <- vect
    hsm.report <- strsplit(execGRASS("r.report", map = hsm.map.name,
        nsteps = 10, units = "c", intern = T)[16:25], split = " ")
    hist.hsm <- c()
    for (i in 1:10) {
        hist.hsm <- c(hist.hsm, as.numeric(gsub("([0-9]+).*$", "\\1", hsm.report[[i]][length(hsm.report[[i]])])))
    }
    intervals <- seq(0, 1, by = 0.1)
    hist.predicted <- hist(vect$predicted, breaks = intervals, plot = FALSE)$counts
    HSM[["validation"]] <- cbind(intervals[-11], intervals[-1], hist.hsm, hist.predicted)
    colnames(HSM[["validation"]])[1:2] <- c("int.inf", "int.sup")
    class(HSM) <- c("hsm")
    return(HSM)
}
