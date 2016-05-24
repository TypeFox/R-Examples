enirg <-
function(presences.table, qtegv.maps, qlegv.maps = NULL, col.w = NULL, scannf = TRUE, 
    nf = 1, method = "normal", load.maps = FALSE, species.name = "species", map.center = NULL, res = NULL) {
    egv.maps <- c(qtegv.maps, qlegv.maps)
    number_maps <- length(egv.maps)
    if(!is.null(map.center)) execGRASS("r.mask", raster = map.center)
    if(is.null(map.center)) map.center <- egv.maps[1]
    cat(paste("\n\nStudy region will be set according to '", map.center, "' map settings:", sep=""))
    execGRASS("g.region", raster = map.center)
    if(!is.null(res)) execGRASS("g.region", res = res, flags = "a")
    region.parameters <- gmeta()
    cat(paste("\nNorth-South resolution set to ", region.parameters$nsres, sep=""))
    cat(paste("\nEast-West resolution set to ", region.parameters$ewres, sep=""))
    null_cells <- NULL
    for(i in egv.maps) {
        null.in.map <- execGRASS("r.univar", map = i, flags = "g", intern = TRUE)
        null.in.map <- agrep(null.in.map, pattern = "null_cells", max.distance = list(all = 0), 
            value = T)
        null.in.map <- as.numeric(strsplit(null.in.map, "=")[[1]][2])
        null_cells <- c(null_cells, null.in.map)
    }
    if(nlevels(factor(null_cells)) != 1) cat("\n\nWARNING: Map information does not seem to match\nThis could cause enirg function to fail.\n\nPlease, check resolution, extension and null values in your data set\n\nStarting calculations ...")
    numberpixel <- as.numeric(region.parameters$cells) - null_cells[1]
    cat(paste("\n\nThe number of pixels for the calculated area is", numberpixel, "\n\n", 
        sep = " "))
    cat(paste("The number of null pixels for the calculated area is", null_cells[1], "\n\n", 
        sep = " "))
    if(species.name %in% list.maps()$vector) execGRASS("g.remove", type = "vector", name = species.name, flags = c("f", "quiet"))
    presences_map <- SpatialPointsDataFrame(presences.table, data = data.frame(presences = presences.table[, 
        3]))
    writeVECT(SDF = presences_map, vname = species.name, v.in.ogr_flags = c("overwrite", "o", "quiet"))
    execGRASS("v.to.rast", input = species.name, output = species.name, type= "point", use = "attr", attribute_column = "presences",
        flags = c("overwrite", "quiet"))
    cat(paste(species.name, " was sucessfully convert into a raster map and loaded into GRASS ...\n\n", 
        sep = ""))
    presences_map <- species.name
    number_presences <- sum(presences.table[, 3])
    combine_maps <- t(combn(1:number_maps, 2))
    map.cor <- cbind(egv.maps[combine_maps[, 1]], egv.maps[combine_maps[, 2]])
    combine_maps2 <- rbind(t(combn(1:number_maps, 2)), cbind(1:number_maps, 1:number_maps))
    map.cor2 <- cbind(egv.maps[combine_maps2[, 1]], egv.maps[combine_maps2[, 2]])
    if (method == "normal") {
        cat("\n\nCalculating the species covariance matrix (Rs) ...\n\n")
        wrapping <- function(tmp) as.numeric(strsplit(tmp, " ")[[1]])
        grass.cor <- NULL
        for (i in 1:nrow(map.cor)) {
            tmp <- execGRASS("r.stats", input = paste(map.cor[i, 1], map.cor[i, 2], 
                presences_map, sep = ","), flags = c("1", "A", "n", "N", "quiet"), intern = TRUE, legacyExec = TRUE)
            tmp <- t(apply(matrix(tmp, ncol = 1), 1, wrapping))
            tmp <- tmp[which(tmp[, 3] > 0), ]
            np <- sum(tmp[, 3])
            tmp <- sum(apply(tmp, 1, prod))
            grass.cor[i] <- tmp/np
        }
        grass.cor2 <- NULL
        for (i in 1:number_maps) {
            tmp <- execGRASS("r.stats", input = paste(egv.maps[i], egv.maps[i], presences_map, 
                sep = ","), flags = c("1", "A","n","N","quiet"), intern = TRUE, legacyExec = TRUE)
            tmp <- t(apply(matrix(tmp, ncol = 1), 1, wrapping))
            tmp <- tmp[which(tmp[, 3] > 0), ]
            np <- sum(tmp[, 3])
            tmp <- sum(apply(tmp, 1, prod))
            grass.cor2[i] <- tmp/np
        }
        Rs <- matrix(1, number_maps, number_maps)
        diag(Rs) <- grass.cor2
        Rs[t(combn(1:number_maps, 2))] <- grass.cor
        Rs[t(combn(1:number_maps, 2))[, 2:1]] <- grass.cor
        rownames(Rs) <- colnames(Rs) <- egv.maps
        cat("\n\nCalculating Ze ...\n\n")
        grass.cor = NULL
        for (i in 1:nrow(map.cor)) {
            tmp <- execGRASS("r.stats", input = paste(map.cor[i, 1], map.cor[i, 2], 
                sep = ","), flags = c("1", "A","n","N","quiet"), intern = TRUE, legacyExec = TRUE)
            tmp <- t(apply(matrix(tmp, ncol = 1), 1, wrapping))
            grass.cor[i] <- sum(apply(tmp, 1, prod))
        }
        Ze <- matrix(1, number_maps, number_maps)
        Ze[t(combn(1:number_maps, 2))] <- grass.cor
        Ze[t(combn(1:number_maps, 2))[, 2:1]] <- grass.cor
        Ze <- Ze/numberpixel
        rownames(Ze) <- colnames(Ze) <- egv.maps
        cat("\n\nCalculating the coordinates of the marginality vector ...\n\n")
        mar <- NULL
        for (i in egv.maps) {
            tmp <- execGRASS("r.stats", input = paste(i, presences_map, sep = ","), 
                flags = c("1","n","N","quiet"), intern = TRUE, legacyExec = TRUE)
            tmp <- (t(apply(matrix(tmp, ncol = 1), 1, wrapping)))
            tmp <- tmp[which(tmp[, 2] > 0), ]
            np <- sum(tmp[, 2])
            tmp <- apply(tmp, 1, prod)
            mar[i] <- sum(tmp)/np
        }
        names(mar) <- egv.maps
        cat("\n\nCalculating the matrix Ge ...\n\n")
        grass.cor <- NULL
        for (i in 1:nrow(map.cor)) {
            tmp <- execGRASS("r.stats", input = paste(map.cor[i, 1], map.cor[i, 2], 
                sep = ","), flags = c("1","n","N","quiet"), intern = TRUE, legacyExec = TRUE)
            tmp <- t(apply(matrix(tmp, ncol = 1), 1, wrapping))
            tmp <- sum(apply(tmp, 1, prod))
            grass.cor[i] <- tmp/numberpixel
        }
        Ge <- matrix(1, number_maps, number_maps)
        Ge[t(combn(1:number_maps, 2))] <- grass.cor
        Ge[t(combn(1:number_maps, 2))[, 2:1]] <- grass.cor
        rownames(Ge) <- colnames(Ge) <- egv.maps
        cat("\n\nCalculating the matrix Se ...\n\n")
        for (i in egv.maps) {
            tmp <- execGRASS("r.mapcalc", expression = paste("pond_", i, "=", i, "*",
                presences_map, "/", number_presences, sep = ""),
                flags = c("overwrite", "quiet"), legacyExec = TRUE)
        }
        grass.cor <- NULL
        for (i in 1:nrow(map.cor2)) {
            tmp <- execGRASS("r.stats", input = paste(map.cor2[i, 1], ",pond_", map.cor2[i, 
                2], sep = ""), flags = c("1", "A","n","N","quiet"), intern = TRUE, legacyExec = TRUE)
            tmp <- t(apply(matrix(tmp, ncol = 1), 1, wrapping))
            grass.cor[i] <- sum(apply(tmp, 1, prod))
        }
        Se <- matrix(1, number_maps, number_maps)
        Se[combine_maps2] <- grass.cor
        Se[combine_maps2[, 2:1]] <- grass.cor
        rownames(Se) <- colnames(Se) <- egv.maps
    }
    if (method == "large") {
        cat("\n\nCalculating the species covariance matrix (Rs) ...\n\n")
        grass.cor.map <- NULL
        for (i in 1:nrow(map.cor)) grass.cor.map[i] <- paste("r.stats -1 -A ", map.cor[i, 
            1], map.cor[i, 2], presences_map, sep = ",")
        grass.cor.p <- NULL
        for (i in 1:nrow(map.cor)) grass.cor.p[i] <- as.numeric(system(paste(grass.cor.map[i], 
            " | awk '$3>0 {sum=sum+($1*$3*$2);np=np+$3} END{print sum/np}'"), intern = T))
        grass.cor.map2 <- NULL
        for (i in 1:number_maps) grass.cor.map2[i] <- paste("r.stats -1 -A ", egv.maps[i], 
            egv.maps[i], presences_map, sep = ",")
        grass.cor.p2 <- NULL
        for (i in 1:number_maps) grass.cor.p2[i] <- as.numeric(system(paste(grass.cor.map2[i], 
            " | awk '$3>0 {sum=sum+($1*$3*$2);np=np+$3} END{print sum/np}'"), intern = T))
        Rs <- matrix(1, number_maps, number_maps)
        diag(Rs) <- grass.cor.p2
        Rs[t(combn(1:number_maps, 2))] <- grass.cor.p
        Rs[t(combn(1:number_maps, 2))[, 2:1]] <- grass.cor.p
        rownames(Rs) <- colnames(Rs) <- egv.maps
        cat("\n\nCalculating Ze ...\n\n")
        grass.cor.map <- NULL
        for (i in 1:nrow(map.cor)) grass.cor.map[i] <- paste("r.stats -1 -A ", map.cor[i, 
            1], map.cor[i, 2], sep = ",")
        grass.cor.p <- NULL
        for (i in 1:nrow(map.cor)) grass.cor.p[i] <- as.numeric(system(paste(grass.cor.map[i], 
            " | awk '{sum=sum+($1*$2)} END{print sum}'"), intern = T))
        Ze <- matrix(1, number_maps, number_maps)
        Ze[t(combn(1:number_maps, 2))] <- grass.cor.p
        Ze[t(combn(1:number_maps, 2))[, 2:1]] <- grass.cor.p
        Ze <- Ze/numberpixel
        rownames(Ze) <- colnames(Ze) <- egv.maps
        cat("\n\nCalculating the coordinates of the marginality vector ...\n\n")
        mar <- NULL
        for (i in egv.maps) mar[i] <- as.numeric(system(paste("r.stats -1 ", i, ",", 
            presences_map, " | awk '{s=s+$1*$2;n=n+$2} END{print s/n}'", sep = ""), 
            intern = T))
        names(mar) <- egv.maps
        cat("\n\nCalculating the matrix Ge ...\n\n")
        grass.cor.map <- NULL
        for (i in 1:nrow(map.cor)) grass.cor.map[i] <- paste("r.stats -1 ", map.cor[i, 
            1], map.cor[i, 2], sep = ",")
        grass.cor.p <- NULL
        for (i in 1:nrow(map.cor)) grass.cor.p[i] <- as.numeric(system(paste(grass.cor.map[i], 
            " | awk '{sum=sum+($1*$2)/", numberpixel, "} END{print sum}'"), intern = T))
        Ge <- matrix(1, number_maps, number_maps)
        Ge[t(combn(1:number_maps, 2))] <- grass.cor.p
        Ge[t(combn(1:number_maps, 2))[, 2:1]] <- grass.cor.p
        rownames(Ge) <- colnames(Ge) <- egv.maps
        cat("\n\nCalculating the matrix Se ...\n\n")
        for (i in egv.maps) system(paste("r.mapcalc --overwrite 'pond_", i, " = ", i, " *", presences_map, 
            " / ", number_presences, "'", sep = ""))
        grass.cor.map <- NULL
        for (i in 1:nrow(map.cor2)) grass.cor.map[i] <- paste("r.stats -1A ", map.cor2[i, 
            1], ",pond_", map.cor2[i, 2], sep = "")
        grass.cor.p <- NULL
        for (i in 1:nrow(map.cor2)) grass.cor.p[i] <- as.numeric(system(paste(grass.cor.map[i], 
            " | awk '{sum=sum+($1*$2)} END{print sum}'"), intern = T))
        Se <- matrix(1, number_maps, number_maps)
        Se[combine_maps2] <- grass.cor.p
        Se[combine_maps2[, 2:1]] <- grass.cor.p
        rownames(Se) <- colnames(Se) <- egv.maps
    }
    eS <- eigen(Se)
    kee <- (eS$values > 1e-09)
    S12 <- eS$vectors[, kee] %*% diag(eS$values[kee]^(-0.5)) %*% t(eS$vectors[, kee])
    W <- S12 %*% Ge %*% S12
    if (is.null(col.w)) 
        col.w <- rep(1, number_maps)
    me <- mar * sqrt(col.w)
    x <- S12 %*% me
    b <- x/sqrt(sum(x^2))
    H <- (diag(ncol(Ze)) - b %*% t(b)) %*% W %*% (diag(ncol(Ze)) - b %*% t(b))
    s <- eigen(H)$values[-number_maps]
    if (scannf) {
        barplot(s)
        cat("Select the number of specialization axes: ")
        nf <- as.integer(readLines(n = 1))
    }
    if (nf <= 0 | nf > (ncol(Ze) - 1) | is.null(nf)) {
        nf <- 1
        cat("\n\nWARNING!\n")
        cat("Number of specialization axes can not be:\n")
        cat("* 0\n")
        cat("* larger than (EGVs - 1)\n\n")
        cat("nf set to 1 ...\n\n")
    }
    co <- matrix(nrow = number_maps, ncol = nf + 1)
    tt <- data.frame((S12 %*% eigen(H)$vectors)[, 1:nf])
    ww <- apply(tt, 2, function(x) x/sqrt(col.w))
    norw <- sqrt(diag(t(as.matrix(tt)) %*% as.matrix(tt)))
    co[, 2:(nf + 1)] <- sweep(ww, 2, norw, "/")
    m <- me/sqrt(col.w)
    co[, 1] <- m/sqrt(sum(m^2))
    m <- sum(m^2)
    marginalities <- matrix(co[, 1], ncol = 1)
    rownames(marginalities) <- egv.maps
    colnames(marginalities) <- "Marginality"
    total_marginality <- sqrt(sum((marginalities/sqrt(col.w))^2))/1.96
    rownames(co) <- egv.maps
    colnames(co) <- c("Mar", paste("Spec", 1:nf, sep = ""))
    specializations <- matrix(co[, 2:(nf + 1)], ncol = nf)
    rownames(specializations) <- egv.maps
    colnames(specializations) <- paste("Spec", 1:nf, sep = "")
    total_specialization <- sqrt(sum(abs(s)))/number_maps
    call <- match.call()
    results_ENFA <- list(call = call, nf = nf, cw = col.w, species = species.name, 
        egvs = egv.maps, qt.egvs = qtegv.maps, ql.egvs = qlegv.maps, presences = presences.table, 
        total.marginality = total_marginality, marginalities = marginalities, total.specialization = total_specialization, 
        specializations = specializations, co = co, mar = mar, m = m, s = s)
    class(results_ENFA) <- c("enirg")
    cat("\n\nCalculating li maps ...\n\n")
    if (method == "normal") {
        for (j in 1:(nf + 1)) {
            execGRASS("r.mapcalc", expression = paste(species.name, "_li_", colnames(co)[j], "=0", sep = ""),
                flags = c("overwrite", "quiet"), legacyExec = TRUE)
            for (i in 1:number_maps) {
                execGRASS("r.mapcalc", expression = paste("temp=", egv.maps[i], "*", results_ENFA$co[i, j], sep = ""),
                    flags = c("overwrite", "quiet"), legacyExec = TRUE)
                execGRASS("r.mapcalc", expression = paste(species.name, "_li_", colnames(co)[j], "=", species.name, "_li_",
                    colnames(co)[j], "+temp", sep = ""), flags = c("overwrite", "quiet"), legacyExec = TRUE)

            }
        }
    }
    if (method == "large") {
        for (j in 1:(nf + 1)) {
            cc <- paste("(", co[, j], ")", sep = "")
            pre1 <- paste(egv.maps, cc, sep = " * ")
            calc <- paste("r.mapcalc --overwrite '", paste(species.name, "_li_", colnames(co)[j], sep = ""), 
                " = 0", sep = "")
            for (i in 1:(number_maps - 1)) calc <- paste(calc, pre1[i], sep = " + ")
            calc <- paste(calc, " + ", pre1[number_maps], "'", sep = "")
            system(calc)
        }
        load.maps <- FALSE
    }
    cat("\n\nExtracting suitability information ...\n\n")
    if (load.maps) {
        for (i in 1:(nf + 1)) {
            results_ENFA[[paste(species.name, "_li_", colnames(co)[i], sep = "")]] <- raster(readRAST(paste(species.name, "_li_", 
                colnames(co)[i], sep = "")))
        }
    }
    for (i in colnames(co)) {
        execGRASS("v.db.addcolumn", map = species.name, columns = paste("li_", 
            i, " double precision", sep = ""), flags = c("quiet"))
        execGRASS("v.what.rast", map = species.name, raster = paste(species.name, "_li_", 
            i, sep = ""), column = paste("li_", i, sep = ""), flags = c("quiet"))
    }
    vectorial.data <- readVECT(species.name)
    results_ENFA[["obs.li"]] <- na.exclude(vectorial.data@data)
    execGRASS("g.remove", type = "raster", name = "temp", flags = c("f", "quiet"))
    if(!is.null(map.center)) execGRASS("r.mask", flags = c("r", "quiet"))
    cat("\n\nCalculations done!\n\n")
    return(results_ENFA)

}
