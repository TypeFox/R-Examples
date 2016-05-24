MapRichness <- function(x, areanames = NA, leg = "continuous", show.occ = F, 
                        lin.col = "grey", lwd = 1, ...) {
    if (!class(x) == "spgeoOUT") {
        stop("this function is only defined for class \"spgeoOUT\"")
    }
    
    num <- data.frame(poly = colnames(x$polygon_table), spec.num = t(x$polygon_table), row.names = NULL, stringsAsFactors = F)
    num$poly <- gsub("[[:punct:]]", " ", as.character(num$poly))
    
    if(all(substring(num$poly, 1, 1) == "X")){
      num$poly <- gsub("X", " ", as.character(num$poly))
    }
    
    if (class(x$polygons) == "SpatialPolygonsDataFrame") {
        if (is.na(areanames)) {
          areanames <- x$areanam
        }else{
          if (!areanames %in% names(x$polygons@data)) {
            stop("cannot find areanames column; check spelling")
         }
        }
        nam.test <- as.vector(unlist(x$polygons@data[areanames]))
        if (length(grep("NA", nam.test)) > 0) {
            msg <- sprintf("The polygondata contain a polygon name %s ;this can cause problems; please rename", dQuote("NA"))
            stop(msg)
        }
        liste <- unique(x$polygons@data[, areanames])
        if (all(!is.na(liste)) == F) {
            stop("area names in polygondata contain missing data (#N/A)")
        }
        nam <- data.frame(ECO_NAME = unique(x$polygons@data[, areanames]))
        nam$ECO_NAME <- gsub("[[:punct:]]", " ", nam$ECO_NAME)
        
        if (suppressWarnings(all(!is.na(as.numeric(nam$ECO_NAME))))) {
            nam$ECO_NAME <- as.numeric(nam$ECO_NAME)
            num$poly <- as.numeric(num$poly)
            warning("areanames are numeric")
        }
        polys.df <- merge(nam, num, sort = F, by.x = "ECO_NAME", by.y = "poly", all = T)
    }
    
    if (class(x$polygons) == "SpatialPolygons") {
        if (length(grep("NA", names(x$polygons))) > 0) {
            stop("the polygondata contain a polygon named \"NA\"; please rename")
        }
        x$polygons <- SpatialPolygonsDataFrame(x$polygons, data.frame(names(x$polygons), row.names = names(x$polygons)))
        names(x$polygons) <- c("ECO_NAME")
        nam <- data.frame(ECO_NAME = x$polygons$ECO_NAME)
        
        polys.df <- merge(nam, num, sort = F, by.x = "ECO_NAME", by.y = "poly", all = T)
        areanames <- "ECO_NAME"
    }
    
    polys.df$spec.num[is.na(polys.df$spec.num)] <- 0
    
    if (leg == "continuous" | leg == "continous") {
        colo <- data.frame(num = c(0:max(polys.df$spec.num)), code = c("#FFFFFFFF", rev(heat.colors(max(polys.df$spec.num)))))
        
    } else {
        colo <- data.frame(num = c(0, sort(unique(polys.df$spec.num))), code = c("#FFFFFFFF", rev(rainbow(length(unique(polys.df$spec.num))))))
        if (colo$num[2] == 0) {
            colo <- colo[-2, ]
        }
    }
    
    polys.df <- merge(polys.df, colo, sort = F, by.x = "spec.num", by.y = "num", all = F)
    polys.df$ord <- pmatch(polys.df$ECO_NAME, nam$ECO_NAME)
    polys.df <- polys.df[order(polys.df$ord), ]
    names(polys.df)[1] <- "sp.count"
    
    dum1 <- x$polygons@data
    dum1$ident.add <- seq(1, dim(dum1)[1])
    dum1[areanames] <- gsub("[[:punct:]]", " ", as.vector(unlist(dum1[areanames])))
    
    if (suppressWarnings(all(!is.na(as.numeric(as.vector(unlist(dum1[areanames]))))))) {
        polys.df$ECO_NAME <- as.numeric(polys.df$ECO_NAME)
        dum1[areanames] <- as.numeric(as.vector(unlist(dum1[areanames])))
    }
    
    dum2 <- merge(dum1, polys.df, by.x = areanames, by.y = "ECO_NAME", sort = F, all = T)
    dum2 <- dum2[order(dum2$ident.add), ]  #,c('sp.count','code','ord')]
    x$polygons@data <- cbind(x$polygons@data, dum2)
    plotpoly <- x$polygons
#     
#     if (lim == "polygons") {
#         limits <- bbox(plotpoly)
#         if (limits[1, 1] < -170 && limits[1, 2] > 170 && max(x$species_coordinates$XCOOR) < -10) {
#             limits[1, 2] <- -10
#         }
#         if (limits[1, 1] < -170 && limits[1, 2] > 170 && min(x$species_coordinates$XCOOR) > 0) {
#             limits[1, 1] <- 0
#         }
#         limits[1, 1] <- max(limits[1, 1] - abs(abs(limits[1, 1]) - abs(limits[1, 2])) * 0.2, -180)
#         limits[1, 2] <- min(limits[1, 2] + abs(abs(limits[1, 1]) - abs(limits[1, 2])) * 0.2, 180)
#         limits[2, 1] <- max(limits[2, 1] - abs(abs(limits[2, 1]) - abs(limits[2, 2])) * 0.2, -90)
#         limits[2, 2] <- min(limits[2, 2] + abs(abs(limits[2, 1]) - abs(limits[2, 2])) * 0.2, 90)
#     }
#     if (lim == "points") {
#         limits <- matrix(ncol = 2, nrow = 2)
#         limits[1, 1] <- max(min(x$species_coordinates_in$XCOOR) - abs(min(x$species_coordinates_in$XCOOR) - max(x$species_coordinates_in$XCOOR)) * 
#             0.2, -180)
#         limits[1, 2] <- min(max(x$species_coordinates_in$XCOOR) + abs(abs(min(x$species_coordinates_in$XCOOR)) - abs(max(x$species_coordinates_in$XCOOR))) * 
#             0.2, 180)
#         limits[2, 1] <- max(min(x$species_coordinates_in$YCOOR) - abs(min(x$species_coordinates_in$YCOOR) - max(x$species_coordinates_in$YCOOR)) * 
#             0.2, -90)
#         limits[2, 2] <- min(max(x$species_coordinates_in$YCOOR) + abs(abs(min(x$species_coordinates_in$YCOOR)) - abs(max(x$species_coordinates_in$YCOOR))) * 
#             0.2, 90)
#     }
#     
        if (length(unique(plotpoly$sp.count)) == 1) {
        leg <- "discrete"
    }
    layout(matrix(c(2, 2, 2, 2, 2, 1), ncol = 6, nrow = 1))
    
    if (leg == "continuous" | leg == "continous") {
        par(mar = c(5, 1, 5, 3))
        ifelse(max(plotpoly@data$sp.count) < 25, leng <- max(plotpoly$sp.count), leng <- 11)
        ticks <- round(seq(min(plotpoly$sp.count), max(plotpoly$sp.count), len = leng), 0)
        scale <- (length(colo$num) - 1)/(max(plotpoly$sp.count) - min(plotpoly$sp.count))
        plot(c(0, 10), c(min(plotpoly$sp.count), max(plotpoly$sp.count)), type = "n", bty = "n", xaxt = "n", xlab = "", yaxt = "n", 
            ylab = "")
        box("plot")
        title("Species\nnumber")
        axis(4, ticks, las = 1)
        for (i in 1:(length(colo$num))) {
            y = (i - 1)/scale + min(plotpoly$sp.count)
            rect(0, y, 10, y + 1/scale, col = as.character(colo$code[i]), border = NA)
        }
    }
    if (leg == "discrete") {
        par(mar = c(5, 1, 5, 3))
        plot(c(0, 10), c(0, dim(colo)[1] + 1), type = "n", bty = "n", xaxt = "n", xlab = "", yaxt = "n", ylab = "", 
            xpd = F)
#         plot(c(0, 10), c(min(colo$num), dim(colo)[1] + 1), type = "n", bty = "n", xaxt = "n", xlab = "", yaxt = "n", ylab = "", 
#              xpd = F)
        title("Species\nnumber")        
        if (length(unique(plotpoly$sp.count)) == 1) {
            rect(0, (dim(colo)[1] + 1)/2 - 3, 5, (dim(colo)[1] + 1)/2 - 1, col = "white", border = "black")
            rect(0, (dim(colo)[1] + 1)/2 + 1, 5, (dim(colo)[1] + 1)/2 + 3, col = as.character(colo$code[unique(plotpoly$sp.count)]), 
                border = NA)
            text(7, (dim(colo)[1] + 1)/2 - 2, labels = "0")
            text(7, (dim(colo)[1] + 1)/2 + 2, labels = unique(plotpoly$sp.count))
        } else {
            y <- seq(0, (dim(colo)[1] - 1), by = 1)
            for (i in 1:length(colo$num)) {
                rect(0, y[i], 5, y[i] + 1, col = as.character(colo$code[i]), border = NA)
                text(8, y[i] + 0.5, colo$num[i])
            }
            rect(0, 0, 5, (max(y) + 1))
        }
    }
#     if(lim == "polygons" | lim == "points"){
#       map("world", xlim = limits[1, ], ylim = limits[2, ])
#     }else{
      map("world", ...) 
#     }
    axis(1)
    axis(2)
    box("plot")
    plot(plotpoly, col = as.character(plotpoly@data$code), border = lin.col, add = T, lwd = lwd)  #, ...)
    
    if (show.occ == T) {
        points(x$species_coordinates_in$XCOOR, x$species_coordinates_in$YCOOR)
    }

    layout(matrix(1, 1, 1))
    return(plotpoly)
} 
