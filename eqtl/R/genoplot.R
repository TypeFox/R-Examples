genoplot <- function (peak.array, cross, etrait.coord, data.gmap, chr.size, save.pict = FALSE, ...) 
{
    require(qtl)
    if (length(class(cross)) < 2 || class(cross)[2] != "cross") stop("Input should have class \"cross\".")
    if (!any(class(peak.array) == "peak.array")) stop("peak.array should have class \"peak.array\" \"data.frame\". ")
    if (!all(c("peak.bp", "inf.bp", "sup.bp") %in% attr(peak.array, "names", exact = TRUE))) 
        stop("Argument peak.array misspecified: peak.array should have columns \"peak.bp\", \"inf.bp\" and \"sup.bp\". You should use the function \"localize.qtl\".")
    if (class(etrait.coord)[1] != "data.frame" || !all(names(etrait.coord) %in% c("etrait.name", "chr", "start", "stop"))) 
        stop("etrait.coord should have class \"data.frame\" and columns names: 'etrait.name','chr','start','stop' \n")
    if (class(data.gmap)[1] != "data.frame" || any(!(names(data.gmap) %in% c("Marker", "chr", "PP")))) 
        stop("data.gmap should have class \"data.frame\" with columns names: 'Marker','chr','PP'")
    if (missing(data.gmap)) stop("Argument data.gmap unspecified")
    if (!is.vector(chr.size) && !is.numeric(chr.size)) stop("chr.size should be a numeric vector\n")
    if (!is.logical(save.pict)) stop("save.pict must be logical: TRUE or FALSE\n")
    
    gmap.size <- NA
    
    for (i in 1:nchr(cross)) {
        last_marker <- names(cross$gen[[i]]$map[length(cross$gen[[i]]$map)])
        bool <- data.gmap$Marker %in% last_marker
        gchr.size <- data.gmap$PP[bool]
        gmap.size <- c(gmap.size, gchr.size)
    }
    gmap.size <- gmap.size[-1]
    
    limit_gst <- c(0, as.integer(cumsum(chr.size))) * 0.001
    limit_qtl <- limit_gst
    print("chromosome size:")
    print(c(0, chr.size) * 0.001)
    print("cumul:")
    print(limit_gst)
    trait <- "NA"
    coord_gst <- NA
    coord_qtl <- NA
    coord_inf <- NA
    coord_sup <- NA
    prev_p_chr.bp <- NA
    prev_g_chr.bp <- NA
    lod <- NA
    ae <- NA
    for (i in 1:nrow(peak.array)) {
        if (!is.na(peak.array$mname.peak[i]))
            bool <- toupper(etrait.coord$etrait.name) %in% toupper(peak.array$trait[i])
        else bool <- FALSE
        if (any(bool)) {
            coord_qtl <- c(coord_qtl, as.numeric(peak.array$peak.bp[i]) * 0.001)
            coord_inf <- c(coord_inf, as.numeric(peak.array$inf.bp[i]) * 0.001)
            coord_sup <- c(coord_sup, as.numeric(peak.array$sup.bp[i]) * 0.001)
            prev_g_chr.bp <- c(prev_g_chr.bp, limit_qtl[as.numeric(as.vector(peak.array$chr[i]))])
            lod <- c(lod, as.numeric(as.vector(peak.array$lod[i])))
            if (any(names(peak.array) %in% "additive.effect")) ae <- c(ae, as.numeric(as.vector(peak.array$additive.effect[i])))
        }
        coord_gst <- c(coord_gst, as.numeric(as.vector(etrait.coord$start[bool])) * 0.001)
        prev_p_chr.bp <- c(prev_p_chr.bp, limit_gst[as.numeric(as.vector(etrait.coord$chr[bool]))])
    }
    coord_gst <- coord_gst[-1]
    coord_qtl <- coord_qtl[-1]
    coord_inf <- coord_inf[-1]
    coord_sup <- coord_sup[-1]
    prev_p_chr.bp <- prev_p_chr.bp[-1]
    prev_g_chr.bp <- prev_g_chr.bp[-1]
    
    if (save.pict) png(filename = "./histogram_controled_etrait.png", width = 1280, height = 1024)
    else {
        opar <- par()
        par(mfrow = c(2, 3))
    }
    
    Break <- NA
    for (i in 1:(length(limit_gst) - 1)) Break <- c(Break, seq(limit_gst[i], limit_gst[i + 1], by = (limit_gst[i + 1] - limit_gst[i])/15))
    Break <- unique(Break[-1])
        
    suppressWarnings(
        hist(
          as.numeric(as.vector(coord_gst)) + prev_p_chr.bp,
          xlim = c(0, max(limit_gst)), main = "Genomic distribution of controlled eTRAIT", 
          bty = "u", axes = FALSE, xlab = "physical e-trait location (Mb)", 
          freq = TRUE, col = "gray")
        )

    axis(1, limit_gst, limit_gst * 0.001, col.axis = "darkblue")
    axis(1, labels = FALSE, Break)
    axis(2)
    abline(v = limit_gst, lty = 2, col = "darkblue")
    box()
    if (save.pict) 
        dev.off()
    if (save.pict) 
        png(filename = "./histogram_qtl.png", width = 1280, height = 1024)
    Break <- NA
    for (i in 1:(length(limit_qtl) - 1)) Break <- c(Break, seq(limit_qtl[i], 
        limit_qtl[i + 1], by = (limit_qtl[i + 1] - limit_qtl[i])/15))
    Break <- unique(Break[-1])
    suppressWarnings(hist(as.numeric(as.vector(coord_qtl)) + 
        prev_g_chr.bp, xlim = c(0, max(limit_qtl)), xlab = "physical eQTL location (Mb)", 
        main = "Genomic QTL distribution", bty = "u", breaks = Break, 
        col = "gray", axes = FALSE, freq = TRUE))
    axis(1, limit_qtl, limit_qtl * 0.001, col.axis = "darkblue")
    axis(1, labels = FALSE, Break)
    axis(2)
    abline(v = limit_qtl, lty = 2, col = "darkblue")
    if (save.pict) 
        dev.off()
    gqplot <- function(col, limit_qtl, limit_gst, prev_g_chr.bp, 
        prev_p_chr.bp, chr.size, file, save.pict) {
        if (save.pict) 
            png(filename = paste(file, "_dotplot_traitxqtl.png", 
                sep = ""), width = 1280, height = 1024)
        plot(as.numeric(as.vector(coord_qtl)) + prev_g_chr.bp, 
            as.numeric(as.vector(coord_gst)) + prev_p_chr.bp, 
            main = "", xlab = "physical eQTL location (Mb)", 
            ylab = "physical e-trait location (Mb)", col = col, 
            pty = "s", axes = FALSE, mar = c(1, 1, 1, 1), xlim = c(0, 
                max(limit_qtl)), ylim = c(0, max(limit_gst)), 
            ...)
        axis(1, c(seq(0, max(limit_qtl), by = 10000), max(limit_qtl)), 
            c(seq(0, max(limit_qtl), by = 10000), max(limit_qtl)) * 
                0.001, )
        for (i in 1:(length(limit_qtl) - 1)) {
            if (i == length(limit_qtl) - 1) 
                end <- gmap.size[i]
            else end <- 0
            axis(3, c(seq(limit_qtl[i], limit_qtl[i + 1], by = 10000), 
                limit_qtl[i + 1]), c(seq(0, chr.size[i] * 0.001, 
                by = 10000), end) * 0.001, col = "blue", col.axis = "blue", 
                lwd = 2)
        }
        axis(3, limit_qtl, c(rep(0, 5), max(limit_qtl)) * 0.001, 
            col = "darkblue", col.axis = "black", lwd = 2)
        abline(v = limit_qtl, lty = 2, )
        axis(2, c(seq(0, max(limit_gst), by = 10000), max(limit_gst)), 
            c(seq(0, max(limit_gst), by = 10000), max(limit_gst)) * 
                0.001)
        abline(h = limit_gst, lty = 2, )
        for (i in 1:(length(limit_gst) - 1)) {
            if (i == length(limit_gst) - 1) 
                end <- chr.size[i]
            else end <- 0
            axis(4, c(seq(limit_gst[i], limit_gst[i + 1], by = 10000), 
                limit_gst[i + 1]), c(seq(0, chr.size[i] * 0.001, 
                by = 10000), end) * 0.001, col = "red", col.axis = "red", 
                lwd = 2)
        }
        axis(4, limit_gst, c(rep(0, 5), max(limit_gst)) * 0.001, 
            col = "darkred", col.axis = "black", lwd = 2)
        box()
        if (save.pict) 
            dev.off()
        if (save.pict) 
            png(filename = paste(file, "_siplot_traitxqtl.png", 
                sep = ""), width = 1280, height = 1024)
        matplot(rbind(t(as.numeric(as.vector(coord_inf)) + prev_g_chr.bp), 
            t(as.numeric(as.vector(coord_sup)) + prev_g_chr.bp)), 
            rbind(t(as.numeric(as.vector(coord_gst)) + prev_p_chr.bp), 
                t(as.numeric(as.vector(coord_gst)) + prev_p_chr.bp)), 
            type = "l", lty = 1, main = "", xlab = "physical eQTL location (Mb)", 
            ylab = "physical e-trait location (Mb)", col = col, 
            pty = "s", axes = FALSE, mar = c(1, 1, 1, 1), xlim = c(0, 
                max(limit_qtl)), ylim = c(0, max(limit_gst)))
        axis(1, c(seq(0, max(limit_qtl), by = 10000), max(limit_qtl)), 
            c(seq(0, max(limit_qtl), by = 10000), max(limit_qtl)) * 
                0.001)
        for (i in 1:(length(limit_qtl) - 1)) {
            if (i == length(limit_qtl) - 1) 
                end <- gmap.size[i]
            else end <- 0
            axis(3, c(seq(limit_qtl[i], limit_qtl[i + 1], by = 10000), 
                limit_qtl[i + 1]), c(seq(0, chr.size[i] * 0.001, 
                by = 10000), end) * 0.001, col = "blue", col.axis = "blue", 
                lwd = 2)
        }
        axis(3, limit_qtl, c(rep(0, 5), max(limit_qtl)) * 0.001, 
            col = "darkblue", col.axis = "black", lwd = 2)
        abline(v = limit_qtl, lty = 2, )
        axis(2, c(seq(0, max(limit_gst), by = 10000), max(limit_gst)), 
            c(seq(0, max(limit_gst), by = 10000), max(limit_gst)) * 
                0.001)
        abline(h = limit_gst, lty = 2, )
        for (i in 1:(length(limit_gst) - 1)) {
            if (i == length(limit_gst) - 1) 
                end <- gmap.size[i]
            else end <- 0
            axis(4, c(seq(limit_gst[i], limit_gst[i + 1], by = 10000), 
                limit_gst[i + 1]), c(seq(0, chr.size[i] * 0.001, 
                by = 10000), end) * 0.001, col = "red", col.axis = "red", 
                lwd = 2)
        }
        axis(4, limit_gst, c(rep(0, 5), max(limit_gst)), col = "darkred", 
            col.axis = "black", lwd = 2)
        box()
        if (save.pict) 
            dev.off()
    }
    col_lod <- rainbow(length(lod[-1]), start = 0.2, end = 0)[order(lod[-1])]
    cat(max(lod[-1]), min(lod[-1]), "\n")
    gqplot(col_lod, limit_qtl, limit_gst, prev_g_chr.bp, prev_p_chr.bp, 
        chr.size, file = "lod", save.pict)
    if (length(ae) > 1) {
        col_ae <- rainbow(length(ae[-1]), start = 0, end = 2/6)[order(ae[-1])]
        cat(max(ae[-1]), min(ae[-1]), "\n")
        gqplot(col_ae, limit_qtl, limit_gst, prev_g_chr.bp, prev_p_chr.bp, 
            chr.size, file = "ae", save.pict)
    }
    if (!save.pict) 
        par(c(1, 1))
    return(list(coord_gst = as.numeric(as.vector(coord_gst)), 
        coord_qtl = as.numeric(as.vector(coord_qtl)), limit = as.numeric(as.vector(limit_gst)), 
        add_gst = prev_p_chr.bp, add_qtl = prev_g_chr.bp))
}
