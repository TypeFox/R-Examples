enirg.plot <-
function(enirg.results, mar.col = "grey", spe.col = "black", method = "extended", 
    plot.egvs = TRUE, asp = FALSE, title = NULL) {
    if (class(enirg.results) != "enirg") 
        stop("The function plot.enirg needs an object of class 'enirg'!")
    if (method == "extended") {
        wrapping <- function(tmp) as.numeric(strsplit(tmp, " ")[[1]])
        x <- execGRASS("r.stats", input = paste(enirg.results$species, "_li_Mar", sep = ""), flags = c("1"), intern = TRUE, 
            legacyExec = TRUE)
        x <- as.numeric(x[which(x != "*")])
        y <- execGRASS("r.stats", input = paste(enirg.results$species, "_li_Spec1", sep = ""), flags = c("1"), intern = TRUE, 
            legacyExec = TRUE)
        y <- as.numeric(y[which(y != "*")])
        range.y <- c(min(y) - abs(min(y)/4), max(y) + abs(max(y)/4))
        range.x <- c(min(x) - abs(min(x)/4), max(x) + abs(max(x)/4))
        if(asp) {
            range.y <- c(min(c(range.y, range.x)), max(c(range.y, range.x)))
            range.x <- range.y
        }
        plot(0, 0, col = mar.col, pch = 19, xlab = "Marginality", ylab = "Specialization 1", 
            type = "n", xlim = range.x, ylim = range.y, main = title)
        hpts <- chull(x, y)
        hpts <- c(hpts, hpts[1])
        polygon(x[hpts], y[hpts], col = mar.col, border = mar.col)
    }
    if (method == "simplified") {
        give.stat <- function(map, stat) {
            result <- execGRASS("r.univar", map = map, flags = "g", intern = TRUE)
            result <- agrep(result, pattern = stat, max.distance = list(all = 0), 
                value = T)
            result <- as.numeric(strsplit(result, "=")[[1]][2])
            return(result)
        }
        mar.range <- c(give.stat(paste(enirg.results$species, "_li_Mar", sep = ""), "min"),
                       give.stat(paste(enirg.results$species, "_li_Mar", sep = ""), "max"))
        spec.range <- c(give.stat(paste(enirg.results$species, "_li_Spec1", sep = ""), "min"),
                        give.stat(paste(enirg.results$species, "_li_Spec1", sep = ""), "max"))
        range.x <- c(min(mar.range) - abs(min(mar.range)/4), max(mar.range) + abs(max(mar.range)/4))
        range.y <- c(min(spec.range) - abs(min(spec.range)/4), max(spec.range) + abs(max(spec.range)/4))
        if(asp) {
            range.y <- c(min(c(range.y, range.x)), max(c(range.y, range.x)))
            range.x <- range.y
        }
        plot(0, 0, xlab = "Marginality", ylab = "Specialization 1", 
            type = "n", xlim = range.x, ylim = range.y, main = title)
        polygon(c(0, mar.range[1], 0, mar.range[2]), c(spec.range[1], 0, spec.range[2], 
            0), col = mar.col)
    }
    data.pres <- enirg.results$obs.li
    points(data.pres[, "li_Mar"], data.pres[, "li_Spec1"], col = spe.col, pch = 19, 
        cex = data.pres[, "presences"]/10)
    abline(v = 0, lty = 2)
    abline(h = 0, lty = 2)
    pmar <- t(enirg.results$mar * enirg.results$cw) %*% as.matrix(enirg.results$co[, 
        1:2])
    dfarr <- enirg.results$co[, 1:2]
    born <- par("usr")
    k1 <- min(dfarr[, 1])/born[1]
    k2 <- max(dfarr[, 1])/born[2]
    k3 <- min(dfarr[, 2])/born[3]
    k4 <- max(dfarr[, 2])/born[4]
    k <- c(k1, k2, k3, k4)
    dfarr <- 0.75 * dfarr/max(k)
    symbols(pmar, circles = 1, fg = "black", bg = "white", inches = 0.03 * 2, add = TRUE)
    if (plot.egvs) {
        if (!is.null(enirg.results$qt.egvs)) {
            qt <- dfarr[which(rownames(dfarr) %in% enirg.results$qt.egvs), ]
            s.arrow(qt, clabel = 1, addaxes = FALSE, add.plot = TRUE)
            arrows(x0 = rep(0, nrow(qt)), y0 = rep(0, nrow(qt)), x1 = qt[, 1], y1 = qt[, 
                2], col = "blue", length = 0.1, angle = 15, lwd = 2)
        }
        if (!is.null(enirg.results$ql.egvs)) {
            ql <- dfarr[which(rownames(dfarr) %in% enirg.results$ql.egvs), ]
            text(ql, labels = rownames(ql), cex = 1.5)
        }
    }
}
