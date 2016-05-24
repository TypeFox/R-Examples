`plot.dlmap` <- 
function(x, chr, max.dist, qcol="light blue", mcol="red", pcol="purple", marker.names=FALSE, ...)
{
    require(wgaim)

    dots <- list(...)
    if (missing(x)) 
        stop("x is a required argument")

    if (attr(x$input, "type")=="other")
    {
	cat("Cannot plot the detected QTL for an association analysis")
	return(invisible())
    }

    parentData <- x$input$genCross
    wchr <- as.character(x$Summary[,1])

    if (missing(chr)) chr <- 1:nchr(parentData) 

    if (!inherits(parentData, "cross")) 
        stop("input$genCross is not of class \"cross\"")
    if (!length(wchr)) {
        warning("There are no significant QTL. Plotting map only...")
        link.map.cross(parentData, chr, max.dist, ...)
        return(invisible())
    } else chr <- unique(wchr)

    if (is.numeric(chr))
	chr <- names(pull.map(parentData))[chr]

    lmap <- link.map.cross(parentData, chr, max.dist, ...)

    map <- lmap$map

    qtlm <- matrix(nrow=nrow(x$Summary), ncol=4)
    qtlm[,1] <- as.vector(x$Summary[,3])
    qtlm[,3] <- as.vector(x$Summary[,4])
    qtlm[,2] <- round(find.markerpos(parentData, qtlm[,1])[,2],2)
    qtlm[,4] <- round(find.markerpos(parentData, qtlm[,3])[,2],2)

    qtlpos <- as.numeric(as.character(x$Summary[,2]))

    # note, this may not work for lme
    trait <- rep(as.character(x$final.model$call$fixed[[2]]), length(wchr))

    qtlm <- cbind(qtlm, as.numeric(factor(trait, levels = unique(trait))))
    if (!missing(chr)) {
        if (any(is.na(wh <- pmatch(wchr, chr, duplicates.ok = TRUE)))) {
            warning("Some QTL exist outside chromosome(s) subset, Omitting QTL....")
            qtlm <- qtlm[!is.na(wh), ]
            wchr <- wchr[!is.na(wh)]
        }
    }
    if (!missing(max.dist)) {
        rml <- as.numeric(qtlm[, 4]) > max.dist
        if (any(rml)) {
            warning("Some QTL regions outside maximum distance specified. Omitting QTL's....")
            qtlm <- matrix(qtlm[!rml, ], nrow = length(rml[!rml]), 
                byrow = FALSE)
            wchr <- wchr[!rml]
        }
    }
    n.chr <- length(map)
    maxlen <- max(unlist(lapply(map, max)))
    chrpos <- lmap$chrpos
    mt <- lmap$mt
    if (!is.na(cind <- pmatch("col", names(dots)))) 
        dots <- dots[-cind]
    if (is.null(dim(qtlm))) 
        qtlm <- matrix(qtlm, nrow = 1, byrow = FALSE)
    qtld <- cbind.data.frame(qtlm)[, 1:4]
    nodup <- !duplicated(do.call("paste", qtld))
    qtls <- qtld[nodup, ]
    whd <- pmatch(do.call("paste", qtld), do.call("paste", qtls), 
        duplicates.ok = TRUE)
    dlis <- split(qtlm[, 5], whd)
    qtlm <- as.matrix(qtls)
    wchr <- wchr[nodup]
    for (i in 1:n.chr) {
        if (as.logical(length(ind <- grep(names(map)[i], 
            wchr)))) {
            for (j in ind) {
                if (marker.names) {
                  wh <- mt[[i]][pmatch(c(as.character(qtlm[j, 
                    1]), as.character(qtlm[j, 3])), names(map[[i]]))]
                  alis <- list(x = chrpos[i] + 0.5, y = wh, labels = names(wh), 
                    adj = c(0, 0.5), col = mcol, cex=.6)
                  do.call("text", c(alis, dots))
                }
                yv <- c(as.numeric(qtlm[j, 2]), as.numeric(qtlm[j,4])) 
                yv2 <- c(qtlpos[j]-1, qtlpos[j]+1)
                yv <- c(yv, rev(yv))
                yv2 <- c(yv2, rev(yv2))
                if (length(dlis[[j]]) > 1) {
                  int <- seq(chrpos[i] - 0.2, chrpos[i] + 0.2, 
                    length = length(dlis[[j]]) + 1)
                  int2 <- seq(chrpos[i] - 0.1, chrpos[i] + 0.1, 
                    length = length(dlis[[j]]) + 1)
                  qcols <- qcol[as.numeric(dlis[[j]])]
                  for (k in 1:length(dlis[[j]])) {
                    xv <- c(rep(int[k], 2), rep(int[k + 1], 2))
                    xv2 <- c(rep(int2[k], 2), rep(int2[k + 1], 2))
                    polygon(xv, y = yv, border = NA, col = qcol[k])
                    polygon(xv2, y = yv2, border = NA, col = pcol)
                  }
                }
                else {
                  xv <- c(rep(chrpos[i] - 0.2, 2), rep(chrpos[i] + 0.2, 2)) 
                  xv2 <- c(rep(chrpos[i] - 0.1, 2), rep(chrpos[i] + 
                    0.1, 2))
                  polygon(xv, y = yv, border = NA, col = qcol)
                  polygon(xv2, y = yv2, border = NA, col = pcol)
                }
            }
        }
        segments(chrpos[i] - 0.2, map[[i]], chrpos[i] + 0.2, 
            map[[i]])
    }
    if (is.na(pmatch("main", names(dots)))) 
        title("Genetic Map with QTL")
}
 
