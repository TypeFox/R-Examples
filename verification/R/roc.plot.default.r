`roc.plot.default` <-
function (x, pred, thresholds = NULL, binormal = FALSE, legend = FALSE, 
    leg.text = NULL, plot = "emp", CI = FALSE, n.boot = 1000, 
    alpha = 0.05, tck = 0.01, plot.thres = seq(0.1, 0.9, 0.1), show.thres = TRUE, 
    main = "ROC Curve", xlab = "False Alarm Rate", ylab = "Hit Rate", 
    extra = FALSE, ...) 
{    
    
	pred <- as.matrix(pred)
	
	f <- function(x){prod(is.finite(x))}
    id <- is.finite(x) & apply(pred, 1, f )
    x <- x[id]
    pred <- pred[id,]
    
    ##### insert checks
    if(!is.null(thresholds)){
    if(min(diff(thresholds))<0){stop("Thresholds must be listed in ascending order")}
    }
    #####
    if(!is.null(plot) && (plot == "binorm" | plot == "both") & binormal == FALSE) {
        stop("binormal must be TRUE in order to create a binormal plot")
    }
    
    pred   <- as.matrix(pred)
    n.forc <- ncol(pred)
   
    #### if thresholds is null
    if (is.null(thresholds)) {
        thresholds <- sort(as.numeric(unique(pred)))
       if(length(thresholds) > 10000 ){warning("More than 10,000 unique predictions are used as thresholds. 100 equally spaced thresholds are used" )
       	
       	thresholds <- 100}

       if(length(thresholds) > 1000 ){warning("Large amount of unique predictions used as thresholds. Consider specifying thresholds." )}
       

    }
    #### insert check on the number of thresholds?
    if (length(thresholds) == 1) {
        n.thres.bins <- thresholds
        t <- seq(0, 1, 1/n.thres.bins)
        thresholds <- quantile(pred, t)
    }
    
    #####
    
    orig <- as.data.frame(roc.int(x, pred, thres = plot.thres, 
        binormal))
    A.boot <- NULL
    if (CI) {
        A.boot <- numeric()
        D <- cbind(x, pred)
        A <- matrix(NA, ncol = 3)
        for (i in 1:n.boot) {
            nr <- nrow(D)
            ind <- sample(1:nr, size = nr, replace = TRUE)
            sub <- D[ind, ]
            A.boot[i] <- roc.area(D[ind, 1], D[ind, 2])$A
            for (j in 1:length(plot.thres)) {
                A <- rbind(A, roc.int(sub[, 1], sub[, 2], plot.thres[j], 
                  binormal = binormal)[2, 1:3])
            } ### 1:length(plot.thres)
        } ## 1:n.boot
        BOOT <- as.data.frame(A[-1, ])
        xleft <- aggregate(BOOT$F, by = list(BOOT$thres), quantile, 
            alpha)$x
        ybot <- aggregate(BOOT$H, by = list(BOOT$thres), quantile, 
            alpha)$x
        xright <- aggregate(BOOT$F, by = list(BOOT$thres), quantile, 
            1 - alpha)$x
        ytop <- aggregate(BOOT$H, by = list(BOOT$thres), quantile, 
            1 - alpha)$x
        box.corners <- cbind(xleft, ybot, xright, ytop)
        row.names(box.corners) <- plot.thres
    }  ### close CI
    
    DAT <- array(NA, dim = c(length(thresholds) + 1, 5, n.forc))
    VOLS <- matrix(nrow = n.forc, ncol = 3)
    binormal.pltpts <- list()
    for (j in 1:n.forc) {
        DAT[, , j] <- roc.int(x, pred[, j], thresholds, binormal = binormal) 

        if (binormal) {
            dat <- as.data.frame(DAT[, , j])
            names(dat) <- c("thres", "proby", "probn", "zH", 
                "zF")
            dat <- dat[is.finite(dat$zH) & is.finite(dat$zF), 
                ]
            new <- as.data.frame(matrix(qnorm(seq(0.005, 0.995, 
                0.005)), ncol = 1))
            names(new) <- "zF"
            mod <- lm(zH ~ zF, data = dat)
            A <-mod$fitted.values
            B <- predict(lm(zH ~ zF, data = dat), newdata = new)
            binormal.pltpts[[j]] <- data.frame(t = new$zF, x = pnorm(new$zF), 
                y = pnorm(B))
 #           binormal.area <- sum(0.005 * pnorm(B), na.rm = TRUE)
 
 aa <- mod$coefficients[1]
 bb <- mod$coefficients[2]
 
 binormal.area <- pnorm(aa/sqrt(1+bb^2) )   ## from Pepe p. 83
        } ## if binormal
        else {
            binormal.area <- NA
        }
        v <- roc.area(x, pred[, j])
        	
        VOLS[j, 1] <- v$A
        VOLS[j, 2] <- v$p.value
        VOLS[j, 3] <- binormal.area
    }  ## close j in 1:n.forc
    VOLS <- data.frame(paste("Model ", seq(1, n.forc)), VOLS)
    names(VOLS) <- c("Model", "Area", "p.value", "binorm.area")
    r <- structure(list(plot.data = DAT, roc.vol = VOLS, binormal.ptlpts = binormal.pltpts, 
        A.boot = A.boot), class = "roc.data")
    if (!is.null(plot)) {
        par(mar = c(4, 4, 4, 1))
        plot(DAT[, 3, ], DAT[, 2, ], type = "n", xlim = c(0, 
            1), ylim = c(0, 1), main = main, xlab = xlab, ylab = ylab, 
            ...)
        abline(h = seq(0, 1, by = 0.1), v = seq(0, 1, by = 0.1), 
            lty = 3, lwd = 0.5, col = "grey")
        abline(0, 1)
        if (length(thresholds) < 16) {
            L <- "b"
        }
        else {
            L <- "l"
        }
    }

    if(!is.null(plot) && (plot == "emp" | plot == "both")) {
        for (i in 1:n.forc) {
        	
            points(DAT[, 3, i], DAT[, 2, i], col = i, lty = i, 
                type = "l", lwd = 2)
            if (!is.null(plot.thres)) {
                if (show.thres) {
                  a1 <- DAT[, 1, 1]
                  b1 <- plot.thres
                  a <- matrix(a1, ncol = length(b1), nrow = length(a1))
                  X <- abs(scale(a, center = b1, scale = FALSE))
                  X[X > 0.5 * max(diff(b1))] <- NA
                  id <- as.numeric(apply(X, 2, which.min))
                  id2 <- is.finite(id)
                  id <- id[id2]
                  rm(a1, b1, a)
                  points(DAT[id, 3, 1], DAT[id, 2, 1], col = 1, 
                    pch = 19)
                  text(DAT[id, 3, i], DAT[id, 2, i], plot.thres[id2], 
                    pos = 4, offset = 2)
                }
            }
        }
    }
    if (!is.null(plot) && (plot == "binorm" || plot == "both")) {
        for (i in 1:n.forc) {
            dat <- binormal.pltpts[[i]]
            points(dat$x, dat$y, col = 2, lty = i, type = "l", 
                lwd = 2)
        }
    }
    if (!is.null(plot) && extra) {
        if (plot == "both") {
            text(0.6, 0.1, "Black lines are the empirical ROC")
            text(0.6, 0.07, "Red lines and symbols are the bi-normal ROC")
            text(0.6, 0.04, "The area under the binormal curve is in parathesis.")
        }
        if (plot == "emp") {
            text(0.6, 0.1, "Black lines are the empirical ROC")
        }
        if (plot == "binorm") {
            text(0.6, 0.1, "Red lines  are the bi-normal ROC")
        }
    }
    if (!is.null(plot) && CI) {
        for (i in 1:nrow(box.corners)) {
            lines(box.corners[i, c(1, 3)], rep(orig$H[i + 1], 
                2), lwd = 1)
            lines(rep(box.corners[i, 1], 2), c(orig$H[i + 1] - 
                tck, orig$H[i + 1] + tck), lwd = 1)
            lines(rep(box.corners[i, 3], 2), c(orig$H[i + 1] - 
                tck, orig$H[i + 1] + tck), lwd = 1)
            lines(rep(orig$F[i + 1], 2), box.corners[i, c(2, 
                4)], lwd = 1)
            lines(c(orig$F[i + 1] - tck, orig$F[i + 1] + tck), 
                rep(box.corners[i, 2], 2), lwd = 1)
            lines(c(orig$F[i + 1] - tck, orig$F[i + 1] + tck), 
                rep(box.corners[i, 4], 2), lwd = 1)
        }
    }
    if (!is.null(plot) && legend) {
        if (is.null(leg.text)) {
            leg.text <- paste("Model ", LETTERS[seq(1, n.forc)])
        }
        if (plot == "emp") {
            leg.text <- paste(leg.text, "  ", formatC(VOLS$Area, 
                digits = 3))
        }
        if (plot == "binorm") {
            leg.text <- paste(leg.text, "  ", formatC(VOLS$binorm.area, 
                digits = 3))
        }
        if (plot == "both") {
            leg.text <- paste(leg.text, "  ", formatC(VOLS$Area, 
                digits = 3), " (", formatC(VOLS$binorm.area, 
                digits = 3), ")")
        }
        legend(list(x = 0.6, y = 0.4), legend = leg.text, bg = "white", 
            cex = 0.6, lty = seq(1, n.forc), col = c("black", 
                "red", "blue"), merge = TRUE)
    }

    invisible(r)
}
