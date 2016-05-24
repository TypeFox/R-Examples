catseye <-
function(DV, grp = NULL, plotFUN = mean, errFUN = c("ci", "se", "sd"),
                conf = .95, xpoints = NULL, grp.names = NULL, tick = FALSE,
                ylim = NULL, col = "gray", len = 0, ...)
{
    se <- function(x) {
        x <- na.omit(x)
        res <- sqrt(var(x)/length(x))
        res
    }
    ci <- function(x) {
        x <- na.omit(x)
        alpha <- 1 - (1 - conf)/2
        res <- qt(alpha, length(x) - 2, lower.tail = T) * se(x)
        res
    }
    if (!is.null(errFUN)) {
        errFUN <- match.arg(errFUN)
    }
    if (is.null(grp)) {
        DV <- na.omit(DV)
        res <- do.call(plotFUN, list(DV))
        if (is.null(errFUN)) {
            plot(res, pch = 19, xaxt = "n", ylim = ylim, ...)
        }
        if (!is.null(errFUN)) {
            e <- do.call(errFUN, list(DV))
            e <- ifelse(res <= 0, -e, e)
            Q <- qnorm(seq(.001, .999, by = .001), se(DV))
            Qscale <- Q * se(DV) + res
            if (is.null(ylim)) {
                lims <- c(min(Qscale), max(Qscale))
            }
            else (lims <- ylim)
            plot(res, ylim = lims, xaxt = "n", ...)
            Qscale2 <- sort(c(Qscale, res - e, res + e))
            poly1 <- which(Qscale2 == res - e)
            poly2 <- which(Qscale2 == res + e)
            polygon(x = c(c(1 - dnorm(Q))[poly1:poly2], c(1 + dnorm(Q))[poly2:poly1]), y = c(Qscale2[poly1:poly2], Qscale2[poly2:poly1]), border = NA, col = col)
            points(x = 1, y = res, pch = 19)
            arrows(1, res + e, 1, res - e, angle = 90, code = 3, length = len)
            lines(x = 1 - dnorm(Q), y = Qscale)
            lines(x = 1 + dnorm(Q), y = Qscale)
        }
        if (is.null(grp.names)) {
            grp.names <- ""
        }
        axis(1, at = 1, labels = grp.names, tick = tick)
    }
    if (!is.null(grp) & !is.list(grp)) {
        dat <- data.frame(DV, grp)
        dat <- dat[complete.cases(dat), ]
        res <- tapply(dat[, 1], dat[, 2], plotFUN)
        if (is.null(xpoints)) {
          places <- 1:length(res)
        }
        else places <- xpoints
        if (is.null(errFUN)) {
            plot(res ~ places, pch = 19, xaxt = "n", xlim = c(0.4, 
                0.4 + places[length(places)]), ylim = ylim, ...)
        }
        if (!is.null(errFUN)) {
            e <- tapply(dat[, 1], dat[, 2], errFUN)
            e <- ifelse(res <= 0, -e, e)
            SEs <- tapply(dat[, 1], dat[, 2], se)
            Q <- matrix(unlist(tapply(dat[, 1], dat[, 2], function(x) qnorm(seq(.001, .999, by = .001), se(x)))), nrow = 999)
            Qscale <- Q * matrix(SEs, nrow = 999, ncol = ncol(Q), byrow = TRUE) + matrix(res, nrow = 999, ncol = ncol(Q), byrow = TRUE)
            if (is.null(ylim)) {
                lims <- c(min(Qscale), max(Qscale))
            }
            else (lims <- ylim)
            plot(res ~ places, pch = 19, xaxt = "n", xlim = c(0.4,
                0.4 + places[length(places)]), ylim = lims, ...)
            Qscale2 <- apply(rbind(Qscale, res - e, res + e), 2, sort) 
            poly1 <- rep(NA, ncol(Qscale2))
            poly2 <- rep(NA, ncol(Qscale2))
            for(i in 1:ncol(Qscale2)) {
              poly1[i] <- which(Qscale2[, i] == res[i] - e[i])
              poly2[i] <- which(Qscale2[, i] == res[i] + e[i])
              polygon(x = c(c(places[i] - dnorm(Q[, i]))[poly1[i]:poly2[i]], c(places[i] + dnorm(Q[, i]))[poly2[i]:poly1[i]]), 
                        y = c(Qscale2[poly1[i]:poly2[i], i], Qscale2[poly2[i]:poly1[i], i]), border = NA, col = col)
              lines(x = places[i] - dnorm(Q[, i]), y = Qscale[, i])
              lines(x = places[i] + dnorm(Q[, i]), y = Qscale[, i])
            }
            points(x = places, y = res, pch = 19)
            arrows(places, res + e, places, res - e, angle = 90, code = 3, length = len)
        }
        if (is.null(grp.names)) {
            grp.names <- 1:length(places)
        }
        axis(1, at = places, labels = grp.names, tick = tick)
    }
    if (is.list(grp)) {
        if (length(unique(unlist(lapply(grp, length)))) != 1) {
            stop("Grouping variables must be the same length.")
        }
        if (length(DV) != lapply(grp, length)[[1]]) {
            stop("DV must be the same length as the grouping variables.")
        }
        dat <- data.frame(DV, matrix(unlist(grp), nrow = length(DV), byrow = FALSE))
        if (sum(is.na(dat)) > 0 ) {
          stop("Please remove missing values in DV and IV first.")
        }
        res <- as.vector(tapply(DV, grp, plotFUN))
        if (is.null(xpoints)) {
          places <- 1:length(res)
        }
        else places <- xpoints
        if (is.null(errFUN)) {
            plot(res ~ places, pch = 19, xaxt = "n", xlim = c(0.4,
                0.4 + places[length(places)]), ylim = ylim, ...)
        }
        if (!is.null(errFUN)) {
            e <- as.vector(tapply(DV, grp, errFUN))
            e <- ifelse(res <= 0, -e, e)
            SEs <- as.vector(tapply(DV, grp, se))
            Q <- matrix(qnorm(seq(.001, .999, by = .001), SEs), nrow = 999, ncol = length(SEs))
            Qscale <- Q * matrix(SEs, nrow = 999, ncol = ncol(Q), byrow = TRUE) + matrix(res, nrow = 999, ncol = ncol(Q), byrow = TRUE)
            if (is.null(ylim)) {
                  lims <- c(min(Qscale), max(Qscale))
            }
            else (lims <- ylim)
            plot(res ~ places, pch = 19, xaxt = "n", xlim = c(0.4,
                0.4 + places[length(places)]), ylim = lims, ...)
            Qscale2 <- apply(rbind(Qscale, res - e, res + e), 2, sort) 
            poly1 <- rep(NA, ncol(Qscale2))
            poly2 <- rep(NA, ncol(Qscale2))
            for(i in 1:ncol(Qscale2)) {
              poly1[i] <- which(Qscale2[, i] == res[i] - e[i])
              poly2[i] <- which(Qscale2[, i] == res[i] + e[i])
              polygon(x = c(c(places[i] - dnorm(Q[, i]))[poly1[i]:poly2[i]], c(places[i] + dnorm(Q[, i]))[poly2[i]:poly1[i]]), 
                        y = c(Qscale2[poly1[i]:poly2[i], i], Qscale2[poly2[i]:poly1[i], i]), border = NA, col = col)
              lines(x = places[i] - dnorm(Q[, i]), y = Qscale[, i])
              lines(x = places[i] + dnorm(Q[, i]), y = Qscale[, i])
            }
            points(x = places, y = res, pch = 19)
            arrows(places, res + e, places, res - e, angle = 90, code = 3, length = len)
        }
        if (is.null(grp.names)) {
            grp.names <- 1:length(places)
        }
        axis(1, at = places, labels = grp.names, tick = tick)
    }
}
