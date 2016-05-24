plottol <- function (tol.out, x, y = NULL, y.hat = NULL, side = c("two", 
    "upper", "lower"), plot.type = c("control", "hist", "both"), 
    x.lab = NULL, y.lab = NULL, z.lab = NULL, ...) 
{
    if (class(tol.out) == "list" & is.null(names(tol.out)) == 
        FALSE & class(x) == "data.frame") {
        temp <- NULL
        for (i in 1:length(tol.out)) temp = c(temp, unlist(tol.out[[i]][, 
            4:5]))
        resp <- comment(tol.out)[1]
        a <- 1 - as.numeric(comment(tol.out)[2])
        P <- as.numeric(comment(tol.out)[3])
        temp.ind <- which(names(x) == resp)
        fact = length(tol.out)
        if (colnames(tol.out[[1]])[4] == "1-sided.lower") {
		par(mfrow = c(round(sqrt(fact)), ceiling(sqrt(fact))))
            for (i in 1:fact) {
                x.axis <- 1:nrow(tol.out[[i]])
                ll <- length(x.axis)
                plot(x.axis, tol.out[[i]][, 1], xlim = range(x.axis) + 
                  c(-0.5, 0.5), ylim = range(c(temp, x[, temp.ind])), 
                  xlab = names(tol.out)[[i]], ylab = resp, xaxt = "n", 
                  pch = 15, cex = 1.5, col = 0)
                points(c(x.axis), tol.out[[i]][, 1], pch = 19, 
                  cex = 1.5)
                axis(side = 1, at = x.axis, labels = rownames(tol.out[[i]]))
                title(main = paste(a * 100, "% / ", P * 100, 
                  "% Tolerance Intervals for ", names(tol.out)[[i]], 
                  sep = ""))
                title(sub = "Red lines are 1-sided upper intervals and green lines are 1-sided lower intervals.", 
                  cex.sub = 0.8)
                for (j in 1:nrow(tol.out[[i]])) {
                  segments(x.axis[j] - 0.1, -1e+100, x.axis[j] - 
                    0.1, tol.out[[i]][j, 5], col = 2)
                  segments(x.axis[j] + 0.1, tol.out[[i]][j, 4], 
                    x.axis[j] + 0.1, 1e+100, col = 3)
                }
            }
        }
        else {
		par(mfrow = c(round(sqrt(fact)), ceiling(sqrt(fact))))
            for (i in 1:length(tol.out)) {
                x.axis <- 1:nrow(tol.out[[i]])
                plot(x.axis, tol.out[[i]][, 1], xlim = range(x.axis) + 
                  c(-0.5, 0.5), ylim = range(c(temp, x[, temp.ind])), 
                  xlab = names(tol.out)[[i]], ylab = resp, xaxt = "n", 
                  pch = 19, cex = 1.5)
                axis(side = 1, at = x.axis, labels = rownames(tol.out[[i]]))
                title(main = paste(a * 100, "% / ", P * 100, 
                  "% Tolerance Intervals for ", names(tol.out)[[i]], 
                  sep = ""))
                for (j in 1:nrow(tol.out[[i]])) {
                  segments(x.axis[j], tol.out[[i]][j, 4], x.axis[j], 
                    tol.out[[i]][j, 5], col = 2)
                }
            }
        }
    }
    else {
        if (length(x) == 1) {
            stop(paste("There are no plots produced for discrete distribution tolerance intervals.", 
                "\n"))
        }
        if (class(tol.out) == "list") {
            y.lim = range(sapply(1:length(tol.out), function(i) tol.out[[i]][, 
                c(ncol(tol.out[[i]]) - 1, ncol(tol.out[[i]]))]))
            temp.tol <- tol.out
            tol.out <- tol.out[[1]]
        }
        else {
            temp.tol <- NULL
            y.lim <- range(tol.out[, c(ncol(tol.out) - 1, ncol(tol.out))])
        }
        if (is.null(x.lab)) 
            x.lab <- "X"
        if (is.null(y.lab)) 
            y.lab <- "Y"
        if (is.null(z.lab)) 
            z.lab <- "Z"
        plot.type <- match.arg(plot.type)
        if (plot.type == "both") {
            par(mfrow = c(1, 2))
        }
        if (is.matrix(x) & is.null(y)) {
            P <- as.numeric(rownames(tol.out))[1]
            a <- 1 - as.numeric(colnames(tol.out))[1]
        }
        else {
            side <- match.arg(side)
            a <- 1 - tol.out[1, 1]
            P <- tol.out[1, 2]
            out <- tol.out
            n.c <- ncol(tol.out)
            n.r <- nrow(tol.out)
            if (max(tol.out[, n.c]) == Inf) 
                tol.out[, n.c] <- max(x)
            if (min(tol.out[, (n.c - 1)]) == -Inf) 
                tol.out[, n.c] <- min(x)
        }
        if (is.null(y)) {
            if (is.matrix(x)) {
                if (ncol(x) == 2) {
                  mu <- apply(x, 2, mean)
                  sigma <- cov(x)
                  es <- eigen(sigma)
                  e1 <- es$vec %*% diag(sqrt(es$val))
                  theta <- seq(0, 2 * pi, len = 1000)
                  r1 <- sqrt(tol.out[1, 1])
                  v1 <- cbind(r1 * cos(theta), r1 * sin(theta))
                  pts <- t(mu - (e1 %*% t(v1)))
                  plot(pts, col = 0, xlab = "", ylab = "")
                  lines(pts, col = 2)
                  invisible(pts)
                  points(x, pch = 19)
                  title(main = paste(a * 100, "% / ", P * 100, 
                    "% Tolerance Region", sep = ""), xlab = x.lab, 
                    ylab = y.lab, ...)
                }
                else if (ncol(x) == 3) {
                  open3d()
                  plot3d(x, size = 3, box = FALSE, xlab = x.lab, 
                    ylab = y.lab, zlab = z.lab)
                  Sigma <- cov(x)
                  Mean <- apply(x, 2, mean)
                  plot3d(ellipse3d(Sigma, centre = Mean, t = sqrt(tol.out)), 
                    col = "red", add = TRUE, alpha = 0.3)
                  title3d(main = paste(a * 100, "% / ", P * 100, 
                    "% Tolerance Region", sep = ""), ...)
                }
            }
            else {
                if ((plot.type == "both") | (plot.type == "control")) {
                  if (side == "upper") {
                    plot(x, ylim = c(min(x), max(x, tol.out[, 
                      n.c])), type = "l", main = "", ylab = x.lab, 
                      ...)
                    points(x, pch = 19)
                    title(main = paste(a * 100, "% / ", P * 100, 
                      "% Upper Tolerance Limit", sep = ""))
                    abline(h = out[, n.c], lty = "dashed", col = 2:(n.r + 
                      1))
                  }
                  else if (side == "lower") {
                    plot(x, ylim = c(min(x, tol.out[, (n.c - 
                      1)]), max(x)), type = "l", main = "", ylab = x.lab, 
                      ...)
                    points(x, pch = 19)
                    title(main = paste(a * 100, "% / ", P * 100, 
                      "% Lower Tolerance Limit", sep = ""))
                    abline(h = out[, (n.c - 1)], lty = "dashed", 
                      col = 2:(n.r + 1))
                  }
                  else {
                    if (colnames(tol.out)[(n.c - 1)] == "1-sided.lower") 
                      print("NOTE: The plot reflects two 1-sided tolerance intervals and NOT a 2-sided tolerance interval!")
                    plot(x, ylim = c(min(x, tol.out[, (n.c - 
                      1)]), max(x, tol.out[, n.c])), type = "l", 
                      main = "", ylab = x.lab, ...)
                    points(x, pch = 19)
                    title(main = paste(a * 100, "% / ", P * 100, 
                      "% Tolerance Limits", sep = ""))
                    abline(h = out[, (n.c - 1)], lty = "dashed", 
                      col = 2:(n.r + 1))
                    abline(h = out[, n.c], lty = "dashed", col = 2:(n.r + 
                      1))
                  }
                }
                if ((plot.type == "both") | (plot.type == "hist")) {
                  if (side == "upper") {
                    hist(x, prob = TRUE, main = "", xlab = x.lab, 
                      xlim = c(min(x), max(x, tol.out[, n.c])))
                    title(main = paste(a * 100, "% / ", P * 100, 
                      "% Upper Tolerance Limit", sep = ""))
                    abline(v = out[, n.c], lty = "dashed", col = 2:(n.r + 
                      1))
                  }
                  else if (side == "lower") {
                    hist(x, prob = TRUE, main = "", xlab = x.lab, 
                      xlim = c(min(x, tol.out[, (n.c - 1)]), 
                        max(x)))
                    title(main = paste(a * 100, "% / ", P * 100, 
                      "% Lower Tolerance Limit", sep = ""))
                    abline(v = out[, (n.c - 1)], lty = "dashed", 
                      col = 2:(n.r + 1))
                  }
                  else {
                    if (colnames(tol.out)[(n.c - 1)] == "1-sided.lower") 
                      print("NOTE: The histogram reflects two 1-sided tolerance intervals and NOT a 2-sided tolerance interval!")
                    hist(x, prob = TRUE, main = "", xlab = x.lab, 
                      xlim = c(min(x, tol.out[, (n.c - 1)]), 
                        max(x, tol.out[, n.c])))
                    title(main = paste(a * 100, "% / ", P * 100, 
                      "% Tolerance Limits", sep = ""))
                    abline(v = out[, (n.c - 1)], lty = "dashed", 
                      col = 2:(n.r + 1))
                    abline(v = out[, n.c], lty = "dashed", col = 2:(n.r + 
                      1))
                  }
                }
            }
        }
        else {
            par(mfrow = c(1, 1))
            if (colnames(out)[3] == "y.hat") {
                out1 <- cbind(x, y)
                out1 <- out1[order(out1[, 2]), ]
                out1 <- cbind(out1[, 1], out[1:length(x), ])
                out1 <- out1[order(out1[, 1]), ]
                plot(out1[, 1], out1[, 5], xlab = x.lab, ylab = y.lab, 
                  ylim = y.lim, pch = 19, ...)
                lines(out1[, 1], out1[, 4])
                if (side == "upper") {
                  title(main = paste(a * 100, "% / ", P * 100, 
                    "% Upper Tolerance Limit", sep = ""))
                  lines(out1[, 1], out1[, 7], col = 2)
                }
                else if (side == "lower") {
                  title(main = paste(a * 100, "% / ", P * 100, 
                    "% Lower Tolerance Limit", sep = ""))
                  lines(out1[, 1], out1[, 6], col = 2)
                }
                else if (side == "two") {
                  if (colnames(out1)[6] == "1-sided.lower") 
                    print("NOTE: The plot reflects two 1-sided tolerance intervals and NOT a 2-sided tolerance interval!")
                  title(main = paste(a * 100, "% / ", P * 100, 
                    "% Tolerance Limits", sep = ""))
                  lines(out1[, 1], out1[, 6], col = 2)
                  lines(out1[, 1], out1[, 7], col = 2)
                }
            }
            if (is.null(temp.tol) == "FALSE" | is.null(y.hat) == 
                "FALSE") {
                plot(x, y, xlab = x.lab, ylab = y.lab, ylim = y.lim, pch = 19,
                  ...)
                if (class(temp.tol) != "list") 
                  temp.tol = list(tol.out)
                out1 <- temp.tol
                len <- length(out1)
                xy.out <- lapply(1:len, function(i) out1[[i]][, 
                  3:7])
                xy.out <- lapply(1:len, function(i) xy.out[[i]][order(xy.out[[i]][, 
                  1]), ])
                lines(xy.out[[1]][, 1], xy.out[[1]][, 3])
                if (side == "upper") {
                  title(main = paste(a * 100, "% / ", P * 100, 
                    "% Upper Tolerance Limit", sep = ""))
                  for (i in 1:len) {
                    lines(xy.out[[i]][, 1], xy.out[[i]][, 5], 
                      col = i + 1)
                  }
                }
                else if (side == "lower") {
                  title(main = paste(a * 100, "% / ", P * 100, 
                    "% Lower Tolerance Limit", sep = ""))
                  for (i in 1:len) {
                    lines(xy.out[[i]][, 1], xy.out[[i]][, 4], 
                      col = i + 1)
                  }
                }
                else if (side == "two") {
                  if (colnames(out1[[1]])[5] == "1-sided.lower") 
                    print("NOTE: The plot reflects two 1-sided tolerance intervals and NOT a 2-sided tolerance interval!")
                  title(main = paste(a * 100, "% / ", P * 100, 
                    "% Tolerance Limits", sep = ""))
                  for (i in 1:len) {
                    lines(xy.out[[i]][, 1], xy.out[[i]][, 4], 
                      col = i + 1)
                    lines(xy.out[[i]][, 1], xy.out[[i]][, 5], 
                      col = i + 1)
                  }
                }
            }
            if (colnames(out)[3] == "y") {
                y.lim <- range(out[, 5:6])
                if (class(x) == "numeric") 
                  x <- matrix(x, ncol = 1)
                if (sum(x[, 1] == 1) != nrow(x)) 
                  print("NOTE: A regression through the origin is fitted!")
                if (sum(x[, 1]) == nrow(x)) {
                  XXX <- x[, 2]
                  plot(x[, 2], y, xlab = x.lab, ylab = y.lab, 
                    ylim = y.lim, pch = 19, ...)
                  reg.out <- lm(y ~ x - 1)
                  temp.x <- seq(min(x[, 2]), max(x[, 2]), length.out = 1000)
                }
                else {
                  plot(x[, 1], y, xlab = x.lab, ylab = y.lab, 
                    ylim = y.lim, pch = 19, ...)
                  reg.out <- lm(y ~ x - 1)
                  temp.x <- seq(min(x[, 1]), max(x[, 1]), length.out = 1000)
                  XXX <- x[, 1]
                }
                b <- reg.out$coef
                out.2 <- out[order(out[, 4]), ]
                out.temp <- apply(is.na(out.2), 1, sum)
                out.2 <- out.2[out.temp == 0, ]
                if (sum(x[, 1]) == nrow(x)) {
                  temp.x = cbind(1, sapply(1:(ncol(x) - 1), function(i) temp.x^i))
                  poly.x = cbind(temp.x[, 2], apply(t(b * t(temp.x)), 
                    1, sum))
                  n.x <- ncol(poly.x)
                }
                else {
                  temp.x <- as.matrix(sapply(1:ncol(x), function(i) temp.x^i), 
                    ncol = ncol(x))
                  poly.x <- cbind(temp.x[, 1], apply(t(b * t(temp.x)), 
                    1, sum))
                  n.x <- ncol(poly.x)
                }
                lines(poly.x)
                x.temp <- cbind(XXX, fitted(reg.out))
                x.temp <- x.temp[order(x.temp[, 2]), ]
                temp <- cbind(x.temp[, 1], out.2[, 4:6])
                temp <- temp[order(temp[, 1]), ]
                if (side == "upper") {
                  title(main = paste(a * 100, "% / ", P * 100, 
                    "% Upper Tolerance Limit", sep = ""))
                  lines(temp[, 1], temp[, 4], col = 2)
                }
                else if (side == "lower") {
                  title(main = paste(a * 100, "% / ", P * 100, 
                    "% Lower Tolerance Limit", sep = ""))
                  lines(temp[, 1], temp[, 3], col = 2)
                }
                else if (side == "two") {
                  if (colnames(tol.out)[5] == "1-sided.lower") 
                    print("NOTE: The plot reflects two 1-sided tolerance intervals and NOT a 2-sided tolerance interval!")
                  title(main = paste(a * 100, "% / ", P * 100, 
                    "% Tolerance Limits", sep = ""))
                  for (i in 1:nrow(out)) {
                    lines(temp[, 1], temp[, 3], col = 2)
                    lines(temp[, 1], temp[, 4], col = 2)
                  }
                }
            }
        }
    }
}
