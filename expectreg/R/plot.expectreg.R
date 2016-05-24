plot.expectreg <-
function (x, rug = TRUE, xlab = NULL, ylab = NULL, ylim = NULL, 
    legend = TRUE, ci = FALSE, ...) 
{
    ask = prod(par("mfcol")) < sum(unlist(x$effects) != "parametric") && 
        dev.interactive()
    if (ask) {
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
    }
    if (inherits(x, "boost") && ci) 
        warning("no confidence intervals are calculated while boosting.")
    yy = x$response
    cov = x$covariates
    Z = x$values
    coefficients = x$coefficients
    formula = x$formula
    intercept = x$intercepts
    m = length(yy)
    types = x$effects
    helper = x$helper
    pp = x$asymmetries
    np <- length(pp)
    if (is.null(ylab)) 
        ylab = attr(yy, "name")
    ylim2 = ylim
    if (identical(pp, seq(0.01, 0.99, by = 0.01))) {
        pp.plot <- c(1, 2, 5, 10, 20, 50, 80, 90, 95, 98, 99)
        row.grid = 3
        col.grid = 4
    }
    else if (identical(pp, c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 
        0.8, 0.9, 0.95, 0.98, 0.99))) {
        pp.plot <- 1:length(pp)
        row.grid = 3
        col.grid = 4
    }
    else {
        if (np > 15) 
            pp.plot = seq(1, np, length = 15)
        else pp.plot <- 1:length(pp)
        row.grid = floor(sqrt(length(pp)))
        col.grid = ceiling(sqrt(length(pp)))
        if (length(pp) > row.grid * col.grid) 
            row.grid = row.grid + 1
    }
    np.plot <- length(pp.plot)
    if (is.null(xlab)) 
        xlab = names(cov)
    else if (length(xlab) < length(types)) 
        xlab = rep(xlab[1], length(types))
    nb = vector()
    for (k in 1:length(types)) {
        if (!inherits(x, "boost")) {
            nb[k] = ncol(x$bases[[k]]$B)
            partbasis = (sum(nb[0:(k - 1)]) + 1):(sum(nb[0:k]))
        }
        if (types[[k]] == "pspline") {
            if (inherits(x, "boost")) {
                ZZZ = Z[[k]][order(cov[[k]])[seq(1, m, length = min(m, 
                  100))], pp.plot, drop = F]
            }
            else {
                ndat = data.frame(seq(min(cov[[k]]), max(cov[[k]]), 
                  length = 100))
                names(ndat) = names(cov)[k]
                Bpred = predict(x$bases[[k]], ndat)
                ZZZ = Bpred %*% coefficients[[k]]
            }
            lower = NA
            upper = NA
            if (ci) {
                lower = matrix(NA, nrow = 100, ncol = np)
                upper = matrix(NA, nrow = 100, ncol = np)
                for (i in 1:np) for (nn in 1:nrow(Bpred)) {
                  deviation = qnorm(0.975) * sqrt(t(c(1, Bpred[nn, 
                    ])) %*% x$covmat[[i]][c(1, partbasis + 1), 
                    c(1, partbasis + 1)] %*% c(1, Bpred[nn, ]))
                  lower[nn, i] = ZZZ[nn, i] - deviation
                  upper[nn, i] = ZZZ[nn, i] + deviation
                }
            }
            for (i in 1:np) {
                ZZZ[, i] = ZZZ[, i] + intercept[i]
                if (ci) {
                  lower[, i] = lower[, i] + intercept[i]
                  upper[, i] = upper[, i] + intercept[i]
                }
            }
            if (rug) {
                if (is.null(ylim)) 
                  ylim2 = range(ZZZ - intercept[1], lower - intercept[1], 
                    upper - intercept[1], na.rm = TRUE)
                matplot(cov[[k]], Z[[k]], type = "n", xlab = xlab[k], 
                  ylab = ylab, ylim = ylim2, ...)
                rug(cov[[k]])
                matlines(seq(min(cov[[k]]), max(cov[[k]]), length = 100), 
                  ZZZ - intercept[1], col = rainbow(np.plot + 
                    1)[1:np.plot], lty = 1, lwd = 2)
                if (ci) {
                  matlines(seq(min(cov[[k]]), max(cov[[k]]), 
                    length = 100), lower - intercept[1], col = rainbow(np.plot + 
                    1)[1:np.plot], lty = 2, lwd = 2)
                  matlines(seq(min(cov[[k]]), max(cov[[k]]), 
                    length = 100), upper - intercept[1], col = rainbow(np.plot + 
                    1)[1:np.plot], lty = 2, lwd = 2)
                }
            }
            else {
                if (is.null(ylim)) 
                  ylim2 = range(yy, Z[[k]], lower, upper, na.rm = TRUE)
                plot(cov[[k]], yy, cex = 0.5, pch = 20, col = "grey42", 
                  xlab = xlab[k], ylab = ylab, ylim = ylim2, 
                  ...)
                matlines(seq(min(cov[[k]]), max(cov[[k]]), length = 100), 
                  ZZZ, col = rainbow(np.plot + 1)[1:np.plot], 
                  lty = 1, lwd = 2)
                if (ci) {
                  matlines(seq(min(cov[[k]]), max(cov[[k]]), 
                    length = 100), lower, col = rainbow(np.plot + 
                    1)[1:np.plot], lty = 2, lwd = 2)
                  matlines(seq(min(cov[[k]]), max(cov[[k]]), 
                    length = 100), upper, col = rainbow(np.plot + 
                    1)[1:np.plot], lty = 2, lwd = 2)
                }
            }
            if (legend) 
                legend(x = "topright", pch = 19, cex = 0.8, col = rev(rainbow(np.plot + 
                  1)[1:np.plot]), legend = rev(pp[pp.plot]), 
                  bg = "white", bty = "n")
        }
        else if (types[[k]] == "markov") {
            z = NULL
            Zspathelp = helper[[k]][[2]]
            bnd = helper[[k]][[1]]
            if (inherits(x, "boost")) {
                for (i in 1:np) {
                  z = cbind(z, coefficients[[k]][, i] + intercept[i])
                }
                if (class(bnd) != "bnd") {
                  plot(seq(0, 1.1 * max(cov[[k]]), length = 10), 
                    seq(0, max(z), length = 10), type = "n", 
                    xlab = "Districts", ylab = "coefficients")
                  matpoints(cov[[k]], Z[[k]], col = rainbow(np.plot + 
                    1)[1:np.plot])
                  if (legend) 
                    legend(x = "right", pch = 19, cex = 1, col = rev(rainbow(np.plot + 
                      1)[1:np.plot]), legend = rev(pp[pp.plot]), 
                      bg = "white", bty = "n")
                }
                else {
                  plot.limits = range(coefficients[[k]])
                  for (i in 1:np.plot) {
                    re = data.frame(attr(bnd, "regions"), coefficients[[k]][, 
                      i])
                    drawmap(re, bnd, regionvar = 1, plotvar = 2, 
                      limits = plot.limits, main = pp[pp.plot[i]], 
                      swapcolors = TRUE, legend = legend)
                  }
                }
            }
            else {
                z = matrix(NA, nrow = nrow(Zspathelp), ncol = np)
                if (ci) {
                  for (i in 1:np) {
                    z[, i] = (Zspathelp %*% (coefficients[[k]][, 
                      i] - qnorm(0.975) * sqrt(diag(x$covmat[[i]][c(partbasis + 
                      1), c(partbasis + 1)]))) > 0) * 1 - 1 * 
                      (Zspathelp %*% (coefficients[[k]][, i] + 
                        qnorm(0.975) * sqrt(diag(x$covmat[[i]][c(partbasis + 
                          1), c(partbasis + 1)]))) < 0)
                  }
                }
                else {
                  coefficients[[k]] = Zspathelp %*% coefficients[[k]]
                  for (i in 1:np) {
                    z[, i] = coefficients[[k]][, i]
                  }
                }
                if (class(bnd) != "bnd") {
                  plot(seq(0, 1.1 * max(cov[[k]]), length = 10), 
                    seq(0, max(z[, pp.plot]), length = 10), type = "n", 
                    xlab = "District", ylab = "coefficients")
                  points(rep(as.numeric(attr(bnd, "regions")), 
                    times = np), z[, pp.plot], col = rainbow(np.plot + 
                    1)[1:np.plot])
                  if (legend) 
                    legend(x = "right", pch = 19, cex = 1, col = rev(rainbow(np.plot + 
                      1)[1:np.plot]), legend = rev(pp[pp.plot]), 
                      bg = "white", bty = "n")
                }
                else {
                  plot.limits = range(z)
                  for (i in 1:np.plot) {
                    re = data.frame(attr(bnd, "regions"), z[, 
                      pp.plot[i]])
                    drawmap(re, bnd, regionvar = 1, plotvar = 2, 
                      limits = plot.limits, main = pp[pp.plot[i]], 
                      swapcolors = TRUE, legend = legend, cex.legend = 1)
                  }
                }
            }
        }
        else if (types[[k]] == "2dspline") {
            if (inherits(x, "boost")) {
                for (i in 1:np) {
                  if (i %in% pp.plot) {
                    gitter = 20
                    x.min = apply(cov[[k]], 2, min, na.rm = TRUE)
                    x.max = apply(cov[[k]], 2, max, na.rm = TRUE)
                    z = x$plotpredict(k)[[i]]
                    if (is.null(ylim)) 
                      ylim2 = range(yy, z, na.rm = TRUE)
                    persp(seq(x.min[1], x.max[1], length = gitter), 
                      seq(x.min[2], x.max[2], length = gitter), 
                      z, ticktype = "detailed", phi = 40, theta = 35, 
                      zlim = ylim2, col = "lightblue", xlab = "X", 
                      ylab = "Y", zlab = ylab, main = pp[pp.plot[i]])
                  }
                }
            }
            else {
                gitter = 20
                x.min = apply(cov[[k]], 2, min, na.rm = TRUE)
                x.max = apply(cov[[k]], 2, max, na.rm = TRUE)
                x.gitter = cbind(rep(seq(x.min[1], x.max[1], 
                  length = gitter), times = gitter), rep(seq(x.min[2], 
                  x.max[2], length = gitter), each = gitter))
                ndat = as.data.frame(x.gitter)
                names(ndat) = rep(xlab[k], 2)
                B.gitter = predict(x$bases[[k]], ndat)
                for (i in 1:np) {
                  if (i %in% pp.plot) {
                    z <- B.gitter %*% coefficients[[k]][, i]
                    if (ci) {
                      for (nn in 1:nrow(B.gitter)) {
                        deviation = qnorm(0.975) * sqrt(t(c(1, 
                          B.gitter[nn, ])) %*% x$covmat[[i]][c(1, 
                          partbasis + 1), c(1, partbasis + 1)] %*% 
                          c(1, B.gitter[nn, ]))
                        z[nn] = 1 * ((z - deviation) > 0) - 1 * 
                          ((z + deviation) < 0)
                      }
                      z = t(matrix(z, nrow = gitter, ncol = gitter))
                      image(seq(x.min[1], x.max[1], length = gitter), 
                        seq(x.min[2], x.max[2], length = gitter), 
                        z, zlim = range(z), main = pp[pp.plot[i]])
                    }
                    else {
                      z = t(matrix(z, nrow = gitter, ncol = gitter))
                      if (is.null(ylim)) 
                        ylim2 = range(z, na.rm = TRUE)
                      persp(seq(x.min[1], x.max[1], length = gitter), 
                        seq(x.min[2], x.max[2], length = gitter), 
                        z, ticktype = "detailed", phi = 40, theta = 35, 
                        zlim = ylim2, col = "lightblue", xlab = "X", 
                        ylab = "Y", zlab = ylab, main = pp[pp.plot[i]])
                    }
                  }
                }
            }
        }
        else if (types[[k]] == "radial") {
            if (inherits(x, "boost")) {
                for (i in 1:np) {
                  if (i %in% pp.plot) {
                    z = x$plotpredict(k)[[i]]
                    if (is.null(ylim)) 
                      ylim2 = range(cbind(yy, z), na.rm = TRUE)
                    persp(seq(x.min[1], x.max[1], length = gitter), 
                      seq(x.min[2], x.max[2], length = gitter), 
                      z, ticktype = "detailed", phi = 40, theta = 35, 
                      zlim = ylim2, col = "lightblue", xlab = "X", 
                      ylab = "Y", zlab = ylab, main = pp[pp.plot[i]])
                  }
                }
            }
            else {
                gitter = 20
                x.min = apply(cov[[k]], 2, min)
                x.max = apply(cov[[k]], 2, max)
                x.gitter = cbind(rep(seq(x.min[1], x.max[1], 
                  length = gitter), times = gitter), rep(seq(x.min[2], 
                  x.max[2], length = gitter), each = gitter))
                cov[[k]] = cov[[k]][order(cov[[k]][, 1]), ]
                knots = helper[[k]]
                B.gitter = matrix(NA, nrow = dim(x.gitter)[1], 
                  ncol = dim(knots)[1])
                for (j in 1:dim(knots)[1]) {
                  r = sqrt(rowSums((x.gitter - matrix(unlist(knots[j, 
                    ]), nrow = nrow(x.gitter), ncol = ncol(knots), 
                    byrow = T))^2))
                  r[r == 0] = 1
                  B.gitter[, j] = r^2 * log(r)
                }
                par(mfrow = (c(row.grid, col.grid)))
                for (i in 1:np) {
                  if (i %in% pp.plot) {
                    z <- B.gitter %*% coefficients[[k]][, i] + 
                      intercept[i]
                    z = t(matrix(z, nrow = gitter, ncol = gitter))
                    if (is.null(ylim)) 
                      ylim2 = range(cbind(yy, z), na.rm = TRUE)
                    persp(seq(x.min[1], x.max[1], length = gitter), 
                      seq(x.min[2], x.max[2], length = gitter), 
                      z, ticktype = "detailed", phi = 40, theta = 35, 
                      zlim = ylim2, col = "lightblue", xlab = "X", 
                      ylab = "Y", zlab = ylab, main = pp[pp.plot[i]])
                  }
                }
            }
        }
        else if (types[[k]] == "krig") {
            gitter = 20
            krig.phi = helper[[k]][[1]]
            x.min = apply(cov[[k]], 2, min)
            x.max = apply(cov[[k]], 2, max)
            x.gitter = cbind(rep(seq(x.min[1], x.max[1], length = gitter), 
                times = gitter), rep(seq(x.min[2], x.max[2], 
                length = gitter), each = gitter))
            cov[[k]] = cov[[k]][order(cov[[k]][, 1]), ]
            knots = helper[[k]][[2]]
            B.gitter = matrix(NA, nrow = dim(x.gitter)[1], ncol = dim(knots)[1])
            for (j in 1:dim(knots)[1]) {
                r = sqrt(rowSums((x.gitter - matrix(unlist(knots[j, 
                  ]), nrow = nrow(x.gitter), ncol = ncol(knots), 
                  byrow = T))^2))/krig.phi
                B.gitter[, j] = exp(-r) * (1 + r)
            }
            for (i in 1:np) {
                if (i %in% pp.plot) {
                  z <- B.gitter %*% coefficients[[k]][, i] + 
                    intercept[i]
                  z = t(matrix(z, nrow = gitter, ncol = gitter))
                  if (is.null(ylim)) 
                    ylim2 = range(cbind(yy, z), na.rm = TRUE)
                  persp(seq(x.min[1], x.max[1], length = gitter), 
                    seq(x.min[2], x.max[2], length = gitter), 
                    z, ticktype = "detailed", phi = 40, theta = 35, 
                    zlim = ylim2, col = "lightblue", xlab = "X", 
                    ylab = "Y", zlab = ylab, main = pp[pp.plot[i]])
                }
            }
        }
        else if (types[[k]] == "random") {
            if (inherits(x, "boost")) {
                matplot(cov[[k]], Z[[k]], col = rainbow(np.plot + 
                  1)[1:np.plot], xlab = xlab[k], ylab = ylab, 
                  pch = 15)
                legend(x = "right", pch = 19, cex = 1, col = rev(rainbow(np.plot + 
                  1)[1:np.plot]), legend = rev(pp[pp.plot]), 
                  bg = "white", bty = "n")
            }
            else {
                plot(seq(0, 1.1 * max(cov[[k]]), length = 10), 
                  seq(0, max(coefficients[[k]] + intercept - 
                    intercept[1]), length = 10), type = "n", 
                  xlab = xlab[k], ylab = "coefficients")
                points(rep(sort(unique(cov[[k]])), times = np.plot), 
                  (coefficients[[k]] + intercept - intercept[1])[, 
                    pp.plot], col = rainbow(np.plot + 1)[1:np.plot])
                if (legend) 
                  legend(x = "right", pch = 19, cex = 1, col = rev(rainbow(np.plot + 
                    1)[1:np.plot]), legend = rev(pp[pp.plot]), 
                    bg = "white", bty = "n")
            }
        }
        else if (types[[k]] == "ridge") {
            plot(seq(0, 1.1 * dim(cov[[k]])[2], length = 10), 
                seq(0, max(coefficients[[k]] + intercept - intercept[1]), 
                  length = 10), type = "n", xlab = xlab[k], ylab = "coefficients")
            points(rep(1:dim(cov[[k]])[2], times = np.plot), 
                (coefficients[[k]] + intercept - intercept[1])[, 
                  pp.plot], col = rainbow(np.plot + 1)[1:np.plot])
            if (legend) 
                legend(x = "right", pch = 19, cex = 1, col = rev(rainbow(np.plot + 
                  1)[1:np.plot]), legend = rev(pp[pp.plot]), 
                  bg = "white", bty = "n")
        }
        else if (types[[k]] == "parametric" && !is.factor(cov[[k]][, 
            1])) {
            if (inherits(x, "boost")) {
                plot(cov[[k]], yy, cex = 0.5, pch = 20, col = "grey42", 
                  xlab = names(cov)[k], ylab = attr(yy, "name"), 
                  ylim = range(cbind(yy, Z[[k]])))
                matlines(sort(cov[[k]]), Z[[k]][order(cov[[k]]), 
                  pp.plot], col = rainbow(np.plot + 1)[1:np.plot], 
                  lty = 1)
                if (legend) 
                  legend(x = "bottomright", pch = 19, cex = 1, 
                    col = rev(rainbow(np.plot + 1)[1:np.plot]), 
                    legend = rev(pp[pp.plot]), bg = "white", 
                    bty = "n")
            }
            else {
                if (is.null(ylim)) 
                  ylim2 = range(yy, Z[[k]], na.rm = TRUE)
                for (i in 1:ncol(cov[[k]])) {
                  if (rug) {
                    matplot(cov[[k]][, i], Z[[k]] - intercept[1], 
                      type = "n", xlab = xlab[k], ylab = ylab, 
                      ylim = range(Z[[k]] - intercept[1], na.rm = TRUE), 
                      ...)
                    rug(cov[[k]][, i])
                    matlines(cov[[k]][, i], Z[[k]][, pp.plot] - 
                      intercept[1], col = rainbow(np.plot + 1)[1:np.plot], 
                      lty = 1, lwd = 2)
                  }
                  else {
                    plot(cov[[k]][, i], yy, cex = 0.5, pch = 20, 
                      col = "grey42", xlab = names(cov)[k], ylab = ylab, 
                      ylim = ylim2)
                    matlines(sort(cov[[k]][, i]), Z[[k]][order(cov[[k]][, 
                      i]), pp.plot], col = rainbow(np.plot + 
                      1)[1:np.plot], lty = 1, lwd = 2)
                  }
                  if (legend) 
                    legend(x = "topright", pch = 19, cex = 1, 
                      col = rev(rainbow(np.plot + 1)[1:np.plot]), 
                      legend = rev(pp[pp.plot]), bg = "white", 
                      bty = "n")
                }
            }
        }
        else if (types[[k]] == "special") {
            if (is.null(ylim)) 
                ylim2 = range(cbind(yy, Z[[k]]), na.rm = TRUE)
            plot(cov[[k]], yy, cex = 0.5, pch = 20, col = "grey42", 
                xlab = xlab[k], ylab = ylab, ylim = ylim2, ...)
            matlines(sort(cov[[k]]), Z[[k]][order(cov[[k]]), 
                pp.plot], col = rainbow(np.plot + 1)[1:np.plot], 
                lty = 1, lwd = 2)
            if (legend) 
                legend(x = "topright", pch = 19, cex = 1, col = rev(rainbow(np.plot + 
                  1)[1:np.plot]), legend = rev(pp[pp.plot]), 
                  bg = "white", bty = "n")
        }
    }
}
