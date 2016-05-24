plotSimDiags <-
function (obj, simvalues = NULL, seed = NULL, types = NULL, which = c(1:3,  5),
              layout = c(4, 1), qqline = TRUE, cook.levels = c(0.5, 1),
              caption = list("Residuals vs Fitted", "Normal Q-Q", "Scale-Location",
              "Cook's distance", "Residuals vs Leverage",
              expression("Cook's dist vs Leverage  " * h[ii]/(1 - h[ii]))), ...)
{
    dropInf <- function(x, h) {
        if (any(isInf <- h >= 1)) {
            x[isInf] <- NaN
        }
        x
    }
    if (!inherits(obj, "lm"))
        stop("use only with \"lm\" objects")
    if (!is.numeric(which) || any(which < 1) || any(which > 6))
        stop("'which' must be in 1:6")
    gphlist <- vector("list", 6)
    names(gphlist) <- c("residVSfitted", "normalQQ", "scaleVSloc",
                        "CookDist", "residVSlev", "CookVSlev")
    isGlm <- inherits(obj, "glm")
    if (is.null(types))
        types <- list(c("p", "smooth"), NULL, c("p", "smooth"),
                      "h", c("p", "smooth"), NULL)
    show <- rep(FALSE, 6)
    show[which] <- TRUE
    numsim <- prod(layout)
    if (is.null(simvalues))
        simvalues <- simulate(obj, nsim = numsim, seed = seed)
    if (ncol(simvalues) != numsim)
        stop(paste("Number of columns of simvalues must", "equal number of panels of layout"))
    okrows <- complete.cases(residuals(obj))
    hat <- fitted(obj)
    nobs <- length(hat)
    df <- as.data.frame(simvalues)
    simnam <- paste("Sim", 1:numsim, sep = "_")
    names(df)[1:numsim] <- simnam
    mmat <- model.matrix(obj)
    regs <- lapply(df, function(y) lm(y ~ mmat))
    hii <- lm.influence(obj, do.coef = FALSE)$hat
    diaglist <- lapply(regs, lmdiags, which = which, hii=hii)
    rslist <- sapply(diaglist, function(x) x[["rs"]])
    objdf <- data.frame(gp = paste("Simulation", rep(1:numsim,
                        rep(nobs, numsim)), sep = "_"))
    if (show[1]) {
        objdf[["r"]] <- as.vector(sapply(diaglist, function(x) x[["r"]]))
        objdf[["yh"]] <- as.vector(sapply(diaglist, function(x) x[["yh"]]))
    }
    if (any(show[2:3])) {
        objdf[["rs"]] <- as.vector(sapply(diaglist, function(x) x[["rs"]]))
        ylab23 <- if (isGlm)
            "Std. deviance resid."
        else "Standardized residuals"
    }
    if (show[3])
        objdf[["yhn0"]] <- as.vector(sapply(diaglist, function(x) x$yhn0))
    if (any(show[c(4, 6)]))
        objdf[["cook"]] <- as.vector(sapply(diaglist, function(x) x$cook))
    if (show[5]){
        objdf[["rsp"]] <- as.vector(sapply(diaglist, function(x) x[["rsp"]]))
        }
    getCaption <- function(k) if (length(caption) < k)
        NA_character_
    else as.graphicsAnnot(caption[[k]])
    if (show[1]) {
        formyx <- r ~ yh | gp
        gph <- xyplot(formyx, type = types[[1]], par.settings = simpleTheme(pch = c(16,
                                                                            16), lty = 2, col = c("black", "gray")), layout = layout,
                      data = objdf, xlab = "Fitted values", ylab = "Residuals",
                      main = getCaption(1), ...)
        gph <- gph + latticeExtra::layer(lattice::panel.abline(h = 0, lty = 3, col = "gray"))
        gphlist[[1]] <- gph
    }
    if (show[2]) {
        gph <- lattice::qqmath(~rs | gp, data = objdf, prepanel = lattice::prepanel.qqmathline,
                      panel = function(x, ...) {
                          lattice::panel.qqmathline(x, lty = 2, ...)
                          lattice::panel.qqmath(x, ...)
                      }, layout = layout, xlab = "Theoretical Quantiles",
                      ylab = ylab23, ...)
        gphlist[[2]] <- gph
    }
    if (show[3]) {
        yl <- as.expression(substitute(sqrt(abs(YL)), list(YL = as.name(ylab23))))
        sqrtabsr <- sqrt(abs(objdf[["rs"]]))
        formyx <- sqrtabsr ~ yhn0 | gp
        gph <- xyplot(formyx, type = types[[3]], par.settings = simpleTheme(pch = c(16,
                                                                            16), lty = 2, col = c("black", "gray")), layout = layout,
                      data = objdf, xlab = "Fitted values", ylab = yl,
                      ...)
        gphlist[[3]] <- gph
    }
    if (show[4]) {
        objdf$x <- rep(1:nobs, numsim)
        gph <- xyplot(cook ~ x | gp, data = objdf, type = types[[4]],
                      layout = layout, xlab = "Obs. number", ylab = "Cook's distance",
                      ...)
        gphlist[[4]] <- gph
    }
    if (show[5]) {
        ylab5 <- if (isGlm)
            "Std. Pearson resid."
        else "Standardized residuals"
        panel5.diag <- function(x, y, ...) {
            lattice::panel.xyplot(x, y, ...)
            lattice::panel.abline(h = 0, lty = 3, col = "gray")
        }

        r.hat <- range(hii, na.rm = TRUE)
        isConst.hat <- all(r.hat == 0) || diff(r.hat) < 1e-10 *
            mean(hii, na.rm = TRUE)
        if (isConst.hat) {
            aterms <- attributes(terms(obj))
            dcl <- aterms$dataClasses[-aterms$response]
            facvars <- names(dcl)[dcl %in% c("factor", "ordered")]
            mf <- model.frame(obj)[facvars]
            if (ncol(mf) > 0) {
                dm <- data.matrix(mf)
                nf <- length(nlev <- unlist(unname(lapply(obj$xlevels,
                                                          length))))
                ff <- if (nf == 1)
                    1
                else rev(cumprod(c(1, nlev[nf:2])))
                facval <- (dm - 1) %*% ff
                levels(facval) <- obj$xlevels[[1L]]
                xlim <- c(-1/2, sum((nlev-1)*ff)+1/2)
                vval <- ff[1L] * (0:(nlev[1L]-1))-1/2
                gph <- xyplot(rsp ~ facval|gp, xlim=xlim, ylab=ylab5, data = objdf,
                              ## panel.groups=function(x,y,...){
                              ##     panel.points(x,y,...)
                              ##     xu <- unique(x)
                              ##     ym <- sapply(split(y,x),mean)
                              ##     browser()
                              ##     panel.points(xu,ym, pch=2, col="red")
                              ## },
                              scales=list(x=list(at=0:0:(nlev[1L]-1), labels=obj$xlevels[[1L]])))
                gph <- gph + latticeExtra::layer(lattice::panel.abline(v=vval, col = "gray", lty = "F4"),
                                   lattice::panel.abline(h = 0, lty = 3, col = "gray"),
                                   data=list(vval=vval))
            }
            else {
                message("hat values (leverages) are all = ",
                        format(mean(r.hat)), "\n and there are no factor predictors; no plot no. 5")
                do.plot <- FALSE
            }
        } else {
        xx <- rep(hii,numsim)
        xx[xx >= 1] <- NA
        p <- length(coef(obj))
        if (length(cook.levels))
            yscale.cpts <- function(lim, ...) {
                ans <- lattice::yscale.components.default(lim = lim, ...)
                ans$right <- ans$left
                ans$right$ticks$at <- c(-rev(sqrt(cook.levels)) *
                                        ymult, sqrt(cook.levels) * ymult)
                ans$right$ticks$tck <- rep(0, 2 * length(cook.levels))
                ans$right$labels$at <- c(-rev(sqrt(cook.levels)) *
                                         ymult, sqrt(cook.levels) * ymult)
                ans$right$labels$labels <- paste(c(rev(cook.levels),
                                                   cook.levels))
                ans
            }
        formyx <- rsp ~ xx | gp
        gph <- xyplot(formyx, type = types[[5]], data = objdf,
                      par.settings = simpleTheme(pch = c(16, 16), lty = 2,
                      col = c("black", "gray")), scales = list(y = list(alternating = 3)),
                      layout = layout, xlab = "Leverage", ylab = ylab5,
                      main = getCaption(5), panel = panel5.diag, yscale.components = yscale.cpts)
        usr <- gph[["x.limits"]]
        xmax <- min(0.99, usr[2L])
        ymult <- sqrt(p * (1 - xmax)/xmax)
        hh <- seq.int(min(r.hat[1L], r.hat[2L]/100), usr[2L],
                      length.out = 101)
        xy <- expand.grid(hh = c(hh, NA), cl.h = cook.levels)
        xy <- within(xy, cl.h <- sqrt(cl.h * p * (1 - hh)/hh))
        xy <- with(xy, data.frame(hh = c(hh, hh), cl.h = c(cl.h,
                                                  -cl.h)))
        aty <- c(-rev(sqrt(cook.levels)) * ymult, sqrt(cook.levels) *
                 ymult)
        laby <- paste(c(rev(cook.levels), cook.levels))
        gph2 <- xyplot(cl.h ~ hh, data = xy, type = "l", lty = 3,
                       col = "red")
        gph <- gph + latticeExtra::as.layer(gph2)
    }
        gphlist[[5]] <- gph
    }
    if (show[6]) {
        g <- with(objdf, dropInf(hii/(1 - hii), hii))
        ymx <- with(objdf, max(cook, na.rm = TRUE) * 1.025)
        athat <- pretty(hii)
        gph <- xyplot(cook ~ g | gp, xlim = c(0, max(g, na.rm = TRUE)),
                      data = objdf, ylim = c(0, ymx))
        p <- length(coef(obj))
        bval <- with(objdf, pretty(sqrt(p * cook/g), 5))
        xmax <- gph[["x.limits"]][2]
        ymax <- gph[["y.limits"]][2]
        panel6 <- function(x, y, ...) {
            lattice::panel.xyplot(x, y, ...)
            for (i in seq_along(bval)) {
                bi2 <- bval[i]^2
                if (ymax > bi2 * xmax) {
                    xi <- xmax
                    yi <- bi2 * xi
                    lattice::panel.abline(0, bi2, lty = 2)
                    lattice::panel.text(xi, yi, paste(bval[i]), adj = c(1.25,
                                                       0.5), cex = 0.75)
                }
                else {
                    yi <- ymax
                    xi <- yi/bi2
                    lattice::panel.lines(c(0, xi), c(0, yi), lty = 2)
                    lattice::panel.text(xi, ymax, paste(bval[i]), adj = c(0.5,
                                                         1.25), cex = 0.75)
                }
            }
        }
        gph <- xyplot(cook ~ g | gp, xlim = c(0, max(g, na.rm = TRUE)),
                      data = objdf, ylim = c(0, ymx), main = getCaption(6),
                      ylab = "Cook's distance", xlab = expression("Leverage  " *
                                                h[ii]), layout = layout, scales = list(x = list(at = athat/(1 -
                                                                                                athat), labels = paste(athat))), panel = panel6)
        gphlist[[6]] <- gph
    }
    gphlist <- gphlist[!sapply(gphlist, is.null)]
    if (length(gphlist) == 1)
        gphlist <- gphlist[[1]]
    gphlist
}
