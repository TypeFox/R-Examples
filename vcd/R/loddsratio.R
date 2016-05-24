## Modifications:
## -- return a dnames component, containing dimnames for the array version of coef
## -- added dim methods: dim.loddsratio, dimnames.loddsratio
## -- added print.loddsratio
## -- handle strata:  result computed correctly, but structure of coef() loses names
##    and confint doesn't work in the 2x2xk or RxCxk case

## -- Fixed problem with strata by setting rownames and colnames for contrast matrix
## DONE: handle multiple strata (|foo:bar)
## -- print.loddsratio now uses drop() for dimensions of length 1
## -- made generic, anticipating a formula method, maybe structable or ftable methods
## DONE: decide which methods should allow a log=FALSE argument to provide exp(lor)

## -- Now handle any number of strata
## -- Added log= argument to print, coef methods, and added confint.loddsratio method,
##    allowing log=FALSE

## -- Incorporated Z code additions, fixing some <FIXME>s
## -- Added as.matrix and as.array methods; had to make as.array S3 generic
## -- Added header to print method
## -- Added as.data.frame method (for use in plots)
## --   "LOR" is renamed "OR" if log=FALSE
## -- Revised as.matrix to drop leading 1:2 dimensions of length 1
## -- Removed as.array generic, now in base

## -- DM: added plot.oddsratio method
## -- DM: added formula interface
## -- DM: add t() and aperm() methdos

loddsratio <- function(x, ...)
    UseMethod("loddsratio")

loddsratio.formula <-
    function(formula, data = NULL, ..., subset = NULL, na.action = NULL)
{
    m <- match.call(expand.dots = FALSE)
    edata <- eval(m$data, parent.frame())

    fstr <- strsplit(paste(deparse(formula), collapse = ""), "~")
    vars <- strsplit(strsplit(gsub(" ", "", fstr[[1]][2]), "\\|")[[1]], "\\+")
    varnames <- vars[[1]]

    condnames <- if (length(vars) > 1) vars[[2]] else NULL

    dep <- gsub(" ", "", fstr[[1]][1])
    if (!dep %in% c("","Freq")) {
        if (all(varnames == ".")) {
            varnames <- if (is.data.frame(data))
                            colnames(data)
                        else
                            names(dimnames(as.table(data)))
            varnames <- varnames[-which(varnames %in% dep)]
        }

        varnames <- c(dep, varnames)
    }


    if (inherits(edata, "ftable") || inherits(edata, "table") || length(dim(edata)) > 2) {
        condind <- NULL
        dat <- as.table(data)
        if(all(varnames != ".")) {
            ind <- match(varnames, names(dimnames(dat)))
            if (any(is.na(ind)))
                stop(paste("Can't find", paste(varnames[is.na(ind)], collapse=" / "), "in", deparse(substitute(data))))

            if (!is.null(condnames)) {
                condind <- match(condnames, names(dimnames(dat)))
                if (any(is.na(condind)))
                    stop(paste("Can't find", paste(condnames[is.na(condind)], collapse=" / "), "in", deparse(substitute(data))))
                ind <- c(ind, condind)
            }
            dat <- margin.table(dat, ind)
        }
        loddsratio.default(dat, strata = if (is.null(condind)) NULL else match(condnames, names(dimnames(dat))), ...)
    } else {
        m <- m[c(1, match(c("formula", "data", "subset", "na.action"), names(m), 0))]
        m[[1]] <- as.name("xtabs")
        m$formula <-
            formula(paste(if("Freq" %in% colnames(data)) "Freq",
                          "~",
                          paste(c(varnames, condnames), collapse = "+")))
        tab <- eval(m, parent.frame())
        loddsratio.default(tab, ...)
    }
}

loddsratio.default <- function(x, strata = NULL, log = TRUE,
                               ref = NULL, correct = any(x == 0), ...)
{
    ## check dimensions
    L <- length(d <- dim(x))
    if(any(d < 2L)) stop("All table dimensions must be 2 or greater")
    if(L > 2L & is.null(strata)) strata <- 3L:L
    if(is.character(strata)) strata <- which(names(dimnames(x)) == strata)
    if(L - length(strata) != 2L) stop("All but 2 dimensions must be specified as strata.")

    ## dimensions of primary R x C table
    dp <- if (length(strata)) d[-strata] else d
    dn <- if (length(strata)) dimnames(x)[-strata] else dimnames(x)
    R <- dp[1]
    C <- dp[2]
                                        # shadow matrix with proper dimnames
    X <- matrix(0, R, C, dimnames=dn)

    ## process reference categories (always return list of length
    ## two with reference for rows/cols, respectively)
    if(is.null(ref)) {
        ref <- list(NULL, NULL)
    } else if(is.character(ref)) {
        if(length(ref) != 2L) stop("'ref' must specify both reference categories")
        ref <- list(match(ref[1L], rownames(x)), match(ref[2L], colnames(x)))
    } else if(is.numeric(ref)) {
        ref <- as.integer(rep(ref, length.out = 2L))
        ref <- list(ref[1L], ref[2L])
    }

    ## compute corresponding indices
    compute_index <- function(n, ref) {
        if(is.null(ref)) return(cbind(1:(n-1), 2:n))
        rval <- cbind(ref, 1:n)
        d <- rval[,2L] - rval[,1L]
        rval <- rbind(
            rval[d > 0, 1:2],
            rval[d < 0, 2:1]
        )
        return(rval[order(rval[,1L]),,drop = FALSE])
    }
    Rix <- compute_index(R, ref[[1L]])
    Cix <- compute_index(C, ref[[2L]])

    ## set up contrast matrix for the primary R x C table
    contr <- matrix(0L, nrow = (R-1) * (C-1), ncol = R * C)
    colnames(contr) <- paste(rownames(X)[as.vector(row(X))], colnames(X)[as.vector(col(X))], sep = ":")
    rownames(contr) <- rep("", (R-1) * (C-1))
    for(i in 1:(R-1)) for(j in 1:(C-1)) {
        rix <- (j-1) * (R-1) + i
        cix <- rep(Rix[i,], 2L) + R * (rep(Cix[j,], each = 2L) - 1L)
        contr[rix, cix] <- c(1L, -1L, -1L, 1L)
        rownames(contr)[rix] <- sprintf("%s/%s",
                                        paste(rownames(X)[Rix[i,]], collapse = ":"),
                                        paste(colnames(X)[Cix[j,]], collapse = ":"))
    }

                                        # handle strata
    if (!is.null(strata)) {
  	if (length(strata)==1) {
            sn <- dimnames(x)[[strata]]
        }
        else {
            sn <- apply(expand.grid(dimnames(x)[strata]), 1, paste, collapse = ":")
        }
        rn <- as.vector(outer( dimnames(contr)[[1]], sn, paste, sep='|'))
        cn <- as.vector(outer( dimnames(contr)[[2]], sn, paste, sep='|'))
        contr <- kronecker(diag(prod(dim(x)[strata])), contr)
        rownames(contr) <- rn
        colnames(contr) <- cn
    }

    ## dimnames for array version
    dn <- list(rep("", R-1), rep("", C-1))
    for(i in 1:(R-1)) dn[[1]][i] <- paste(rownames(x)[Rix[i,]], collapse = ":")
    for(j in 1:(C-1)) dn[[2]][j] <- paste(colnames(x)[Cix[j,]], collapse = ":")
    if (!is.null(strata)) dn <- c(dn, dimnames(x)[strata])
    if (!is.null(names(dimnames(x)))) names(dn) <- names(dimnames(x))

    ## point estimates

    if (is.logical(correct)) {
        add <- if(correct) 0.5 else 0
    }
    else if(is.numeric(correct)) {
        add <- as.vector(correct)
        if (length(add) != length(x))
            stop("array size of 'correct' does not conform to the data")
    }
    else stop("correct is not valid")

    ##coef <- drop(contr %*% log(as.vector(x) + add))
    ##FIXME: 0 cells mess up the matrix product, try workaround:
    mat <- log(as.vector(x) + add) * t(contr)
    nas <- apply(contr != 0 & is.na(t(mat)), 1, any)
    coef <- apply(mat, 2, sum, na.rm = TRUE)
    coef[nas] <- NA
    ## covariances
    ##vcov <- crossprod(diag(sqrt(1/(as.vector(x) + add))) %*% t(contr))
    tmp <- sqrt(1/(as.vector(x) + add)) * t(contr)
    tmp[is.na(tmp)] <- 0
    vcov <- crossprod(tmp)
    vcov[nas,] <- NA
    vcov[,nas] <- NA

    rval <- structure(list(
        coefficients = coef,
        dimnames = dn,
        dim = as.integer(sapply(dn, length)),
        vcov = vcov,
        contrasts = contr,
        log = log
    ), class = "loddsratio")
    rval
}

## dim methods
dimnames.loddsratio <- function(x, ...) x$dimnames
dim.loddsratio <- function(x, ...) x$dim

## t/aperm-methods
t.loddsratio <- function(x)
    aperm(x)

aperm.loddsratio <- function(a, perm = NULL, ...)
{
    d <- length(a$dim)
    if(is.null(perm)) {
        perm <- if (d < 3) 2L : 1L else c(2L : 1L, d : 3L)
    } else {
        if (any(perm[1:2] > 2L) || (d > 2L) && any(perm[-c(1:2)] < 2L))
            stop("Mixing of strata and non-strata variables not allowed!")
    }
    nams <- names(a$coefficients)
    a$coefficients <- as.vector(aperm(array(a$coef, dim = a$dim),
                                      perm, ...))
    nams <- as.vector(aperm(array(nams, dim = a$dim), perm, ...))
    names(a$coefficients) <- nams
    a$dimnames <- a$dimnames[perm]
    a$dim <- a$dim[perm]
    a$vcov <- a$vcov[nams, nams]
    a$contrasts <- a$contrasts[nams,]
    a
}

## straightforward methods
coef.loddsratio <- function(object, log = object$log, ...)
    if(log) object$coefficients else exp(object$coefficients)

vcov.loddsratio <- function(object, log = object$log, ...)
    if(log) object$vcov else `diag<-`(object$vcov, diag(object$vcov) * exp(object$coefficients)^2)

confint.loddsratio <-
    function(object, parm, level = 0.95, log = object$log, ...) {
        if (log) confint.default(object, parm = parm, level = level, ... )
        else {
            object$log = TRUE
            exp(confint.default(object, parm = parm, level = level, ... ))
        }
    }


make_header <- function(x)
{
    vn <- names(dimnames(x))
    header <- c(if(x$log) "log" else "",
                "odds ratios for", vn[1], "and", vn[2],
                if (length(vn)>2) c("by", paste(vn[-(1:2)], collapse=', ')), "\n\n")
    paste(header, sep = " ")
}

## print method
print.loddsratio <- function(x, log = x$log, ...) {
    cat(make_header(x))
    print(drop(array(coef(x, log = log), dim = dim(x), dimnames = dimnames(x)), ...))
    invisible(x)
}

summary.loddsratio <- function(object, ...)
    lmtest::coeftest(object, ...)

## reshape coef() methods
as.matrix.loddsratio <- function (x, log=x$log, ...) {
    Coef <- coef(x, log = log)
    if (length(dim(x))==2) matrix(Coef, ncol = dim(x)[2], dimnames=dimnames(x))
    else {  # drop leading dimensions with length 1, then reshape
        ddim <- which(dim(x)[1:2]==1)
        dim(Coef) <- dim(x)[-ddim]
        dimnames(Coef) <- dimnames(x)[-ddim]
        if (length(dim(Coef))==1) Coef
        else
            matrix(Coef, ncol = prod(dim(Coef)[-1]),
                   dimnames=list(dimnames(Coef)[[1]], apply(expand.grid(dimnames(Coef)[[-1]]), 1, paste, collapse = ":")))
    }
}


as.array.loddsratio <- function (x, log=x$log, ...) {
    res <- array(coef(x, log = log), dim = dim(x), dimnames=dimnames(x))
    drop(res)
}

as.data.frame.loddsratio <- function(x, row.names = NULL, optional, log=x$log, ...) {
    df <-data.frame(expand.grid(dimnames(x)),
                    LOR = coef(x, log=log),
                    ASE = sqrt(diag(vcov(x, log=log))), row.names=row.names,  ...
                    )
    if (!log) colnames(df)[ncol(df)-1] <- "OR"
    df
}

image.loddsratio <-
    function(x, interpolate = NULL, legend = legend_fixed,
             gp = shading_Friendly, gp_args = NULL,
             labeling = labeling_values("residuals", suppress = 0),
             perm = NULL, ...)
{
    a <- as.array(x)
    if (!is.null(dim(a))) {
        if (is.null(perm)) {
            d <- seq_along(dim(a))
            perm <- c(d[-c(1:2)], 1:2)
        }
        a <- aperm(a, perm)
    } else {
        a <- as.table(a)
        names(dimnames(a)) <- names(dimnames(x))[1]
    }
    if (is.null(interpolate))
        interpolate <- seq(0.1, max(abs(a), length.out = 4))
    if (is.null(gp_args))
        gp_args <- list(interpolate = interpolate)
    tmp <- a
    tmp[] <- 1
    mosaic(tmp, type = "expected",
           residuals = a, shade = TRUE,
           gp = shading_Friendly, gp_args = gp_args,
           legend = legend,
           labeling = labeling, ...)
}

tile.loddsratio <-
    function(x, interpolate = NULL, legend = legend_fixed,
             gp = shading_Friendly, gp_args = NULL,
             labeling = labeling_values("residuals", suppress = 0),
             halign = "center", valign = "center", perm = NULL, ...)
{
    a <- as.array(x)
    if (!is.null(dim(a))) {
        if (is.null(perm)) {
            d <- seq_along(dim(a))
            perm <- c(d[-c(1:2)], 1:2)
        }
        a <- aperm(a, perm)
    } else {
        a <- as.table(a)
        names(dimnames(a)) <- names(dimnames(x))[1]
    }

    if (is.null(interpolate))
        interpolate <- seq(0.1, max(abs(a), length.out = 4))
    if (is.null(gp_args))
        gp_args <- list(interpolate = interpolate)

    tile(abs(a), halign = halign, valign = valign,
         residuals = a, shade = TRUE,
         gp = shading_Friendly, gp_args = gp_args,
         legend = legend,
         labeling = labeling, ...)
}

"plot.loddsratio" <-
    function(x,
             baseline = TRUE,
             gp_baseline = gpar(lty = 2),
             lines = TRUE,
             lwd_lines = 3,
             confidence = TRUE,
             conf_level = 0.95,
             lwd_confidence = 2,
             whiskers = 0,
             transpose = FALSE,
             col = NULL,
             cex = 0.8,
             pch = NULL,
             bars = NULL,
             gp_bars = gpar(fill = "lightgray", alpha = 0.5),
             bar_width = unit(0.05, "npc"),

             legend = TRUE, legend_pos = "topright", legend_inset = c(0, 0),
             legend_vgap = unit(0.5, "lines"),
             gp_legend_frame = gpar(lwd = 1, col = "black"),
             gp_legend_title = gpar(fontface = "bold"),
             gp_legend = gpar(), legend_lwd = 1,
             legend_size = 1,

             xlab = NULL,
             ylab = NULL,
             xlim = NULL,
             ylim = NULL,

             main = NULL,
             gp_main = gpar(fontsize = 12, fontface = "bold"),
             newpage = TRUE, pop = FALSE, return_grob = FALSE,
             prefix = "",
             ...)
{
    ## handle default values, limits etc.
    LOG <- x$log
    values <- as.array(x)
    d <- dim(values)
    if (is.null(bars))
        bars <- is.null(d)
    oddsrange <- range(values, na.rm = TRUE)

    if(confidence) {
        CI  <- confint(x, log = LOG, level = conf_level)
        lwr <- CI[,1]
        upr <- CI[,2]
        oddsrange <- if (baseline)
                         c(min(0, lwr, na.rm = TRUE), max(0, upr, na.rm = TRUE))
                     else
                         c(min(lwr, na.rm = TRUE), max(upr, na.rm = TRUE))
    }

    if (is.null(main))
        main <- paste(make_header(x), collapse = " ")

    if (is.null(xlim))
        xlim <- if (is.null(d))
                    c(1, length(values))
                else
                    c(1, d[1])

    if (is.null(ylim))
        ylim <- oddsrange

    ylimaxis <- ylim + c(-1, 1) * diff(ylim) * 0.04
    xlimaxis <- xlim + c(-1, 1) * diff(xlim) * 0.04

    ncols <- if (is.null(d)) 1 else prod(d[-1])
    if (is.null(col))
        col <- rainbow_hcl(ncols, l = 50)

    if (is.null(pch))
        pch <- c(19,15,17, 1:14, 16, 18, 20:25)

    labs <- if (is.null(d)) names(values) else dimnames(values)[[1]]
    if (is.null(xlab))
        xlab <- if (is.null(d)) names(dimnames(x))[3] else names(dimnames(values))[1]
    if (is.null(ylab))
        ylab <- paste(if (LOG) "L" else "",
                      "OR(", paste(names(dimnames(x))[1:2], collapse = " / "), ")", sep = "")


    if (newpage) grid.newpage()
    if (transpose) {
        ## set up plot region, similar to plot.xy()
        pushViewport(plotViewport(xscale = ylimaxis, yscale = xlimaxis,
                                  default.units = "native",
                                  name = paste(prefix,"oddsratio_plot")))
        grid.yaxis(name = "yaxis", seq_along(labs), labs,
                   edits = gEdit("labels", rot = 90, hjust = .5, vjust = 0))
        grid.xaxis()
        grid.text(ylab, y = unit(-3.5, "lines"))
        grid.text(xlab, x = unit(-3, "lines"), rot = 90)
        grid.text(main, y = unit(1, "npc") + unit(1, "lines"), gp = gp_main)
        pushViewport(viewport(xscale = ylimaxis, yscale = xlimaxis,
                              default.units = "native", clip = "on"))

        ## baseline
        if (baseline)
            grid.lines(unit(c(1,1) - LOG, "native"),
                       unit(c(0,1), "npc"),
                       gp = gp_baseline)

        # workhorse for one stratum
        draw_one_stratum <- function(vals, pch = "o", col = "black", offset = 0,
                                     jitter = 0) {
            if (bars) {
                if (any(vals > !LOG))
                    grid.rect(unit(vals[vals > !LOG], "native"),
                              unit(seq_along(vals)[vals > !LOG], "native"),
                              height = bar_width,
                              width = unit(vals[vals > !LOG] - !LOG, "native"),
                              just = "right",
                              gp = gp_bars
                              )
                if (any(vals < !LOG))
                    grid.rect(unit(vals[vals < !LOG], "native"),
                              unit(seq_along(vals)[vals < !LOG], "native"),
                              height = bar_width,
                              width = unit(abs(vals[vals < !LOG] - !LOG), "native"),
                              just = "left",
                              gp = gp_bars
                              )
            }
            if (lines)
                grid.lines(unit(vals, "native"),
                           unit(seq_along(vals), "native"),
                           gp = gpar(col = col, lwd = lwd_lines),
                           default.units = "native"
                           )
            grid.points(unit(vals, "native"),
                        unit(seq_along(vals), "native"),
                        pch = pch,
                        size = unit(cex, "char"),
                        gp = gpar(col = col, lwd = lwd_lines),
                        default.units = "native"
                        )
            if (confidence)
                for (i in seq_along(vals)) {
                    ii <- i + jitter
                    grid.lines(unit(c(lwr[offset + i], upr[offset + i]), "native"),
                               unit(c(ii, ii), "native"),
                               gp = gpar(col = col, lwd = lwd_confidence))
                    grid.lines(unit(c(lwr[offset + i], lwr[offset + i]), "native"),
                               unit(c(ii - whiskers/2, ii + whiskers/2), "native"),
                               gp = gpar(col = col, lwd = lwd_confidence))
                    grid.lines(unit(c(upr[offset + i], upr[offset + i]), "native"),
                               unit(c(ii - whiskers/2, ii + whiskers/2), "native"),
                               gp = gpar(col = col, lwd = lwd_confidence))
                }
        }
    } else {
        ## set up plot region
        pushViewport(plotViewport(xscale = xlimaxis, yscale = ylimaxis,
                                  default.units = "native",
                                  name = "oddsratio_plot"))
        grid.xaxis(seq_along(labs), labs)
        grid.yaxis()
        grid.text(xlab, y = unit(-3.5, "lines"))
        grid.text(ylab, x = unit(-3, "lines"), rot = 90)
        grid.text(main, y = unit(1, "npc") + unit(1, "lines"), gp = gp_main)
        pushViewport(viewport(xscale = xlimaxis, yscale = ylimaxis,
                              default.units = "native", clip = "on"))

        ## baseline
        if (baseline)
            grid.lines(unit(c(0,1), "npc"),
                       unit(c(1,1) - LOG, "native"),
                       gp = gp_baseline)

        ## workhorse for one stratum
        draw_one_stratum <- function(vals, pch = "o", col = "black", offset = 0,
                                     jitter = 0) {
            if (bars) {
                if (any(vals > !LOG))
                    grid.rect(unit(seq_along(vals)[vals > !LOG], "native"),
                              unit(vals[vals > !LOG], "native"),
                              width = bar_width,
                              height = unit(vals[vals > !LOG] - !LOG, "native"),
                              just = "top",
                              gp = gp_bars
                              )
                if (any(vals < !LOG))
                    grid.rect(unit(seq_along(vals)[vals < !LOG], "native"),
                              unit(vals[vals < !LOG], "native"),
                              width = bar_width,
                              height = unit(abs(vals[vals < !LOG] - !LOG), "native"),
                              just = "bottom",
                              gp = gp_bars
                              )
            }
            if (lines)
                grid.lines(unit(seq_along(vals), "native"),
                           unit(vals, "native"),
                           gp = gpar(col = col, lwd = lwd_lines),
                           default.units = "native"
                           )
            grid.points(unit(seq_along(vals), "native"),
                        unit(vals, "native"),
                        pch = pch,
                        size = unit(cex, "char"),
                        gp = gpar(col = col, lwd = lwd_lines),
                        default.units = "native"
                        )
            if (confidence)
                for (i in seq_along(vals)) {
                    ii <- i + jitter
                    grid.lines(unit(c(ii, ii), "native"),
                               unit(c(lwr[offset + i], upr[offset + i]), "native"),
                               gp = gpar(col = col, lwd = lwd_confidence))
                    grid.lines(unit(c(ii - whiskers/2, ii + whiskers/2), "native"),
                               unit(c(lwr[offset + i], lwr[offset + i]), "native"),
                               gp = gpar(col = col, lwd = lwd_confidence))
                    grid.lines(unit(c(ii - whiskers/2, ii + whiskers/2), "native"),
                               unit(c(upr[offset + i], upr[offset + i]), "native"),
                               gp = gpar(col = col, lwd = lwd_confidence))
                }
        }

    }

    if (is.null(d))
        draw_one_stratum(values, pch[1], col[1])
    else {
        jitt <- scale(seq_len(prod(d[-1])), scale = 25 * prod(d[-1]))
        for (i in 1 : prod(d[-1]))
            draw_one_stratum(values[(i - 1) * d[1] + seq(d[1])],
                             pch[(i - 1 ) %% length(pch) + 1],
                             col[i],
                             offset = (i - 1) * d[1],
                             jitt[i])
        if (legend)
            grid_legend(legend_pos,
                        labels = apply(expand.grid(dimnames(values)[-1]),
                                       1, paste, collapse = "|"),
                        pch = pch[1 : prod(d[-1])],
                        col = col,
                        lwd = legend_lwd,
                        lty = "solid",
                        size = legend_size,
                        vgap = legend_vgap,
                        gp = gp_legend,
                        gp_frame = gp_legend_frame,
                        inset = legend_inset,
                        title = paste(names(dimnames(values)[-1]), collapse = " x "),
                        gp_title = gp_legend_title, ...)
    }

    grid.rect(gp = gpar(fill = "transparent"))
    if (pop) popViewport(2)
    if (return_grob)
        invisible(grid.grab())
    else
        invisible(NULL)
}

