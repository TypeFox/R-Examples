binreg_plot <-
function(model, main = NULL, xlab = NULL, ylab = NULL,
         xlim = NULL, ylim = NULL, pred_var = NULL, pred_range = c("data", "xlim"),
         group_vars = NULL, base_level = NULL, subset,
         type = c("response", "link"), conf_level = 0.95, delta = FALSE,
         pch = NULL, cex = 0.6, jitter_factor = 0.1,
         lwd = 5, lty = 1, point_size = 0, col_lines = NULL, col_bands = NULL,
         legend = TRUE, legend_pos = NULL, legend_inset = c(0, 0.1),
         legend_vgap = unit(0.5, "lines"),
         labels = FALSE, labels_pos = c("right", "left"),
         labels_just = c("left","center"),
         labels_offset = c(0.01, 0),
         gp_main = gpar(fontface = "bold", fontsize = 14),
         gp_legend_frame = gpar(lwd = 1, col = "black"),
         gp_legend_title = gpar(fontface = "bold"),
         newpage = TRUE, pop = FALSE, return_grob = FALSE)
{
    if (!inherits(model, "glm"))
        stop("Method requires a model of class 'glm'.")
    type <- match.arg(type)
    labels_pos <- match.arg(labels_pos)
    if (is.character(pred_range))
        pred_range <- match.arg(pred_range)

    ## extract data from model
    mod <- model.frame(model)
    term <- terms(mod)
    data.classes <- attr(term, "dataClasses")
    nam <- names(data.classes)

    ## determine response
    r <- attr(term, "response")
    resp <- nam[r]
    data.classes <- data.classes[-r]
    nam <- nam[-r]

    ## determine numeric predictor (take first)
    if (is.null(pred_var)) {
        fac <- data.classes %in% c("factor","logical")
        pred_var_model <- names(data.classes[!fac][1])
        pred_var <-
            names(unlist(sapply(all.vars(term),
                                grep, pred_var_model)))[1]

    } else pred_var_model <- pred_var

    ## filter observed data using model (to account for models fitted with subset=...)
    dat <- model$data[row.names(mod),]
    ## sort observations using order of numeric predictor
    o <- order(dat[,pred_var])
    mod <- mod[o,]
    dat <- dat[o,]

    ## apply subset argument, if any
    if (!missing(subset)) {
        e <- substitute(subset)
        i <- eval(e, dat, parent.frame())
        i <- i & !is.na(i)
        dat <- dat[i,]
        mod <- mod[i,]
    }

    ## determine conditioning variables. Remove all those with only one level observed.
    if (is.null(group_vars)) {
        group_vars <- nam[data.classes %in% "factor"]
        sing <- na.omit(sapply(dat, function(i) all(i == i[1])))
        if (any(sing))
            group_vars <- setdiff(group_vars, names(sing)[sing])
        if(length(group_vars) < 1)
            group_vars <- NULL
    } else
        if (is.na(group_vars) || is.logical(group_vars) && !group_vars[1])
            group_vars <- NULL

    ## set y axis limits - either probability or logit scale
    if(is.null(ylim))
        ylim <- if (type == "response")
                    c(0,1)
                else
                    range(predict(model, dat, type = "link"))
    ## allow for some cosmetic extra space
    ylimaxis <- ylim + c(-1, 1) * diff(ylim) * 0.04

    if(is.null(xlim))
        xlim <- if (is.numeric(pred_range))
                    range(pred_range)
                else
                    range(dat[,pred_var])
    xlimaxis <- xlim + c(-1, 1) * diff(xlim) * 0.04

    ## set default base level ("no effect") of response to first level/0
    if (is.null(base_level))
        base_level <- if(is.matrix(mod[,resp]))
                          2
                      else if(is.factor(mod[,resp]))
                          levels(mod[,resp])[1]
                      else
                          0
    if (is.matrix(mod[,resp]) && is.character(base_level))
        base_level <- switch(base_level, success =, Success = 1, failure =, Failure = 2)

    ## determine labels of conditioning variables, if any
    if (is.null(group_vars)) {
        labels <- legend <- FALSE
    } else {
        ## compute cross-factors for more than two conditioning variables
        if (length(group_vars) > 1) {
            cross <- paste(group_vars, collapse = " x ")
            dat[,cross] <- factor(apply(dat[,group_vars], 1, paste, collapse = " : "))
            group_vars <- cross
        }
        lev <- levels(dat[,group_vars])
    }

    ## set x- and y-lab
    if (is.null(xlab))
        xlab <- pred_var
    if (is.null(ylab))
        ylab <- if (type == "response") {
                    if (is.matrix(mod[,resp]))
                        paste0("P(",
                               c("Failure","Success")[base_level],
                               ")")
                    else
                        paste0("P(", resp, ")")
                } else {
                    if (is.matrix(mod[,resp]))
                        paste0("logit(",
                               c("Failure","Success")[base_level],
                               ")")
                    else
                        paste0("logit(", resp, ")")
                }

    ## rearrange default plot symbol palette
    if (is.null(pch))
        pch <- c(19,15,17, 1:14, 16, 18, 20:25)

    ## determine normal quantile for confidence band
    quantile <- qnorm((1 + conf_level) / 2)

    ## determine default legend position, given the curve's slope
    ## (positive -> topleft, negative -> topright)
    if (is.null(legend_pos))
        legend_pos <-
            if (coef(model)[grep(pred_var, names(coef(model)))[1]] > 0)
                "topleft"
            else
                "topright"

    ## work horse for drawing points, fitted curve and confidence band
    draw <- function(ind, colband, colline, pch, label) {
        ## plot observed data as points on top or bottom
        ycoords <- if (is.matrix(mod[,resp])) {
            tmp <- prop.table(mod[ind,resp], 1)[,switch(base_level, 2, 1)]
            if (type == "link")
                family(model)$linkfun(tmp)
            else
                tmp
        } else
            jitter(ylim[1 + (mod[ind, resp] != base_level)], jitter_factor)
        if (cex > 0)
            grid.points(unit(dat[ind, pred_var], "native"),
                        unit(ycoords, "native"),
                        pch = pch, size = unit(cex, "char"), gp = gpar(col = colline),
                        default.units = "native"
                        )

        ## confidence band and fitted values
        typ <- if (type == "response" && !delta) "link" else type
        if (is.character(pred_range)) {
            if (pred_range == "data") {
                D <- dat[ind,]
                P <- D[,pred_var]
            } else {
                P <- seq(from = xlim[1L], to = xlim[2L], length.out = 100L)
                D <- dat[ind,][rep(1L, length(P)),]
                D[,pred_var] <- P
            }
        } else {
            P <- pred_range
            D <- dat[ind,][rep(1L, length(P)),]
            D[,pred_var] <- P
        }
        pr <- predict(model, D, type = typ, se.fit = TRUE)
        lower <- pr$fit - quantile * pr$se.fit
        upper <- pr$fit + quantile * pr$se.fit
        if (type == "response" && !delta) {
            lower <- family(model)$linkinv(lower)
            upper <- family(model)$linkinv(upper)
            pr$fit <- family(model)$linkinv(pr$fit)
        }
        if (type == "response") { ## cut probs at unit interval
            lower[lower < 0] <- 0
            upper[upper > 1] <- 1
        }
        grid.polygon(unit(c(P, rev(P)), "native"),
                     unit(c(lower, rev(upper)), "native"),
                     gp = gpar(fill = colband, col = NA))
        grid.lines(unit(P, "native"),
                   unit(pr$fit, "native"),
                   gp = gpar(col = colline, lwd = lwd, lty = lty))
        if (point_size > 0)
            grid.points(unit(P, "native"),
                        unit(pr$fit, "native"), pch = pch,
                        size = unit(point_size, "char"),
                        gp = gpar(col = colline))

        ## add labels, if any
        if (labels) {
            x = switch(labels_pos,
                       left = P[1],
                       right = P[length(P)])
            y = switch(labels_pos,
                       left = pr$fit[1],
                       right = pr$fit[length(pr$fit)])
            grid.text(x = unit(x, "native") + unit(labels_offset[1], "npc"),
                      y = unit(y, "native") + unit(labels_offset[2], "npc"),
                      label = label,
                      just = labels_just,
                      gp = gpar(col = colline))
        }
    }

    ## determine colors and plot symbols
    llev <- if (is.null(group_vars)) 1 else length(lev)
    pch <- rep(pch, length.out = llev)
    if (is.null(col_bands))
        col_bands <- colorspace::rainbow_hcl(llev, alpha = 0.2)
    if (is.null(col_lines))
        col_lines <- colorspace::rainbow_hcl(llev, l = 50)

    ## set up plot region, similar to plot.xy()
    if (newpage) grid.newpage()
    pushViewport(plotViewport(xscale = xlimaxis, yscale = ylimaxis, default.units = "native", name = "binreg_plot"))
    grid.rect(gp = gpar(fill = "transparent"))
    grid.xaxis()
    grid.yaxis()
    grid.text(xlab, y = unit(-3.5, "lines"))
    grid.text(ylab, x = unit(-3, "lines"), rot = 90)
    grid.text(main, y = unit(1, "npc") + unit(2, "lines"), gp = gp_main)
    pushViewport(viewport(xscale = xlimaxis, yscale = ylimaxis, default.units = "native", clip = "on"))

    ## draw fitted curve(s)
    if (is.null(group_vars)) {
        ## single curve
        draw(1:nrow(dat),
             col_bands,
             col_lines,
             pch[1])
    } else {
        ## multiple curves
        for (i in seq_along(lev)) {
            ind <- dat[,group_vars] == lev[i]
            draw(ind, col_bands[i], col_lines[i], pch[i], lev[i])
        }

        if (legend)
            grid_legend(legend_pos,
                        labels = lev,
                        col = col_lines,
                        lty = "solid",
                        lwd = lwd,
                        vgap = legend_vgap,
                        gp_frame = gp_legend_frame,
                        inset = legend_inset,
                        title = group_vars,
                        gp_title = gp_legend_title)
    }

    if (pop) popViewport(2)
    if (return_grob)
        invisible(grid.grab())
    else
        invisible(NULL)
}

###########
grid_abline <- function(a, b, ...)
{
    ## taken from graphics::abline()
    if (is.object(a) || is.list(a)) {
        p <- length(coefa <- as.vector(coef(a)))
        if (p > 2)
            warning(gettextf("only using the first two of %d regression coefficients",
                p), domain = NA)
        islm <- inherits(a, "lm")
        noInt <- if (islm)
            !as.logical(attr(stats::terms(a), "intercept"))
        else p == 1
        if (noInt) {
            a <- 0
            b <- coefa[1L]
        }
        else {
            a <- coefa[1L]
            b <- if (p >= 2)
                coefa[2L]
            else 0
        }
    }
    grid.abline(a, b, ...)
}
