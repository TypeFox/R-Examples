ancovaplot <- function(object, ...)
  UseMethod("ancovaplot")

ancovaplot.formula <-
  function(object, data, groups=NULL, x=NULL, ...,
           formula=object,
           col=rep(tpg$col,
             length=length(levels(as.factor(groups)))),
           pch=rep(c(15,19,17,18,16,20, 0:14),
             length=length(levels(as.factor(groups)))),
           slope, intercept,
           layout=c(length(levels(cc)), 1),
           col.line=col, lty=1,
           superpose.panel=TRUE,
           between=if (superpose.panel)
                      list(x=c(rep(0, length(levels(cc))-1), 1))
                   else
                      list(x=0),
           ## ignore.groups.panel=FALSE, ## this option is not yet written
           ## ## it will allow an additional panel with a regression line ignoring the groups
           col.by.groups=FALSE ## ignored unless groups= is specified
           ) {
    groups <- substitute(groups)
    x <- substitute(x)
    if (is.call(groups)) groups <- eval(groups)
    if (is.call(x)) x <- eval(x)
    tpg <- trellis.par.get()$superpose.symbol

    ## formula[[1]] ## ~
    ## formula[[2]] ## y
    ## formula[[3]] ## x + g || x || g || x * g
    ## formula[[3]][[1]] ## + *  ## class=="name"
    ## formula[[3]][[2]] ## x ## integer==class(data[[as.character((x[[3]][[2]]))]])
    ## formula[[3]][[3]] ## g ## factor ==class(data[[as.character((x[[3]][[3]]))]])

    formula.xyplot <- formula
    environment(formula.xyplot) <- environment()

    y <- data[[as.character((formula[[2]]))]]

    ## form of equation
    if (class(formula.xyplot[[3]])=="call") {  ## y ~ x + c or y ~ c + x; either with optional groups=
      formula.xyplot[[3]][[1]] <- as.name("|")
      f2 <- data[[as.character(formula[[3]][[2]])]]
      f3 <- data[[as.character(formula[[3]][[3]])]]
      if (is.factor(f2)) {
        cc <- f2 ## cc is condition
        xx <- f3
        x.name <- as.character((formula[[3]][[3]]))
        cc.name <- as.character((formula[[3]][[2]]))
        formula.xyplot[[3]][[2]] <- as.name(x.name)
        formula.xyplot[[3]][[3]] <- as.name(cc.name)

      } else {
        cc <- f3
        xx <- f2
        x.name <- as.character((formula[[3]][[2]]))
      }
      if (is.null(groups)) {
        groups <- cc
##        col.line <- NULL
        groups.cc.incompatible <- FALSE
      }
      else {
        if (is.name(groups))
          groups <- data[[as.character(groups)]]
        else
          if (length(groups) != length(y))
            stop("groups length not equal to y length", call.=FALSE)
        groups.cc.incompatible <- TRUE
      }
      if (!is.null(x) && !(is.name(x) && x == x.name))
        stop("x specified and not equal to covariate name", call.=FALSE)
    }


    if (class(formula[[3]])=="name" &&
        is.numeric(data[[as.character(formula[[3]])]])) { ## y ~ x
      groups.cc.incompatible <- FALSE
      x.name <- as.character(formula[[3]])
      xx <- data[[x.name]]
      if (is.null(groups)) stop("groups must be specified for model y ~ x", call.=FALSE)
      if (is.name(groups)) {
        cc.name <- as.character(groups)
        groups <- data[[cc.name]]
      }
      else {
        if (length(groups) != length(y))
          stop("groups length not equal to y length", call.=FALSE)
        cc.name <- "groups"
      }
      if (!is.null(x) && !(is.name(x) && x == x.name))
        stop("x= specified and not equal to covariate name", call.=FALSE)
      cc <- groups
##      col.line <- NULL
      yx.lm <- lm(y ~ xx)
      if (missing(intercept)) intercept <- coef(yx.lm)[1]
      if (missing(slope)) slope <- coef(yx.lm)[2]
      formula.xyplot[[3]] <-  (~ xx | cc)[[2]]
      formula.xyplot[[3]][[2]] <- as.name(x.name)
    }

    if (class(formula[[3]])=="name" &&
        !is.numeric(data[[as.character(formula[[3]])]])) { ## y ~ c
      cc.name <- as.character(formula[[3]])
      cc <- data[[cc.name]]
      if (is.null(x)) stop("x must be specified for model y ~ groups", call.=FALSE)
      if (is.name(x)) {
        x.name <- as.character(x)
        xx <- data[[x.name]]
      }
      else {
        xx <- x
        if (length(x) != length(y))
          stop("x length not equal to y length", call.=FALSE)
        x.name <- "x"
      }
      if (is.null(groups)) {
        groups <- cc
##        col.line <- NULL
        groups.cc.incompatible <- FALSE
      }
      else
        {
          if (is.name(groups))
            groups <- data[[as.character(groups)]]
          else
            if (length(groups) != length(y))
              stop("groups length not equal to y length", call.=FALSE)
          groups.cc.incompatible <- TRUE
        }
      if (missing(intercept)) intercept <- tapply(y, cc, mean)
      if (missing(slope)) slope <- 0
      formula.xyplot[[3]] <-  (~ xx | cc)[[2]]
      formula.xyplot[[3]][[2]] <- as.name(x.name)
    }

    if (class(formula[[3]])=="call" &&
        as.character(formula[[3]][[1]])=="+") { ## y ~ x + c
      yxc.lm <- aov(y ~ -1 + xx + cc)
      if (missing(intercept)) intercept <- coef(yxc.lm)[-1]
      if (missing(slope)) slope <- coef(yxc.lm)[1]
    }

    if (class(formula[[3]])=="call" &&
        as.character(formula[[3]][[1]])=="*") { ## y ~ x * c
      yxc.coef <- sapply(levels(cc), function(cci) coef(aov(y ~ xx, subset=(cc==cci))))
      if (missing(intercept)) intercept <- yxc.coef[1,]
      if (missing(slope)) slope <- yxc.coef[2,]
    }

    intercept <- rep(intercept, length=length(levels(cc)))
    slope <-  rep(slope, length=length(levels(cc)))

    ## if (is.null(col.line))
    ##   xyplot(formula.xyplot, groups=groups, data=data, ...,
    ##          panel=panel.ancova.superpose,
    ##          panel.groups=panel.ancova.groups,
    ##          slope=slope,
    ##          intercept=intercept,
    ##          layout=layout,
    ##          ...,
    ##          col=col, pch=pch
    ##          )
    ## else {

    groups <- as.factor(groups)
    ## if (is.null(col.line))
    ##   col.line <- col
    ## else
    ##   col.line <- rep(col.line, length=length(levels(cc)))

    if (superpose.panel) {
      data2 <- rbind(
        cbind(data, copy="original"),
        cbind(data, copy="superpose"))

      new.cc.name <- paste(as.character(formula.xyplot[[3]][[3]]), "2", sep=".")
      formula.xyplot[[3]][[3]] <- as.name(new.cc.name)

      data2[[new.cc.name]] <- factor(c(as.character(cc),
                                       rep("Superpose", nrow(data))),
                                     levels=c(levels(cc), "Superpose"))

      data <- data2
      groups <- factor(c(groups, groups), labels=levels(groups)) ## this is new, also check whether cc also needs doubling
      layout[1] <- layout[1]+1
    }

    result <-
      xyplot(formula.xyplot, groups=groups, data=data, ...,
             panel=panel.ancova.superpose,
             panel.groups=panel.points,
             slope=slope,
             intercept=intercept,
             layout=layout,
             col=col, pch=pch,
             col.line=rep(col.line, length=layout[1]),
             lty=lty,
             superpose.panel=superpose.panel,
             col.by.groups=col.by.groups,
             condition.factor=cc,
             groups.cc.incompatible=groups.cc.incompatible,
             between=between
             )
    class(result) <- c("ancovaplot", class(result))
    result
  }




panel.ancova.superpose <-
  function(x, y, subscripts, groups,
           slope, intercept,
           col,
           pch,
           ...,
           col.line, lty,
           superpose.panel,
           col.by.groups,
           condition.factor,
           groups.cc.incompatible,
           plot.resids=FALSE,
           print.resids=FALSE,
           mean.x.line=FALSE,
           col.mean.x.line="gray80") {

    gc <- groups

    if (!col.by.groups && groups.cc.incompatible)
      {
        gc <- interaction(groups, condition.factor)
        col <- rep(col, each=length(levels(groups)))
      }

    if (panel.number() <= length(intercept)) {
      panel.superpose(x, y, subscripts, gc,
                      col=col,
                      pch=pch, ...)

      panel.abline(a=intercept[panel.number()],
                   b=slope[panel.number()],
                   col=col.line[panel.number()],
                   lty=lty)

      if (plot.resids) {
        panel.segments(x, y,
                       x, intercept[panel.number()] + slope[panel.number()]*x,
                       col=col.line[panel.number()])
        if (print.resids)
          cat(panel.number(), y-(intercept[panel.number()] + slope[panel.number()]*x), "\n")
        if (is.numeric(mean.x.line))
          panel.abline(v=mean.x.line, lty=2, col=col.mean.x.line)
      }
    }


    if (superpose.panel && panel.number() > length(intercept)) {

      panel.superpose(x, y, subscripts, gc,
                      col=col,
                      pch=pch, ...)

      for (cci in 1:length(intercept)) {
        panel.abline(a=intercept[cci],
                     b=slope[cci],
                     col=col.line[cci],
                     lty=lty)
        if (plot.resids)
          panel.segments(x[groups==levels(groups)[cci]],
                         y[groups==levels(groups)[cci]],
                         x[groups==levels(groups)[cci]],
                         intercept[cci] + slope[cci]*x[groups==levels(groups)[cci]],
                         col=col.line[cci])

      }

    }
  }
