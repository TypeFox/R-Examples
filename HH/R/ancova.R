"ancova" <- function(formula, data.in=NULL,
                     ..., x, groups,
                     transpose=FALSE,
                     display.plot.command=FALSE,
                     superpose.level.name="superpose",
                     ignore.groups=FALSE, ignore.groups.name="ignore.groups",
                     blocks, blocks.pch=letters[seq(levels(blocks))],
                     layout, between, main,
                     pch=trellis.par.get()$superpose.symbol$pch) {
  ## on.exit(browser())
  if (missing(data.in)) stop("Please explicitly name a data.frame.")
  a.aov <- aov(formula, data=data.in)
  a.aov$call$formula <- substitute(formula)
  a.aov$call$data <- substitute(data.in)

  ## determine sequence:  x+a or x*a  gives c(FALSE,TRUE)
  ##                      a+x or a*x  gives c(TRUE,FALSE)
  ##                      a, x=x      gives c(TRUE)
  ##                      x, groups=a gives c(FALSE)
  ## anything else is an error
  tl <- attr(a.aov$terms,"term.labels")
  if (length(tl)==3) tl <- tl[1:2]
  classes <- sapply(data.in[,tl, drop=FALSE], is.factor)
                  ## sapply(tl,
                  ##   function(cc, data) {
                  ##     ccc <- class(data[,cc])
                  ##     if (is.null(ccc)) FALSE
                  ##     else ccc=="factor"},
                  ##   data.in)
  if ((length(classes)==2 && classes[1]==classes[2])## two factors or two numeric
      ||
      (length(classes)==1 && !classes  ## group as numeric
       && ((!missing(groups) &&
           !is.factor(data.in[,deparse(substitute(groups))]))
           ||
           (!missing(x))
           )
       )
      ||
      (length(classes)==1 && classes   ## x as factor
       && ((!missing(x) && is.factor(data.in[,deparse(substitute(x))]))
           ||
           (!missing(groups))
           ))
      )
    stop("ancova requires exactly one factor and exactly one numeric variable.")

  formula.plot <- formula

  if (length(formula[[3]]) == 3) { ## (y ~ x | a) or (y ~ x | a)
    formula.plot[[3]][[1]] <- as.name("|")
    formula.plot[[3]][2:3] <- formula[[3]][2+classes] ## y ~ x | a
    coef.aov <- coef(a.aov)
  }
  else { ## (y ~ a, x=x) or (y ~ x, groups=a)
    formula.plot[[3]] <- (~  x | g)[[2]]
    if (!missing(x)) { ## (y ~ a, x=x)
#      formula.plot[[3]][[2]] <- as.name(deparse(substitute(x))) ## x
      formula.plot[[3]][[2]] <- substitute(x) ## x
      formula.plot[[3]][[3]] <- formula[[3]]                    ## a
      classes <- c(FALSE, classes)
      tl <- c(deparse(substitute(x)), tl)
      coef.aov <- c(coef(a.aov)[1], x=0, coef(a.aov)[-1])
    }
    if (!missing(groups)) { ## y ~ x, groups=a)
      formula.plot[[3]][[2]] <- formula[[3]]                         ## x
#      formula.plot[[3]][[3]] <- as.name(deparse(substitute(groups))) ## a
      formula.plot[[3]][[3]] <- substitute(groups) ## a
      classes <- c(classes, TRUE)
      tl <- c(tl, deparse(substitute(groups)))
      coef.aov <-
        c(coef(a.aov),
          rep(0, length(levels(data.in[[deparse(substitute(groups))]]))-1))
    }
    if (missing(groups) == missing(x)) stop("Invalid formula")
  }

  ## xyplot(formula.plot, data=data.in, ...) ## constructed
  m <- match.call()
  m[[1]] <- as.name("xyplot")
  m$formula <- formula.plot
  names(m)[3] <- "data"
  if (length(formula[[3]]) == 1) m$x <- NULL
  m$display.plot.command <- NULL
  m$main <- list(label=deparse(substitute(formula)))

  if (!missing(main)) {
    if (is.list(main)) {
      label.location <- which(names(main) == "")
      if (length(label.location)==1)
        names(main)[label.location] <- "label"
    }
    if (is.atomic(main)) main <- list(label=main)
    m$main[names(main)] <- main
  }
  m$coef <- coef.aov
  if (missing(groups)) {
    m$groups <- data.in[[tl[classes]]]
  }
  else {
    m$groups <- data.in[[as.character(m$groups)]]
  }
  m$contrasts <- contrasts(m$groups)
  m$classes <- classes
  m$panel <- "panel.ancova"
  a.labels <- dimnames(m$contrasts)[[1]]

  tpgs <- trellis.par.get("superpose.symbol")

  tpgs$pch[] <- pch
  if (is.null(m$pch)) m$pch <- pch
  tpgl <- trellis.par.get("superpose.line")

  m$key <- list(text=list(a.labels),   ## treatment key
                points = Rows(tpgs, 1:length(a.labels)),
                lines = Rows(tpgl, 1:length(a.labels)),
                border=TRUE,
                space="right",
                title=as.character(formula.plot[[3]][[3]]))

  if (!missing(blocks)) {
    blocks.char <- deparse(substitute(blocks))
    bl <- data.in[[blocks.char]]
    if (!is.null(bl)) {
      blocks <- bl
      m$blocks <- blocks
    }
    m$blocks.pch <- blocks.pch
    m$blocks.cex <- m$key$points$cex
    m$key$points <- NULL
    m$sub <- paste(blocks.char, ": ",
                   paste(blocks.pch, collapse=" "),
                   sep="")
  }
  levels.a <- levels(data.in[[as.character(formula.plot[[3]][[3]])]])
  data.in[[as.character(formula.plot[[3]][[3]])]] <-
    factor(data.in[[as.character(formula.plot[[3]][[3]])]],
           levels=c(levels.a, superpose.level.name,
             if(ignore.groups) ignore.groups.name))
## create a groups column that duplicates the classes column
  data.in <- cbind(data.in, new.groups=data.in[[tl[classes]]])
  data2 <- data.in
  data2[[tl[classes]]] <- superpose.level.name
  data3 <- rbind(data.in, data2)
  if (ignore.groups) {
    data2[[tl[classes]]] <- ignore.groups.name
    data3 <- rbind(data3, data2)
  }
  m$data <- data3
  m$groups <- data3$new.groups
  if (missing(layout))
    m$layout <- c(length(levels.a)+1, 1)
  else m$layout <- layout
  if (missing(between))
    m$between <- list(x=c(rep(0,length(levels.a)-1),2), y=0)
  else m$between <- between
  m$superpose.level.name <- NULL

  ## print or evaluate the xyplot call
  if (display.plot.command) print(m)
  m$transpose <- NULL
  if.R(r={names(m)[[match("formula", names(m))]] <- "x"},
       s={})
  if (transpose)
    attr(a.aov,"trellis") <- t(eval(m))
  else
    attr(a.aov,"trellis") <- eval(m)
  oldClass(a.aov) <- c("ancova", oldClass(a.aov))
  a.aov
}

"anova.ancova" <-
function(object, ...)
  NextMethod("anova")

"predict.ancova" <-
function(object, ...)
  NextMethod("predict")

"print.ancova" <-
function(x, ...) {
  print(anova(x, ...))
  print(attr(x,"trellis"))
  invisible(x)
}

"model.frame.ancova" <-
  if.R(r={function(formula, ...)
            NextMethod("model.frame")},
       s={function(formula, data = NULL, na.action = na.fail, ...)
            NextMethod("model.frame")})

"summary.ancova" <-
function(object, ...)
  NextMethod("summary")

"plot.ancova" <-
  if.R(r={
    function(x, y, ...) {
      x.full <- x
      attr(x, "trellis") <- NULL
      NextMethod("plot")
      invisible(x.full)
    }
  },s={
    function(x, ...) {
      x.full <- x
      attr(x, "trellis") <- NULL
      NextMethod("plot")
      invisible(x.full)
    }
  })

"coef.ancova" <-
function(object, ...)
  NextMethod("coef")

## "coefficients.ancova" <-
## function(object, ...)
##   NextMethod("coef")

"panel.ancova" <-
function(x, y, subscripts, groups, transpose=FALSE, ...,
                         coef, contrasts, classes, ignore.groups,
                         blocks, blocks.pch, blocks.cex, pch) {
##  contrasts <- contrasts[-nrow(contrasts), -ncol(contrasts)]
  n.contr <- ncol(contrasts)
  if (length(classes)==1)
    coefs.a <- (1:n.contr) + 1
  else
    coefs.a <- (1:n.contr) + (1:2)[classes]
  a <- coef[1] + contrasts %*% coef[coefs.a]

  if (length(classes)==1)
    b <- rep(0, length(a))
  else {
    b <- coef[2+c(0,n.contr)[!classes]]
    if (n.contr != (length(coef)-2))
      b <- b + contrasts %*% coef[-(1:(n.contr+2))]
    else
      b <- rep(b, length(a))
  }

  if (transpose) {
    a.untransposed <- a
    a <- ifelse (b==0, 0, -a/b)
    b <- 1/b  ## if (b==0) Inf
  }

  tpgs <- trellis.par.get("superpose.symbol")
  tpgs$pch[] <- pch
  tpgl <- trellis.par.get("superpose.line")

  ## browser()
  cell <-
    ## if.R(r=
         panel.number()
         ## ,
         ##       s=get("n", frame=sys.parent()))

  if (cell == length(a)+1) {

    if (missing(blocks))
      panel.superpose(x, y, subscripts=subscripts, groups=groups, ..., pch=pch)
    else
      panel.superpose(x, y, subscripts=subscripts, groups=rep(blocks,2),
                        pch=blocks.pch, cex=blocks.cex, ...)
    for (i in seq(along=a)) {
      if (abs(b[i]) == Inf)
        panel.abline(v=a.untransposed[i],
                     col=tpgl$col[i], lty=tpgl$lty[i], lwd=tpgl$lwd[i])
      else
        panel.abline(a=a[i], b=b[i],
                     col=tpgl$col[i], lty=tpgl$lty[i], lwd=tpgl$lwd[i])
    }
  }
  else
    if (cell == length(a)+2) {
      if (missing(blocks))
        panel.superpose(x, y, subscripts=subscripts, groups=groups, ..., pch=pch)
      else
        panel.superpose(x, y, subscripts=subscripts, groups=rep(blocks,3),
                        pch=blocks.pch, cex=blocks.cex, ...)
      if (transpose) {
        tmp.lm <- coef(lm(x ~ y))       # transpose=TRUE interchanged the names
                                        # we must un-transpose here
        ## isolate coefficients
        a <- tmp.lm["(Intercept)"]
        b <- tmp.lm["y"]
        ## transpose coefficients
        a[] <- ifelse (b==0, 0, -a/b)   # keep "(Intercept)" name
        b <- 1/b  ## if (b==0) Inf
        names(b) <- "x"                 # keep "x" name
      }
      else {
        tmp.lm <- coef(lm(y ~ x))
        a <- tmp.lm["(Intercept)"]
        b <- tmp.lm["x"]
      }
      if (transpose && tmp.lm["y"]==0)
        panel.abline(v=tmp.lm["(Intercept)"])
      else
        panel.abline(a, b)
    }
    else {
      if (missing(blocks))
          panel.xyplot(x, y, col=tpgl$col[cell], pch=tpgs$pch[cell], ...)
        else
          panel.superpose(x, y, subscripts=subscripts, groups=blocks,
                        pch=blocks.pch, cex=blocks.cex, ...)
      if (abs(b[cell]) == Inf)
        panel.abline(v=a.untransposed[cell],
                     col=tpgl$col[cell],
                     lty=tpgl$lty[cell],
                     lwd=tpgl$lwd[cell])
      else
        panel.abline(a=a[cell], b=b[cell],
                     col=tpgl$col[cell],
                     lty=tpgl$lty[cell],
                     lwd=tpgl$lwd[cell])
    }
}

setOldClass(c("ancova", "aov", "lm"))
