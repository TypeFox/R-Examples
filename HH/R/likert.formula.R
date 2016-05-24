plot.likert.formula <- function(x, data, ReferenceZero=NULL, value, levelsName="",

                                scales.in=NULL,   ## use scales=
                                between=list(x=1 + (horizontal), y=.5 + 2*(!horizontal)),
                                auto.key.in=NULL, ## use auto.key=
                                panel.in=NULL,    ## use panel=
                                horizontal=TRUE,
                                par.settings.in=NULL, ## use par.settings=
                                ...,
                                as.percent = FALSE,
                                ## titles
                                ylab= if (horizontal) {
                                  if (length(x)==3)
                                    deparse(x[[2]])
                                  else
                                    "Question"
                                }
                                else
                                if (as.percent != FALSE) "Percent" else "Count",

                                xlab= if (!horizontal) {
                                  if (length(x)==3)
                                    deparse(x[[2]])
                                  else
                                    "Question"
                                }
                                else
                                if (as.percent != FALSE) "Percent" else "Count",

                                main = x.sys.call,

                                ## right axis
                                rightAxisLabels = rowSums(data.list$Nums),
                                rightAxis = !missing(rightAxisLabels),
                                ylab.right = if (rightAxis) "Row Count Totals" else NULL,
                                xlab.top = NULL,
                                right.text.cex =
                                  if (horizontal) { ## lazy evaluation
                                    if (!is.null(scales$y$cex)) scales$y$cex else .8
                                  }
                                  else
                                    {
                                      if (!is.null(scales$x$cex)) scales$x$cex else .8
                                    },

                                ## scales
                                xscale.components = xscale.components.top.HH,
                                yscale.components = yscale.components.right.HH,
                                xlimEqualLeftRight = FALSE,
                                xTickLabelsPositive = TRUE,

                                ## row sequencing
                                as.table=TRUE,
                                positive.order=FALSE,
                                data.order=FALSE,
                                reverse=ifelse(horizontal, as.table, FALSE),

                                ## resizePanels arguments
                                h.resizePanels=sapply(result$y.used.at, length),
                                w.resizePanels=sapply(result$x.used.at, length),

                                ## color options
                                reference.line.col="gray65",
                                key.border.white=TRUE,
                                col=likertColor(Nums.attr$nlevels,
                                  ReferenceZero=ReferenceZero,
                                  colorFunction=colorFunction,
                                  colorFunctionOption=colorFunctionOption),
                                colorFunction="diverge_hcl",
                                colorFunctionOption="lighter"
                                ) {
  ## force(rightAxis)
  rightAxisMissing <- missing(rightAxis)  ## needed by as.percent
  if (positive.order) data.order <- FALSE  ## both TRUE is not meaningful.  positive.order wins.

if (!missing(value)) {
  x.sys.call <- deparse(match.call()[1:4], width.cutoff = 500L)
  varNamesUsedLong <- c(getVarNames(x, data), list(Value=value)) ## list(Question, Conditions, Nums, Value)
  levelsName <- varNamesUsedLong$LevelNames[[1]]
  data.list.list <- getLikertDataLong(x, data, varNamesUsedLong)
  data.list <- data.list.list$data.list  ## reshaped
  varNamesUsed <- data.list.list$varNamesUsed ## Nums names to match reshaped data, Value no longer needed
  x <- data.list.list$x ## or FormulaString
} else {
  x.sys.call <- deparse(match.call()[1:3], width.cutoff = 500L)

  varNamesUsed <- getVarNames(x, data)
  ## list(QuestionName, CondNames, LevelNames) ## Subset of data columns actually used

  data.list <- getLikertData(data, varNamesUsed) ## list(Question, Conditions, Nums)
}


  if (as.percent != FALSE) {
    Nums.pct <- data.list$Nums / rowSums(data.list$Nums) * 100
    Nums.pct[data.list$Nums == 0] <- 0
    Nums.lik <- as.likert(Nums.pct, ReferenceZero=ReferenceZero)
    if (rightAxisMissing && as.percent != "noRightAxis" ) {
      rightAxis <- TRUE
      if (missing(ylab.right))
        ylab.right <- "Row Count Totals"
    }
    ## else
    ##  rightAxis <- FALSE
  } else {
    Nums.lik <- as.likert(data.list$Nums, ReferenceZero=ReferenceZero)
  }

  par.settings <- par.settings.in
  if (rightAxis)  {
    par.settings$clip$panel <- "off"
    if (horizontal) {
      par.settings$layout.widths$ylab.right <-
        max(6, par.settings$layout.widths$ylab.right, na.rm=TRUE)
      par.settings$layout.widths$axis.key.padding <-
        max(2, par.settings$layout.widths$axis.key.padding, na.rm=TRUE)
    }
    else { ## vertical
      ##par.settings$layout.heights$xlab.top <-
      ##  max(8, par.settings$layout.heights$xlab.top, na.rm=TRUE)
     par.settings$layout.heights$main.key.padding <- max(2, par.settings$layout.heights$main.key.padding, na.rm=TRUE)
     par.settings$layout.heights$key.axis.padding <- max(1.5, par.settings$layout.heights$key.axis.padding, na.rm=TRUE)
 ##    par.settings$layout.heights$axis.top <-
 ##      max(2, par.settings$layout.heights$axis.top, na.rm=TRUE)
     }
  }

  Nums.attr <- attributes(Nums.lik)
## recover()
  scales <- list(x=list(alternating=1), y=list(alternating=1))
  if (!missing(scales.in)) {
    scales.x <- scales$x
    scales.y <- scales$y
    if (is.null(scales.in$x) && is.character(scales.in$y))
      scales.in$y <- list(relation=scales.in$y)
    if (is.null(scales.in$y) && is.character(scales.in$x))
      scales.in$x <- list(relation=scales.in$x)
    if (is.character(scales.in))
      scales.in <- list(x=list(relation=scales.in$x), y=list(relation=scales.in$y))
    scales.x[names(scales.in$x)] <- scales.in$x
    scales.y[names(scales.in$y)] <- scales.in$y
    scales[names(scales.in)] <- scales.in
    scales$x <- scales.x
    scales$y <- scales.y
  }

  lim <- NULL
  if (horizontal) {
    xlim <- list(...)$xlim
    if (!is.null(xlim) && is.null(scales.in$x$limits))
      lim <- xlim
    if (!is.null(scales.in$x$limits))
      lim <- scales.in$x$limits
  }
  else {
    ylim <- list(...)$ylim
    if (!is.null(ylim) && is.null(scales.in$y$limits))
      lim <- ylim
    if (!is.null(scales.in$y$limits))
      lim <- scales.in$y$limits
  }

  if (((horizontal && is.null(scales$x$at) && is.null(scales$x$labels)) ||
       (!horizontal && is.null(scales$y$at) && is.null(scales$y$labels)))
      &&
      (xlimEqualLeftRight || xTickLabelsPositive)) {

    if (is.null(lim)) {
      tmp <- Nums.lik
      tmp[Nums.lik < 0] <- 0
      data.max <- max(rowSums(tmp))
      tmp <- Nums.lik
      tmp[Nums.lik > 0] <- 0
      data.min <- min(rowSums(tmp))
      lim <- c(data.min, data.max)
      lim <- lim + c(-.04, .04)*diff(lim)
    }

    if (xlimEqualLeftRight)
      lim <- c(-1, 1)*max(abs(lim))

    at <- pretty(lim)

    if (horizontal && !is.null(scales.in$x$at))
      at <- scales.in$x$at
    if (!horizontal && !is.null(scales.in$y$at))
      at <- scales.in$y$at

    if (xTickLabelsPositive)
      labels <- abs(at)
    else
      labels <- at

    if (horizontal && !is.null(scales.in$x$labels))
      at <- scales.in$x$labels
    if (!horizontal && !is.null(scales.in$y$labels))
      at <- scales.in$y$labels


    if (horizontal) {
      scales$x$limits <- lim
      scales$x$at <- at
      scales$x$labels <- labels
    }
    else {
      scales$y$limits <- lim
      scales$y$at <- at
      scales$y$labels <- labels
    }

  }


  ## between <- list(x=1, y=.5)
  ## if (!missing(between.in)) between[names(between.in)] <- between.in

  if (horizontal)
    FormulaString <- with(varNamesUsed,
                        paste("`", QuestionName, "` ~ .value",
                              if (is.null(CondNames))
                              NULL
                              else
                              paste(" |", paste(CondNames, collapse=" + ")),
                              sep=""))
  else
     FormulaString <- with(varNamesUsed,
                        paste(".value ~ `", QuestionName, "`",
                              if (is.null(CondNames))
                              NULL
                              else
                              paste(" |", paste(CondNames, collapse=" + ")),
                              sep=""))

  if (is.logical(auto.key.in) && length(auto.key.in) == 1 && auto.key.in == FALSE)
    auto.key <- FALSE
  else
    {
      auto.key=list(
        title=levelsName,
        text=Nums.attr$original.levels,
        columns=ifelse(horizontal, Nums.attr$nlevels, 1),
        space=ifelse(horizontal, "bottom", "right"),
        reverse.rows=ifelse(horizontal, FALSE, TRUE),
        size=2, cex=.8,
        between=0.6,
        points=FALSE,
        rectangles=FALSE,
        rect=list(
          col=col,
          border=if (key.border.white) "white" else col
          ))

      if (!missing(auto.key.in)) auto.key[names(auto.key.in)] <- auto.key.in
    }

  data2 <- with(data.list,
                data.frame(rightAxisLabels=rightAxisLabels, Question, Conditions, Nums.lik, check.names=FALSE))

  names(rightAxisLabels) <- data.list$Question[[1]]

  {

    if (positive.order || data.order) {
      if (positive.order) {
        if (reverse)
          data2 <- data2[Nums.attr$positive.order, ]      ## 1 1 0 and 1 1 1
        else
          data2 <- data2[rev(Nums.attr$positive.order), ] ## 1 0 0 and 1 0 1
      } else { ## data.order
        do <- 1:nrow(data2)
        if (reverse)
          data2 <- data2[do, ]                            ## 0 1 1
        else
          data2 <- data2[rev(do), ]                       ## 0 0 1
      }
      newQ <- factor(data2[[varNamesUsed$QuestionName]],
                     levels=unique(data2[[varNamesUsed$QuestionName]]))

      new.order <- do.call(order,
                           data.frame(data2[, varNamesUsed$CondNames, drop=FALSE],
                                      newQ=newQ,
                                      check.names=FALSE))

      data2[[varNamesUsed$QuestionName]] <-
        factor(data2[[varNamesUsed$QuestionName]], levels=unique(newQ[new.order]))
    }


    if (!positive.order && reverse && !data.order)
      data2[[varNamesUsed$QuestionName]] <-               ## 0 1 0
        factor(data2[[varNamesUsed$QuestionName]],
               levels=rev(levels(data2[[varNamesUsed$QuestionName]])))

    ## if (!positive.order && !reverse && !data.order)    ## default
    ##   ## data2 <- data2                                ## 0 0 0


    ##                 positive.order  reverse  data.order
    ## previous        1               1        0
    ##                 1               0        0
    ##                 0               1        0
    ##                 0               0        0

    ## new             0               1        1
    ##                 0               0        1

    ## ==> 1 1 0       1               1        1
    ## ==> 1 0 0       1               0        1
   }

  ## data2.melt <- melt(data2[,-1],
  ##                    id.vars=unlist(varNamesUsed[1:2]))
  #### data2.melt <- reshape2::melt((data2[match(unique(names(data2)), names(data2))])[,-1],
  ####                    id.vars=unique(unlist(varNamesUsed[1:2])),
  ####                    variable_name=".variable")
  if (rightAxis)
    data2.melt <- reshape2::melt((data2[match(unique(names(data2)), names(data2))]),
                                id.vars=c(unique(unlist(varNamesUsed[1:2])), "rightAxisLabels"),
                                variable.name=".variable")
  else
    data2.melt <- reshape2::melt((data2[match(unique(names(data2)), names(data2))])[,-1],
                                id.vars=unique(unlist(varNamesUsed[1:2])),
                                variable.name=".variable")
  names(data2.melt)[ncol(data2.melt)] <- ".value"  ## avoid name conflict with "value"

  panel <- function(x, y, subscripts, ..., horizontal=horizontal,
                    rightAxis=rightAxis, rightAxisLabels=rightAxisLabels,
                    reference.line.col=reference.line.col,
                    right.text.cex=.8) {
    if (horizontal)
      panel.abline(v=0, col=reference.line.col)
    else
      panel.abline(h=0, col=reference.line.col)
    panel.barchart(x, y, subscripts=subscripts, ..., horizontal=horizontal)

    if (rightAxis) {
      if (horizontal) {
        at.which <- match(levels(y), y)
        labels <- (rightAxisLabels[subscripts])[at.which]
        panel.axis("right",
                   at=seq(along=levels(y)), labels=labels,
                   outside=TRUE, half=FALSE,
                   text.cex=right.text.cex)
      } else {
        at.which <- match(levels(x), x)
        labels <- (rightAxisLabels[subscripts])[at.which]
        panel.axis("top",
                   at=seq(along=levels(x)), labels=paste(labels, "\n", sep=""),
                   outside=TRUE, half=FALSE, rot=0,
                   text.cex=right.text.cex)
      }
    }

  }
  if (!is.null(panel.in)) panel <- panel.in

  ## if (!horizontal) {
  ##   tmp <- xlab.top
  ##   xlab.top <- ylab.right
  ##   ylab.right <- tmp
  ## }
  result <-
  barchart(as.formula(FormulaString), groups=data2.melt$.variable,
           data=data2.melt,
           as.table=as.table,
           xlab=xlab, ylab=ylab,
           ylab.right=ylab.right,
           xlab.top=xlab.top,
           main=main, horizontal=horizontal,
           stack=TRUE,
           reference=TRUE,
           col=col[Nums.attr$color.seq],
           border=col[Nums.attr$color.seq],
           panel=panel,
           scales=scales,
           right.text.cex=right.text.cex,
           between=between,
           auto.key=auto.key,
           par.settings=par.settings,
           ...,
           xscale.components=xscale.components,
           yscale.components=yscale.components,
           rightAxis=rightAxis,
           rightAxisLabels=data2.melt$rightAxisLabels,
           subscripts=TRUE
           )


  if (length(h.resizePanels) > 0) {
    if (is.character(h.resizePanels) && h.resizePanels=="rowSums")
      h.resizePanels <- rowSums(data.list$Nums)
    result <- resizePanels(result, h=h.resizePanels)
  }

  if (length(w.resizePanels) > 0) {
    if (is.character(w.resizePanels) && w.resizePanels=="rowSums")
      w.resizePanels <- rowSums(data.list$Nums)
    result <- resizePanels(result, w=w.resizePanels)
  }

  result
}



getVarNames <- function(x, data) {
  switch(length(x),
         stop(paste("formula must be in form",
                    "  Question ~ . | Subtable",
                    "or",
                    "  Question ~ A+B+C | Subtable",
                    "or",
                    "  Question ~ .",
                    "or",
                    "  Question ~ A+B+C",
                    "or",
                    "  ~ . | Subtable",
                    "or",
                    "  ~ A+B+C | Subtable",
                    "or",
                    "  ~ .",
                    "or",
                    "  ~ A+B+C",
                    sep="\n"),
              call.=FALSE),
         {
           QuestionName <- 'rownames(data)'
           x <- as.formula(paste('` `', deparse(x, width.cutoff = 500L))) ## length(x) is now 3
         },
         QuestionName <- as.character(x[[2]]))

  x3 <- x[[3]]
  DOT <- as.call(~ .)[[2]]

  if (x3 == DOT) {
    CondNamesFormula <- NULL
    LevelNamesFormula <- NULL
  }
  else { ## ((x3 != DOT)
    if (length(x3) > 1 && x3[[1]] == '|' && x3[[2]] == DOT) {
      CondNamesFormula <- x3[[3]]
      LevelNamesFormula <- NULL
    }
    else {
      if (length(x3) > 1 && x3[[1]] == '|' && x3[[2]] != DOT) {
        LevelNamesFormula <- x3[[2]]
        CondNamesFormula <- x3[[3]]
      }
      else { ## length(x3) == 1 || x3[[1]] != '|'
        LevelNamesFormula <- x3
        CondNamesFormula <- NULL
      }
    }
  }

  if (is.null(CondNamesFormula))
    CondNames <- NULL
  else {
    CondNamesRaw <- strsplit(deparse(CondNamesFormula, width.cutoff = 500L), "+", fixed=TRUE)[[1]]
    CondNames <- sub('^[[:space:]]+', '', sub('[[:space:]]+$', '', CondNamesRaw )) ## remove leading and trailing white space, POSIX-style
    CondNames <- sub('^\"', '', sub('\"$', '', CondNames )) ## remove leading and trailing '\"' character
  }

  if (is.null(LevelNamesFormula))
    LevelNames <- NULL
  else {
    LevelNamesRaw <- strsplit(deparse(LevelNamesFormula, width.cutoff = 500L), "+", fixed=TRUE)[[1]]
    LevelNames <- sub('^[[:space:]]+', '', sub('[[:space:]]+$', '', LevelNamesRaw )) ## remove leading and trailing white space, POSIX-style
    LevelNames <- sub('^\"', '', sub('\"$', '', LevelNames )) ## remove leading and trailing '\"' character
  }

  list(QuestionName=QuestionName, CondNames=CondNames, LevelNames=LevelNames)
}


getLikertData <- function(data, varNamesUsed) {
  if (varNamesUsed$QuestionName == 'rownames(data)') {
    Question <- data.frame('rownames(data)'=factor(rownames(data), levels=unique(rownames(data))),
                           check.names=FALSE)
    RemainingNames <- names(data)
  }
  else {
    Question <- data[, varNamesUsed$QuestionName, drop=FALSE]
    if (!is.factor(Question[[1]]))
      Question[[1]] <- factor(Question[[1]], levels=unique(Question[[1]]))
    RemainingNames <- names(data)[! names(data) %in% varNamesUsed$QuestionName]
  }

  if (length(varNamesUsed$CondNames) > 0) {
    Conditions <- data.frame(
      lapply(data[, varNamesUsed$CondNames, drop=FALSE],
             function(ff) {
               if (is.factor(ff))
                 ff
               else
                 factor(ff, levels=unique(ff))
             }))
    RemainingNames <- RemainingNames[! RemainingNames %in% varNamesUsed$CondNames]
  }
  else
    Conditions <- data.frame(matrix(NA, nrow(Question), 0))

  Nums <- if (is.null(varNamesUsed$LevelNames)) {
    tmp <- data[, RemainingNames, drop=FALSE]
    tmp[, sapply(tmp, is.numeric)]
  }
  else
    data[, varNamesUsed$LevelNames, drop=FALSE]

  Nums <- data.matrix(Nums)
  names(dimnames(Nums)) <- if (is.null(attr(data, "names.dimnames")))
    c(varNamesUsed$QuestionName, "Level Names")
  else
    attr(data, "names.dimnames")
  list(Question=Question, Conditions=Conditions, Nums=Nums)

}

getLikertDataLong <- function(x, data, varNamesUsedLong) {
  if (class(x[[3]]) == "call") {
    cond <- deparse(x[[3]][[3]])
    aaa <- strsplit(cond, " ", fixed=TRUE)[[1]]
    aaa[aaa=='+'] <- ' + '
    bbb <- paste(rev(aaa), collapse="")
    y <- paste(bbb, ' + ', as.character(x[[2]]), ' ~ ', x[[3]][[2]], sep='')
    levelsName <- as.character(x[[3]][[2]])
    x[[3]][[2]] <- as.name(".")
  } else {
    y <- x
    levelsName <- as.character(x[[3]])
    x[[3]] <- as.name(".")
  }

  varNamesUsed <- c(varNamesUsedLong[c("QuestionName","CondNames")],
                    Nums=levels(data[[varNamesUsedLong$LevelNames]]))
  data2 <- reshape2::dcast(y, data=data[unlist(varNamesUsedLong)], value.var=varNamesUsedLong$Value)
  list(data.list=getLikertData(data2, varNamesUsed), varNamesUsed=varNamesUsed, x=x)
}

likertStripDefault <- strip.custom(
  bg="gray97", ## col.strip.background,
  par.strip.text=list(cex=.6))

## if (FALSE) {

##   require(HH)
##   require(reshape)

##   data(SFF8121)
##   SFF <- as.likertDataFrame(SFF8121)

##   ## These two statements are equivalent to each other and will be
##   ## equivalent to the next two in the special case that there are no
##   ## repetitions in the rownames(data).

##   likert( ~ . | Subtable, data=SFF, layout=c(2,1))

##   likert( ~
##          "Strongly Disagree" + Disagree + Neutral + Agree + "Strongly Agree" | Subtable,
##          data=SFF, layout=c(2,1))

##   ## These two statements are equivalent
##   likert(Question ~ . | Subtable, data=SFF, layout=c(2,1))

##   likert(Question ~
##          "Strongly Disagree" + Disagree + Neutral + Agree + "Strongly Agree" | Subtable,
##          data=SFF, layout=c(2,1))

##   ## fancy
##   fancy <-
##   update(par.strip.text=list(cex=.7), xlab="Percent", main="Student Evaluations", scales=list(x=list(limits=c(-20, 100))),
##          resizePanels(useOuterStrips(combineLimits(
##            likert(Question ~ . | Subtable+fake,
##                   data=cbind(SFF,
##                     fake=rep(
##                       factor(c("Instructor","Course","Instructor","Course"),
##                              levels=c("Instructor","Course")),
##                       c(8,3,8,3))),
##                   scales=list(y=list(relation="free")), between=list(y=.5, x=.5))
##            )), h=c(8,3))
##          )
##   fancy

##   ## others
##   tmpvert <-
##   likert(Question ~ . | Subtable, data=SFF, ylab="Percent", layout=c(1,2), horizontal=FALSE,
##          between=list(y=1), scales=list(x=list(rot=90), y=list(limits=c(-30,110), alternating=1)))
##   tmpvert

##   likert(Question ~ . | Subtable, data=SFF, ylab="Percent", layout=c(1,2), horizontal=FALSE,
##          between=list(y=1), scales=list(x=list(rot=90)), ylim=c(-30,110))

##   tmp <- likert(Question ~ . | Subtable, data=SFF, layout=c(2,1),
##                 scales=list(x=list(limits=c(-30, 110), alternating=1, at=seq(-20, 100, 20), labels=abs(seq(-20, 100,20)))))
##   tmp


##   tmp2 <- likert(Question ~ . | Subtable, data=SFF, layout=c(2,1),
##                  scales=list(x=list(limits=c(-30, 110), alternating=1)))
##   tmp2

##   ## ## if (!is.numeric(tmp2$x.scales$labels)) {
##   ## ##   if (!is.numeric(tmp2$x.scales$at))
##   ## ##     tmp2$x.scales$at <- pretty(tmp2$x.scales$limits)
##   ## ##   tmp2$x.scales$labels <- abs(tmp2$x.scales$at)
##   ## ## }
##   ## ## tmp2

##   ## TickLabelsPositive <- function(trellis, horizontal=TRUE) { ## trellis is a trellis object
##   ##   scales <-if (horizontal) trellis$x.scales else trellis$y.scales
##   ##   if (!is.numeric(scales$labels)) {
##   ##     if (!is.numeric(scales$at))
##   ##       scales$at <- pretty(scales$limits)
##   ##     scales$labels <- abs(scales$at)
##   ##   }
##   ##   else
##   ##     scales$labels <- abs(scales$labels)
##   ##   if (horizontal) trellis$x.scales <- scales else trellis$y.scales <- scales
##   ##   trellis
##   ## }

##   ## TickLabelsPositive(tmp2)

##   ## TickLabelsPositive(tmpvert, FALSE)



##   likert(Question ~ . | Subtable, data=SFF, layout=c(2,1))
##   likert(Question ~ . | Subtable, data=SFF, layout=c(2,1), xlimEqualLeftRight=TRUE)
##   likert(Question ~ . | Subtable, data=SFF, layout=c(2,1), xTickLabelsPositive=FALSE, xlim=c(-80,110))
##   likert(Question ~ . | Subtable, data=SFF, layout=c(2,1), xlim=c(-80,110))
##   likert(Question ~ . | Subtable, data=SFF, layout=c(2,1), xlimEqualLeftRight=TRUE, xTickLabelsPositive=FALSE)

##   likert(Question ~ . | Subtable, data=SFF, layout=c(2,1))
##   likert(Question ~ . | Subtable, data=SFF, layout=c(2,1), as.percent=TRUE)  ## rightAxisLabels doesn't work with layout=c(2,1)
##   likert(Question ~ . | Subtable, data=SFF, layout=c(1,2), as.percent=TRUE)  ## rightAxisLabels doesn't work with layout=c(1,2)

##   HH:::plot.likert.formula.old(Question ~ . | Subtable, data=SFF, layout=c(2,1), as.percent=TRUE)  ## Works but needs horizontal stretch to be useful
##   HH:::plot.likert.formula.old(Question ~ . | Subtable, data=SFF, layout=c(1,2), as.percent=TRUE)  ## Works

##   likert(Question ~ . | Subtable, data=SFF, layout=c(1,2), as.percent=TRUE)

## }
