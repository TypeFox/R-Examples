library(shiny)
library(OpenMx)
library(rpf)
library(digest)

# drm == dichotomous help TODO
# do something for minItemsPerScore

#options(shiny.trace=TRUE)
#options(shiny.reactlog=TRUE)

verbose <- FALSE

# Use stupid nouns here to encourage folks to change them
# to more sensible ability labels.
sillyFactorName <- c('teacup','puppy','oink','sun','moon')

AppStateTables <- c("recodeTable", "permuteTable", "itemModel", "bayesianPrior")
AppStateSliderInput <- c("numFactors")
AppStateTextInput <- paste0("nameOfFactor", 1:5)

itemToSpec <- function(item) {
#  cat("itemToSpec", item$model, item$outcomes, item$factors, fill=T)
  if (is.null(item$model)) browser()
  outcomes <- item$outcomes
  factors <- item$factors
  # Tc doesn't really matter here
  switch(item$model,
         'drm' = rpf.drm(factors = factors),
         'grm' = rpf.grm(outcomes, factors = factors),
         'nrm' = rpf.nrm(outcomes, factors = factors, T.a=item$Ta,
                         T.c=ifelse(item$Tc == 'partial credit', 'id', item$Tc)))
}

mergeDataAndModel <- function(col, outcomes, factors, item) {
  if (length(outcomes) != 1) browser()
  if (outcomes < 2) return(item)  # will complain later
  if (is.null(item) ||
        (outcomes > 2 && item$model == 'drm') ||
        (outcomes == 2 && item$model == 'nrm')) {
    #    if (verbose) cat("set default model", fill=T)
    spec <- rpf.grm(outcomes, factors = factors)
    np <- rpf.numParam(spec)
    return(list(name=col, model='grm', outcomes=outcomes, factors=factors,
                starting=rpf.rparam(spec), labels=rep(NA, np), free=rep(TRUE, np)))
  }
  if (item$model == 'nrm') {
    if (is.null(item$Ta)) item$Ta <- "trend"
    if (is.null(item$Tc)) item$Tc <- "trend"
  } else {
    item$Ta <- NULL
    item$Tc <- NULL
  }
  item$name <- col
  item$outcomes <- outcomes
  item$factors <- factors
  spec <- itemToSpec(item)
  np <- rpf.numParam(spec)
  if (length(item$starting) != np) {
    newSV <- rpf.rparam(spec)
    preserveMap <- match(names(newSV), item$starting)
    newInd <- which(!is.na(preserveMap))
    preserveMap <- preserveMap[!is.na(preserveMap)]
    oldItem <- item
    item$starting <- newSV
    item$free <- rep(TRUE, np)
    item$labels <- rep(NA, np)
    item$starting[newInd] <- oldItem$starting[preserveMap]
    item$free[newInd] <- oldItem$free[preserveMap]
    item$labels[newInd] <- oldItem$labels[preserveMap]
  }
  item
}

changeItemModel <- function(itemModel, input, outcomes, fi, newModel) {
  oldModel <- isolate(itemModel[[fi]]$model)
  if (is.null(oldModel)) return()
  if (newModel != oldModel) {
    if (verbose) cat("change",fi,"from", oldModel,"to", newModel,fill=T)
    item <- isolate(itemModel[[fi]])
    item$model <- newModel
    numFactors <- isolate(input$numFactors)
    itemModel[[fi]] <- mergeDataAndModel(fi, outcomes, numFactors, item)
  }
}

getFactorNames <- function(input) {
  numFactors <- isolate(input$numFactors)
  if (numFactors == 0) return(c())
  sapply(1:numFactors, function (x) input[[ paste0("nameOfFactor", x) ]])
}

getFocusedItems <- function(input, rawData, permuteTable) {
  dcols <- dataColumnNames(input, rawData, permuteTable)
  range <- sort(match(c(input$focusedItemStart, input$focusedItemEnd), dcols))
  if (length(range) != 2) return(c())
  dcols[seq(range[1], range[2])]
}

buildParameterTable <- function(input, rawData, permuteTable, itemModel, attr) {
  dcols <- getFocusedItems(input, rawData, permuteTable)
  if (length(dcols) == 0) return(NULL)
  massign <- lapply(dcols, function (col) {
    im <- itemModel[[col]]
    v <- im[[attr]]
    names(v) <- names(im$starting)
    v
  })
  names(massign) <- dcols
  tbl <- mxSimplify2Array(massign)
  fnames <- getFactorNames(input)
  if (length(fnames)) {
    rownames(tbl)[1:length(fnames)] <- fnames
  }
  tbl
}

getFocusedItem <- function(input, rawData, permuteTable, itemModel) {
  inames <- getFocusedItems(input, rawData, permuteTable)
  if (length(inames) == 0) return()
  
  im <- itemModel[[ inames[1] ]]
  im
}

getFocusedParameterNames <- function(input, im) {
  pname <- names(im$starting)
  fnames <- getFactorNames(input)
  if (length(fnames)) {
    pname[1:length(fnames)] <- fnames
  }
  pname
}

computeFreeSelected <- function(im, fx) {
  sel <- "Free"
  isFree <- im$free[fx]
  if (!isFree) {
    if (is.finite(im$starting[fx])) {
      sel <- im$starting[fx]
    } else {
      sel <- "inf"
    }
  }
  sel
}

maybeUpdateFree <- function(input, itemModel, im, pname) {
  fx <- match(pname, names(im$starting))
  if (is.na(fx)) return()
  spec <- itemToSpec(im)
  sel <- computeFreeSelected(im, fx)
  if (sel == input$focusedParameterFree) return()
  
  if (input$focusedParameterFree == "Free") {
    if (verbose) cat("free", pname, "of", im$name, fill=TRUE)
    im$free[fx] <- TRUE
    im$starting[fx] <- rpf.rparam(spec)[fx]
  } else {
    if (verbose) cat("fixed", pname, "of", im$name, "to",
                     input$focusedParameterFree, fill=TRUE)
    im$free[fx] <- FALSE
    pi <- rpf.paramInfo(spec, fx)
    if (pi$type == 'bound' && input$focusedParameterFree == 'inf') {
      if (names(im$starting[fx]) == 'g') {
        im$starting[fx] <- qlogis(0)
      } else {
        im$starting[fx] <- qlogis(1)
      }
    } else {
      im$starting[fx] <- as.numeric(input$focusedParameterFree)
    }
  }
  itemModel[[ im$name ]] <- im
}

maybeUpdateLabel <- function(input, itemModel, im, pname) {
  proposal <- isolate(input$focusedParameterLabel)
  newLabel <- ifelse(nchar(proposal), mxMakeNames(proposal), proposal)
  fx <- match(pname, names(im$starting))
  if (is.na(fx)) return()
  sel <- im$labels[fx]
  if (is.na(sel)) sel <- "none"
  if (sel == newLabel) return()
  
  if (newLabel == "none" || !nzchar(newLabel)) {
    im$labels[fx] <- NA
  } else {
    im$labels[fx] <- newLabel
  }
  itemModel[[ im$name ]] <- im
}

toScript <- function(input, rawData, recodeTable, permuteTable, itemModel, bayesianPrior) {
  if (is.null(rawData$name)) {
    loadData <- paste(c("# load some demonstration data",
                        as.character(rawData$loadDemo)), collapse="\n")
  } else {
    loadData <- c("# Adjust the path in the next statement to load your data\n",
                  paste0("data <- read.csv(file='",rawData$name,"'"))
    if (!input$dataHeader) loadData <- c(loadData, ",header=FALSE")
    if (input$dataSep != ",") {
      loadData <- c(loadData, paste0(",sep='", input$dataSep, "'"))
    }
    if (input$dataQuote != '"') {
      loadData <- c(loadData, paste0(',quote="', input$dataQuote, '"'))
    }
    loadData <- c(loadData, ",stringsAsFactors=FALSE,check.names=FALSE)\n",
                  "colnames(data) <- mxMakeNames(colnames(data), unique=TRUE)")
  }
  if (input$freqColumnName != '-') {
    fc <- input$freqColumnName
    loadData <- c(loadData, paste0("\ndata[['", fc, "']] <- as.numeric(data[['", fc, "']])"))
  }
  
  loadData <- c(loadData,
                sapply(rawData$exclude, function(col) paste0("\ndata[['", col, "']] <- NULL  # excluded")))
  
  data <- rawData$val
  dcols <- includedColumnNames(input, rawData, permuteTable)
  mkSpec <- list("spec <- list()")

  mker <- sapply(dcols, function(col) {
    im <- itemModel[[col]]
    str <- c("rpf.", im$model, "(factors=numFactors")
    if (im$model != "drm") {
      str <- c(str, ",outcomes=", im$outcomes)
    }
    if (im$model == "nrm") {
      tc <- im$Tc
      if (tc == 'partial credit') {
        tc <- paste0("lower.tri(diag(", im$outcomes-1, "),TRUE) * -1")
      } else {
        tc <- paste0("'",tc,"'")
      }
      str <- c(str, ",T.a='", im$Ta, "', T.c=", tc)
    }
    paste0(c(str, ")"), collapse="")
  })

  code <- rle(mker)
  val <- code$values
  starts <- cumsum(c(1, code$lengths))[c(rep(TRUE,length(val)),FALSE)]
  lens <- code$lengths
  if (length(starts) == 1) {
    mkSpec <- c(mkSpec, paste0("spec[1:", length(mker),"] <- ",mker[1]))
  } else {
    mkSpec <- c(mkSpec, mapply(function(start, len, v1) {
      paste0("spec[",start,":",(start+len),"] <- ",v1)
    }, starts, lens - 1L, val, SIMPLIFY=FALSE))
  }
  
  mkSpec <- c(mkSpec,
              paste0("names(spec) <- ", paste(deparse(dcols), collapse="\n  ")))
  mkSpec <- paste(mkSpec, collapse="\n")
  
  fnames <- getFactorNames(input)
  
  starting <- mxSimplify2Array(lapply(dcols, function(col) { itemModel[[col]]$starting }))
  free <- mxSimplify2Array(lapply(dcols, function(col) { itemModel[[col]]$free }))
  free[is.na(free)] <- TRUE  # can ignore these
  labels <- mxSimplify2Array(lapply(dcols, function(col) { itemModel[[col]]$labels }))
  labels[!free] <- NA
  
  rnames <- paste0('p', 1:nrow(starting))
  if (length(fnames)) {
    rnames[1:length(fnames)] <- fnames
  }
  
  itemInit <- list()
  for (rx in 1:nrow(free)) {
    rname <- paste0("'", rnames[rx], "'")
    fcode <- rle(!free[rx,])
    mask <- fcode$values
    starts <- cumsum(c(1, fcode$lengths))[c(mask,FALSE)]
    lens <- fcode$lengths[mask]
    if (length(starts) == 0) {
      next
    }
    if (length(starts) == 1 && lens == ncol(free)) {
      itemInit <- c(itemInit, paste0("imat$free[",rname,",] <- FALSE"))
    } else {
      itemInit <- c(itemInit, mapply(function(start, len) {
        paste0("imat$free[",rname,",",start,":",(start+len),"] <- FALSE")
      }, starts, lens - 1L, SIMPLIFY=FALSE))
    }
    srow <- starting[rx,]
    srow[free[rx,]] <- NA
    code <- rle(srow)
    mask <- !is.na(code$value)
    itemInit <- c(itemInit, mapply(function(start, len) {
      paste0("imat$values[",rname,",",start,":",(start+len),"] <- ", srow[start])
    }, cumsum(c(1, code$lengths))[c(mask,FALSE)], code$lengths[mask] - 1L, SIMPLIFY=FALSE))
  }
  
  for (rx in 1:nrow(free)) {
    rname <- paste0("'", rnames[rx], "'")
    code <- rle(labels[rx,])
    mask <- !is.na(code$values)
    val <- code$values[mask]
    starts <- cumsum(c(1, code$lengths))[c(mask,FALSE)]
    lens <- code$lengths[mask]
    if (length(starts) == 0) next
    if (length(starts) == 1 && lens == ncol(free)) {
      itemInit <- c(itemInit, paste0("imat$labels[",rname,",] <- '",val[1],"'"))
    } else {
      itemInit <- c(itemInit, mapply(function(start, len, v1) {
        paste0("imat$labels[",rname,",",start,":",(start+len),"] <- '",v1,"'")
      }, starts, lens - 1L, val, SIMPLIFY=FALSE))
    }
  }
  
  ltbl <- table(labels)
  if (any(ltbl > 1)) {
    itemInit <- c(itemInit,
                  "hasLabel <- !is.na(imat$labels)",
                  paste0("for (lab1 in unique(imat$labels[hasLabel])) {
  imat$values[hasLabel & imat$labels==lab1] <- 
    sample(imat$values[hasLabel & imat$labels==lab1], 1)
}"))
  }
  
  itemInit <- paste(itemInit, collapse="\n")
  
  topModel <- "itemModel"
  fitIFAgroup <- "m1Fit"

  ul <- unique(labels)
  ul <- ul[!is.na(ul)]
  bpList <- match(ul, names(bayesianPrior$map))
  bpList <- sort(bpList[!is.na(bpList)])
  priorInit <- c()
  if (length(bpList)) {
    topModel <- "container"
    fitIFAgroup <- "m1Fit$itemModel"
    bmode <- bayesianPrior$map[bpList]
    bmodeNames <- names(bmode)
    names(bmode) <- NULL
    priorInit <- c(priorInit,
                   "",
                   paste0("priorLabels <- ", paste0(deparse(bmodeNames), collapse="")),
                   "priorMode <- rep(NA, length(priorLabels))")
    code <- rle(as.character(bmode))
    val <- code$values
    starts <- cumsum(c(1, code$lengths))[c(rep(TRUE, length(val)),FALSE)]
    lens <- code$lengths
    if (length(starts) == 1) {
      priorInit <- c(priorInit, paste0("priorMode[1:",length(bmode),"] <- ",bmode[1]))
    } else {
      priorInit <- c(priorInit, mapply(function(start, len, v1) {
        paste0("priorMode[",start,":",(start+len),"] <- ",v1)
      }, starts, lens - 1L, val, SIMPLIFY=FALSE))
    }

    if (input$boundPriorForm == "Beta") {
      priorInit <- c(priorInit,
                     "priorModel <- univariatePrior('beta', priorLabels, priorMode)")
    } else if (input$boundPriorForm == "Logit-normal") {
      priorInit <- c(priorInit,
                     "priorModel <- univariatePrior('logit-norm', priorLabels, priorMode)")
    } else { browser() }
  }
  if (length(priorInit)) {
    priorInit <- c(priorInit,
                   "container <- mxModel(model='container', itemModel, priorModel,
  mxFitFunctionMultigroup(groups=c('itemModel.fitfunction', 'univariatePrior.fitfunction')))")
    priorInit <- paste(priorInit, collapse="\n")
  }
  
  writeFactorNames <- ''
  if (length(fnames)) {
    writeFactorNames <- "rownames(startingValues)[1:numFactors] <- factors"
  }
  
  freqCol <- input$freqColumnName
  maybeCompress <- list()
  if (freqCol == '-') {
    maybeCompress <- "data <- compressDataFrame(data)"
    freqCol <- "freq"
  }
  freqExpectationArgs <- paste0(", weightColumn='", freqCol,"'")
  freqDataArgs <- paste0(", numObs=sum(data[['", freqCol,"']]), sort=FALSE")
  numExtraCol <- 1
  
  getRefModels <- ""
  if (input$fitReferenceModels) {
    getRefModels <- ", refModels=mxRefModels(m1Fit, run = TRUE)"
  }
  
  emArgs <- paste0("'itemModel.expectation', 'scores',
  mxComputeNewtonRaphson(), verbose=", ifelse(input$showFitProgress, "2L", "0L"))
  
  if (input$infoMethod == "*none*") {
    computePlan <- paste0("computePlan <- mxComputeEM(", emArgs, ")")
  } else if (input$infoMethod == "Meat") {
    computePlan <- paste0("emStep <- mxComputeEM(", emArgs,")\n",
                          "computePlan <- mxComputeSequence(list(EM=emStep,
         IM=mxComputeOnce('fitfunction', 'information', 'meat'),
         HQ=mxComputeHessianQuality(),
         SE=mxComputeStandardError()))")
  } else if (input$infoMethod == "Oakes") {
    emArgs <- paste0(emArgs, ",\n  information='oakes1999', infoArgs=list(fitfunction='fitfunction')")
    computePlan <- paste0("emStep <- mxComputeEM(", emArgs,")\n",
                          "computePlan <- mxComputeSequence(list(EM=emStep,
         HQ=mxComputeHessianQuality(),
         SE=mxComputeStandardError()))")
  } else if (input$infoMethod == "Agile SEM") {
    emArgs <- paste0(emArgs, ",\n  information='mr1991',
                     infoArgs=list(fitfunction='fitfunction', semMethod='agile')")
    computePlan <- paste0("emStep <- mxComputeEM(", emArgs,")\n",
                          "computePlan <- mxComputeSequence(list(SE=emStep,
         HQ=mxComputeHessianQuality(),
         SE=mxComputeStandardError()))")
  } else {
    browser()
  }
  
  doPlots <- c()
  if (length(fnames)) {
	  doPlots <- c("",
		       "```{r,fig.height=2}
map1 <- itemResponseMap(m1Grp, factor=1)
ggplot(map1, aes(x=score, y=item, label=outcome)) +
  geom_text(size=4, position=position_jitter(h=.25))
```",
    "",
    "```{r,fig.height=3}
pl <- lapply(names(sfit), function(item) { SitemPlot(sfit, item) })
for (px in 1:length(pl)) {
  print(pl[[px]])
}

basis <- rep(0, length(factors))
basis[1] <- 1
plotInformation(m1Grp, width=5, basis=basis)
```")
  }

  paste(
    "---",
    paste0('title: "',rawData$stem, '"'),
    paste0('date: "', format(Sys.time(), "%d-%b-%Y"), '"'),
    paste0("output: html_document"),
    "---\n",
    "```{r}",
    "options(width=120, scipen=2, digits=2)",
    "suppressPackageStartupMessages(library(OpenMx))",
    "suppressPackageStartupMessages(library(rpf))",
    "suppressPackageStartupMessages(library(ifaTools))",
    "library(xtable)",
    "options(xtable.type='html')",
    "",
    paste0(loadData, collapse=""),
    "",
    paste0("factors <- ", paste0(deparse(fnames)), collapse=""),
    "numFactors <- length(factors)",
    mkSpec,
    "",
    "missingColumns <- which(is.na(match(names(spec), colnames(data))))",
    "if (length(missingColumns)) {",
    "  stop(paste('Columns missing in the data:', omxQuotes(names(spec)[missingColumns])))",
    "}",
    "",
    genRecodeOutcomesCode(input, rawData, recodeTable, permuteTable),
    "",
    "#set.seed(1)   # uncomment to get the same starting values every time",
    "startingValues <- mxSimplify2Array(lapply(spec, rpf.rparam))",
    "rownames(startingValues) <- paste0('p', 1:nrow(startingValues))",
    writeFactorNames,
    "",
    "imat <- mxMatrix(name='item', values=startingValues, free=!is.na(startingValues))",
    itemInit,
    maybeCompress,
    paste0("itemModel <- mxModel(model='itemModel', imat,
           mxData(observed=data, type='raw'", freqDataArgs, "),
           mxExpectationBA81(ItemSpec=spec", freqExpectationArgs, "),
           mxFitFunctionML())"),
    priorInit,
    "",
    computePlan,
    paste0(topModel, " <- mxModel(", topModel, ", computePlan)"),
    "",
    paste0("m1Fit <- mxRun(", topModel, ")"),
    "",
    paste0("m1Grp <- as.IFAgroup(", fitIFAgroup ,", minItemsPerScore=1L)"),
    "```",
    "",
    "An item factor model was fit with `r length(factors)`
factors (`r ifelse(length(factors), factors, '-')`), -2LL=$`r m1Fit$output$fit`$.",
    "The condition number of the information matrix was `r round(m1Fit$output$conditionNumber)`.",
    "",
    "```{r,fig.height=2}
got <- sumScoreEAPTest(m1Grp)
df <- data.frame(score=as.numeric(names(got$observed)),
            expected=got$expected, observed=got$observed)
df <- melt(df, id='score', variable.name='source', value.name='n')
ggplot(df, aes(x=score, y=n, color=source)) + geom_line()
```",
    "",
    "```{r,results='asis'}
ct <- ChenThissen1997(m1Grp)
print(xtable(ct$pval, paste('Log p-value of local dependence between item pairs.')))
```",
    "",
    "```{r,results='asis'}
sfit <- SitemFit(m1Grp)
tbl <- t(sapply(sfit, function(r) c(n=r$n, df=r$df, stat=r$statistic, pval=r$pval)))
print(xtable(tbl, paste0('Sum-score item-wise fit')))
```",
      paste(doPlots, collapse="\n"),
    "",
    paste0("```{r}
summary(m1Fit", getRefModels, ")
```"),
    "",
    "```{r,results='asis'}
citation('OpenMx')
```",
    sep="\n")
}

execRecodeRule <- function(rc, outcomes) {
  if (rc$type == "item") {
    ix <- match(rc$name, names(outcomes))
    if (is.na(ix)) return(outcomes)
    p1 <- setdiff(outcomes[[ix]], rc$from)
    if (rc$to != '') p1 <- union(p1, rc$to)
    outcomes[[ix]] <- sort(p1)
  } else if (rc$type == "outcomeSet") {
    oSetNames <- sapply(outcomes, function(oc) digest(oc, ascii=TRUE))
    for (ix in which(rc$nameHash == oSetNames)) {
      p1 <- setdiff(outcomes[[ix]], rc$from)
      if (rc$to != '') p1 <- union(p1, rc$to)
      outcomes[[ix]] <- sort(p1)
    }
  }
  outcomes
}

recodeOutcomes <- function(input, rawData, recodeTable, permuteTable) {
  dat <- rawData$val
  if (is.null(dat)) return(NULL)
  ch <- dataColumnNames(input, rawData, permuteTable)
  outcomes <- lapply(dat[,ch], function(col) {
    sort(unique(col))
  })
  if (!is.null(recodeTable$val) && nrow(recodeTable$val)) for (rcX in 1:nrow(recodeTable$val)) {
    rc <- recodeTable$val[rcX,]
    outcomes <- execRecodeRule(rc, outcomes)
  }
  outcomes
}

flushFactorTransform <- function(fa1, col) {
  args <- list(col," <- mxFactor(", col, ", ",
               "levels=", paste(deparse(fa1$levels), collapse=""))
  if (any(fa1$levels != fa1$labels)) {
    args <- c(args, ", labels=", paste(deparse(fa1$labels), collapse=""))
  }
  if (length(fa1$exclude)) {
    args <- c(args, ", exclude=", paste(deparse(fa1$exclude), collapse=""))
  }
  if (any(duplicated(fa1$labels))) {
    args <- c(args, ", collapse=TRUE")
  }
  args <- c(args, ")")
  do.call(paste0, args)
}

trackRecodeRule <- function(fa1, ix, rc, col) {
  xf1 <- c()
  from <- as.character(rc$from)
  to <- as.character(rc$to)
  if (from != '' && to != '') {
    fx <- match(from, fa1$levels)
    if (is.na(fx)) {
      xf1 <- c(xf1, flushFactorTransform(fa1, paste0("data[[",col,"]]")))
      newLevels <- sort(fa1$labels)
      fa1 <- list(levels=newLevels, labels=newLevels, exclude=c())
    }
    fx <- match(from, fa1$levels)
    fa1$labels[fx] <- to
  } else if (from == '' && to != '') {
    fa1$levels <- append(fa1$levels, to)
    fa1$labels <- append(fa1$labels, to)
  } else if (from != '' && to == '') {
    fx <- match(from, fa1$levels)
    if (is.na(fx)) {
      xf1 <- c(xf1, flushFactorTransform(fa1, paste0("data[[",col,"]]")))
      fa1 <- list(levels=fa1$labels, labels=fa1$labels, exclude=c())
    }
    fa1$exclude <- append(fa1$exclude, from)
    fx <- match(from, fa1$levels)
    fa1$levels <- fa1$levels[-fx]
    fa1$labels <- fa1$labels[-fx]
  }
  list(fa1=fa1, xf1=xf1)
}

genRecodeOutcomesCode <- function(input, rawData, recodeTable, permuteTable) {
  dat <- rawData$val
  if (is.null(dat)) return(NULL)
  ch <- includedColumnNames(input, rawData, permuteTable)
  outcomes <- lapply(dat[,ch], function(col) {
    sort(unique(col))
  })
  farg <- lapply(outcomes, function(lev) list(levels=lev, labels=lev, exclude=c()))

  xform <- c()
  if (!is.null(recodeTable$val) && nrow(recodeTable$val)) for (rcX in 1:nrow(recodeTable$val)) {
    rc <- recodeTable$val[rcX,]
    
    xf1 <- c()
    if (rc$type == "item") {
      ix <- match(rc$name, names(outcomes))
      if (!is.na(ix)) {
        trr <- trackRecodeRule(farg[[ix]], rc$name, rc, paste0("'",rc$name,"'"))
        farg[[ix]] <- trr$fa1
        xf1 <- c(xf1, trr$xf1)
      }
    } else if (rc$type == "outcomeSet") {
      oSetNames <- sapply(outcomes, function(oc) digest(oc, ascii=TRUE))
      xf2 <- c()
      loop <- names(which(rc$nameHash == oSetNames))
      for (ix in loop) {
        trr <- trackRecodeRule(farg[[ix]], ix, rc, "col")
        farg[[ix]] <- trr$fa1
        xf2 <- trr$xf1   # only need 1 copy
      }
      if (length(xf2)) {
        # tricky and hard to test TODO
        xf1 <- c(xf1, paste0("for (col in ",paste("    ", deparse(loop), collapse="\n"),") {\n",
                             xf2, "}\n"))
      }
    }
    xform <- c(xform, xf1)
    outcomes <- execRecodeRule(rc, outcomes)
  }
  
  farg <- lapply(farg, function(fa) {
    perm <- order(fa$labels)
    fa$levels <- fa$levels[perm]
    fa$labels <- fa$labels[perm]
    fa
  })
  
  oSetNames <- sapply(outcomes, function(oc) digest(oc, ascii=TRUE))
  for (pe in names(permuteTable$val)) {
    if (all(pe != oSetNames)) next
    perm <- permuteTable$val[[pe]]
    for (name in names(oSetNames[pe == oSetNames])) {
      fa <- farg[[name]]
      fa$levels <- fa$levels[perm]
      fa$labels <- fa$labels[perm]
      farg[[name]] <- fa
    }
  }
  
  for (pe in intersect(names(permuteTable$val), names(outcomes))) {
    perm <- permuteTable$val[[pe]]
    fa <- farg[[pe]]
    fa$levels <- fa$levels[perm]
    fa$labels <- fa$labels[perm]
    farg[[pe]] <- fa
  }

  fSetNames <- sapply(farg, function(oc) digest(oc, ascii=TRUE))
  fSetTbl <- table(fSetNames)
  if (length(fSetTbl) == 1) {
    pick <- names(farg)[1]
    xform <- c(xform, flushFactorTransform(farg[[pick]], "data[names(spec)]"))
  } else {
    for (fs in names(fSetTbl)) {
      loop <- names(fSetNames[fs == fSetNames])
      l1 <- loop[1]
      xform <- c(xform, paste0("for (col in ",paste("    ", deparse(loop), collapse="\n"),") {\n  ",
                            flushFactorTransform(farg[[l1]], "data[[col]]"),
                            "\n}"))
      
    }
  }
  
  toRev <- intersect(ch, permuteTable$reversed)
  if (length(toRev)) {
    xform <- c(xform, paste0("# reverse levels of reversed items\n",
                             "for (col in ",paste("    ", deparse(toRev), collapse="\n"),") {\n  ",
                             "data[[col]] <- mxFactor(data[[col]], rev(levels(data[[col]])))",
                             "\n}"))
  }
  
  paste(xform, collapse="\n")
}

getFocusedOutcomeDetailUnordered <- function(input, rawData, recodeTable, permuteTable) {
  outcomes <- recodeOutcomes(input, rawData, recodeTable, permuteTable)  # could accept this as an argument TODO
  if (input$focusedOutcomeItem != '-') {
    fi <- input$focusedOutcomeItem
    iout <- outcomes[[fi]]
    perm <- permuteTable$val[[ digest(iout, ascii=TRUE) ]]
    if (!is.null(perm) && length(perm) == length(iout)) {
      iout <- iout[perm]
    }
    return(iout)
  }
  if (input$focusedOutcomeSet != '-') {
    fi <- input$focusedOutcomeSet
    return(outcomes[[fi]])
  }
  return(NULL)
}

getFocusedOutcomeDetail <- function(input, rawData, recodeTable, permuteTable) {
  outcomes <- getFocusedOutcomeDetailUnordered(input, rawData, recodeTable, permuteTable)
  
  if (input$focusedOutcomeItem != '-') {
    fi <- input$focusedOutcomeItem
  } else if (input$focusedOutcomeSet != '-') {
    fi <- digest(outcomes, ascii=TRUE)
  } else {
    return(outcomes)
  }
  perm <- permuteTable$val[[fi]]
  if (!is.null(perm) && length(perm) == length(outcomes)) {
    outcomes <- outcomes[perm]
  }
  outcomes
}

dataColumnNamesUnsorted <- function(input, rawData) {
  dcol <- colnames(rawData$val)
  if (input$freqColumnName != '-') {
    dcol <- setdiff(dcol, input$freqColumnName)
  }
  dcol
}

dataColumnNames <- function(input, rawData, permuteTable) {
  dcol <- dataColumnNamesUnsorted(input, rawData)
  perm <- permuteTable$items
  if (!is.null(perm) && length(perm) == length(dcol)) dcol <- dcol[perm]
  dcol
}

includedColumnNames <- function(input, rawData, permuteTable) {
  setdiff(dataColumnNames(input, rawData, permuteTable), rawData$exclude)
}

setupItemModels <- function(input, rawData, itemModel, outcomes) {
  numFactors <- isolate(input$numFactors)
  dcol <- isolate(dataColumnNamesUnsorted(input, rawData))
  for (col in dcol) {
    itemModel[[col]] <- mergeDataAndModel(col, length(outcomes[[col]]), numFactors,
                                          isolate(itemModel[[col]]))
  }
}

setupFreqColumnName <- function(session, rawData, freqCol) {
  dcols <- isolate(colnames(rawData$val))
  if (!missing(freqCol) && is.na(match(freqCol, dcols))) browser()
  if (missing(freqCol)) freqCol <- "-"
  updateSelectInput(session, "freqColumnName", selected=freqCol,
                    choices=c("-", dcols))
}

updatePriorMode <- function(input, bayesianPrior, itemModel, im, pname, set) {
  fx <- match(pname, names(im$starting))
  if (is.na(fx)) return()
  oldLabel <- im$labels[fx]
  map <- isolate(bayesianPrior$map)
  
  if (set) {
    if (!(im$model == 'drm' && (pname == 'g' || pname == 'u')))
      return("Can only set prior on dichotomous (drm) bound parameters.")
    prior <- isolate(input$focusedParameterPrior)
    if (is.na(oldLabel) || !nzchar(oldLabel)) {
      oldLabel <- paste0(im$name,"_",pname)
      im$labels[fx] <- oldLabel
      itemModel[[im$name]] <- im
    }
    if (pname == 'g') {
      mode <- paste0("logit(1/",prior,")")
    } else if (pname == 'u') {
      mode <- paste0("logit(",(as.integer(prior)-1),"/",prior,")")
    }
    if (verbose) cat("set bayesian prior", im$name, pname, "to", mode,"using label", oldLabel, fill=TRUE)
    map[[oldLabel]] <- mode
    bayesianPrior$map <- map
  } else {
    if (!is.na(oldLabel) && !is.null(map[[oldLabel]])) {
      if (verbose) cat("clear bayesian prior", im$name, pname, "on label", oldLabel, fill=TRUE)
      map[[oldLabel]] <- NULL
      bayesianPrior$map <- map
    }
  }
  return('')
}

# -----------------------------------------------------------------------------------------
shinyServer(function(input, output, session) {
  feedback <- reactiveValues(newOutcomeAction="", resetRecodeAction="",
                             focusedOutcomeMapAction="", codingFile="",
                             focusedParameterPrior="", parseFile="")
  rawData <- reactiveValues(val=NULL, name=NULL, exclude=NULL)
  recodeTable <- reactiveValues(val=NULL)
  permuteTable <- reactiveValues(val=NULL, items=NULL, reversed=NULL)
  itemModel <- reactiveValues()  # colname to list(model, Ta, Tc, startingValues, free, labels)
  bayesianPrior <- reactiveValues(map=list())  # label -> value map

  observe({
      if (input$exampleDataKCT == 0) return()

      loader <- parse(text=c(
        'utils::data("kct", package="rpf")',
        'data <- kct.people[20:37]',
        'colnames(data) <- mxMakeNames(kct.items[["NAME"]], unique=TRUE)',
        'data <- as.data.frame(lapply(data, as.character), stringsAsFactors=FALSE)',
        'rownames(data) <- kct.people[[19]]'))

      eval(loader)
      rawData$datapath <- NULL
      rawData$name <- NULL
      rawData$loadDemo <- loader
      rawData$stem <- "kct"
      rawData$val <- data
      rawData$exclude <- c("X1x1x4", "X2x2x3", "X3x1x2x4", "X18x4x1x3x4x2x1x4")
      
      recodeTable$val <- 
        structure(list(type = structure(c(1L, 1L, 1L, 1L), .Label = "outcomeSet", class = "factor"), 
                       name = structure(c(1L, 2L, 1L, 1L), .Label = c("X1x1x4", "X18x4x1x3x4x2x1x4"), class = "factor"),
                       nameHash = structure(c(1L,  2L, 3L, 7L), .Label = c("f99c1d57d88bd654baacdcc798a016b6", "8ee3f19dd735a96f0ff59fb514b8efb8", "b695085c9330f407121bb28e859ad508", 
                                                                           "21c943db9d7a66561c37be2e751f7c85", "19eca0e5d319402184b1a1b20eaaa987", 
                                                                           "2572e3e6d229a7405d0ed1839c58f54c", "b867678184d2befa2e1b57c4529dbdc3" ), class = "factor"),
                       action = structure(c(1L, 1L, 2L, 2L), .Label = c("add", "recode"), class = "factor"),
                       from = structure(c(1L, 1L, 2L, 3L), .Label = c("", "0", "1"), class = "factor"), 
                       to = structure(c(1L, 2L, 4L, 3L), .Label = c("0", "1", "correct", "incorrect"), class = "factor")),
                  .Names = c("type", "name", "nameHash", "action", "from", "to"), row.names = c(NA, 4L), class = "data.frame")
      
      permuteTable$val <- 
        structure(list("9b162432e63482dee83b819368c73fb4" = c(2L, 1L)),
                  .Names = "9b162432e63482dee83b819368c73fb4")
      
      isolate({
        outcomes <- isolate(recodeOutcomes(input, rawData, recodeTable, permuteTable))
        setupItemModels(input, rawData, itemModel, outcomes)
        for (name in names(outcomes)) {
          itemModel[[name]]$labels[1] <- 'slope'
        }
      })
      setupFreqColumnName(session, rawData)
  })
  
  observe({
    if (input$exampleDataLSAT6 == 0) return()
    
    loader <- parse(text=c(
      'utils::data("LSAT6", package="rpf")',
      'data <- LSAT6'
    ))
    eval(loader)
    
    rawData$datapath <- NULL
    rawData$name <- NULL
    rawData$loadDemo <- loader
    rawData$stem <- "LSAT6"
    rawData$val <- data

    recodeTable$val <- 
      structure(list(type = structure(c(1L, 1L), .Label = "outcomeSet", class = "factor"), 
                     name = structure(c(1L, 1L), .Label = "Item_1", class = "factor"), 
                     nameHash = structure(1:2, .Label = c("b66901a10dfd2349673b611ab231b535", "b867678184d2befa2e1b57c4529dbdc3"), class = "factor"),
                     action = structure(c(1L,  1L), .Label = "recode", class = "factor"),
                     from = structure(1:2, .Label = c("0", "1"), class = "factor"),
                     to = structure(1:2, .Label = c("incorrect", "correct"), class = "factor")),
                .Names = c("type", "name", "nameHash", "action", "from", "to"), row.names = 1:2, class = "data.frame")
    permuteTable$val <- 
      structure(list("9b162432e63482dee83b819368c73fb4" = c(2L, 1L)), .Names = "9b162432e63482dee83b819368c73fb4")
    
    setupFreqColumnName(session, rawData, "Freq")
  })
  
  observe({
    if (input$exampleDataScience == 0) return()

    loader <- parse(text=c(
      'utils::data("science", package="rpf")',
      'data <- sfpf[20:44]',
      'colnames(data) <- mxMakeNames(colnames(data), unique=TRUE)',
      'data <- as.data.frame(lapply(data, as.character), stringsAsFactors=FALSE)',
      'rownames(data) <- sfpf[[19]]'
    ))
    eval(loader)
    
    rawData$loadDemo <- loader
    rawData$datapath <- NULL
    rawData$name <- NULL
    rawData$stem <- "science"
    rawData$val <- data
    rawData$exclude <- 'GOTOMUSEUM'
    
    recodeTable$val <- 
      structure(list(type = structure(1L, .Label = "outcomeSet", class = "factor"), 
                     name = structure(1L, .Label = "GOTOMUSEUM", class = "factor"), 
                     nameHash = structure(1L, .Label = "7fcb27479582cc445bad8fa1ce4c4c2a", class = "factor"), 
                     action = structure(1L, .Label = "add", class = "factor"), 
                     from = structure(1L, .Label = "", class = "factor"), to = structure(1L, .Label = "dislike", class = "factor")),
                .Names = c("type", "name", "nameHash", "action", "from", "to"), row.names = 1L, class = "data.frame")
    permuteTable$val <- 
      structure(list("37cba13974a597e56737c53035b1a6f0" = c(1L, 3L, 2L)),
                .Names = "37cba13974a597e56737c53035b1a6f0")
    
    setupFreqColumnName(session, rawData)
  })
  
  observe({
    inFile <- input$file1
    
    if (is.null(inFile)) return(NULL)
    
    rawData$datapath <- inFile$datapath
    rawData$name <- inFile$name
  })
  
  output$parseFileFeedback <- renderText(feedback[["parseFile"]])
  
  observe({
    if (is.null(rawData$datapath) || !file.exists(rawData$datapath)) return()
    
    if (0) {
      cat("read.csv with options", input$dataRowNames,
          input$dataSep, input$dataQuote, fill=T)
    }
    args <- list(rawData$datapath, header=input$dataHeader,
                 sep=input$dataSep, quote=input$dataQuote,
                 stringsAsFactors=FALSE, check.names=FALSE)
    
    if (input$dataRowNames) {
      args$row.names=1L
    }
    dat <- try(do.call(read.csv, args), silent = TRUE)
    
    if (inherits(dat, "try-error")) {
      feedback[["parseFile"]] <- paste("Something went wrong trying to load your file:\n", dat)
      return()
    }

    hint <- "Try different quote and separator options."
    if (nrow(dat) == 0) {
      feedback[["parseFile"]] <- paste("Dataset appears to have 0 rows.", hint)
      return()
    }
    if (ncol(dat) == 0) {
      feedback[["parseFile"]] <- paste("Dataset appears to have 0 columns.", hint)
      return()
    }

    colnames(dat) <- mxMakeNames(colnames(dat), unique=TRUE)
    rawData$stem <- sub('\\..{3}$', '', rawData$name, perl=TRUE)
    rawData$val <- dat
    
    setupFreqColumnName(session, rawData)
    feedback[["parseFile"]] <- ""
  })
  
  output$unparsedDataContents <- renderText({
    if (is.null(rawData$datapath)) return()
    paste(readLines(rawData$datapath, n=6L, warn=FALSE), collapse="\n")
  })
  
  output$dataContents <- renderTable({
    dat <- rawData$val
    #    updateTabsetPanel(session, "dataPreviewTabset", "front")  dunno why doesn't work TODO
    head(dat)
  })
  
  output$nameOfDataFile <- renderText({ rawData$name })

  observe({
    ch <- dataColumnNamesUnsorted(input, rawData)
    if (length(ch) == 0) return()
    updateSelectInput(session, "focusedItemStart", choices=ch, selected=ch[[1]])
    updateSelectInput(session, "focusedItemEnd", choices=ch, selected=ch[[ length(ch) ]])
  })
  
  observe({
    perm <- permuteTable$items
    ch <- isolate(dataColumnNamesUnsorted(input, rawData))
    if (!is.null(perm) && length(perm) == length(ch)) ch <- ch[perm]
    updateSelectInput(session, "focusedItemStart", choices=ch,
                      selected=isolate(input$focusedItemStart))
    updateSelectInput(session, "focusedItemEnd", choices=ch,
                      selected=isolate(input$focusedItemEnd))
  })
  
  observe({
    hits <- input$selectAllItemsAction
    if (hits == 0) return()
    
    ch <- isolate(dataColumnNames(input, rawData, permuteTable))
    if (length(ch) == 0) return()
    updateSelectInput(session, "focusedItemStart", selected=ch[[1]])
    updateSelectInput(session, "focusedItemEnd", selected=ch[[ length(ch) ]])
  })
  
  output$numberOfDataRows <- renderText({
    nrow(rawData$val)
  })
  
  output$dataSummary <- renderTable({
    dat <- rawData$val
    if (is.null(dat)) return(NULL)
    
    ch <- dataColumnNames(input, rawData, permuteTable)
    tbl <- sapply(rawData$val[,ch,drop=FALSE], function(col) c(Outcomes=length(unique(col)),
                                                    Missing=sum(is.na(col))))
    tbl <- t(tbl)
    rownames(tbl) <- ch
    tbl
  })
  
  # ------------------------------------------------------------------ Outcomes
  
  observe({
    prevFocus <- isolate(input$focusedOutcomeSet)
    outcomes <- recodeOutcomes(input, rawData, recodeTable, permuteTable)
    oSetHash <- sapply(outcomes, function(oc) digest(oc, ascii=TRUE))
    ch <- names(oSetHash[!duplicated(oSetHash)])
    if (length(ch) == 0) return()
    updateSelectInput(session, "focusedOutcomeSet",
                      choices=c('-', ch),
                      selected=ifelse(!is.na(match(prevFocus, ch)), prevFocus, ch[1]))
  })
  
  observe({
    ch <- dataColumnNames(input, rawData, permuteTable)
    updateSelectInput(session, "focusedOutcomeItem", choices=c("-",ch))
  })
  
  observe({
    if (input$focusedOutcomeItem == '-') return()
    updateSelectInput(session, "focusedOutcomeSet", selected="-")
  })
  
  observe({
    if (input$focusedOutcomeSet == '-') return()
    updateSelectInput(session, "focusedOutcomeItem", selected="-")
  })
  
  output$recodeTable <- renderTable({
    tbl <- recodeTable$val
    if (is.null(tbl) || nrow(tbl) == 0) return()
    tbl[,c("type","name","action","from","to")]
  })

  observe({
    outcomes <- getFocusedOutcomeDetail(input, rawData, recodeTable, permuteTable)
    if (length(outcomes) == 0) outcomes <- "No item selected"
    updateSelectInput(session, "focusedOutcomeMapFrom", choices=outcomes)
    updateSelectInput(session, "focusedOutcomeMapTo", choices=c("<NA>", outcomes, "<Rename>"))
  })
  
  output$addNewOutcomeActionFeedback <- renderText(feedback[["newOutcomeAction"]])

  observe({
    if (input$addNewOutcomeAction == 0) return()
    newOutcome <- isolate(input$newOutcomeName)
    
    if (nchar(newOutcome) == 0) {
      feedback[["newOutcomeAction"]] <- "Enter the name for a new outcome."
      return()
    }

    itemOutcomes <- isolate(getFocusedOutcomeDetailUnordered(input, rawData, recodeTable, permuteTable))

    if (any(newOutcome == itemOutcomes)) {
      feedback[["newOutcomeAction"]] <-
        paste(omxQuotes(newOutcome), "is already an outcome.")
      return()
    }
    
    focusedOutcomeItem <- isolate(input$focusedOutcomeItem)
    focusedOutcomeSet <- isolate(input$focusedOutcomeSet)
    origRecodeTable <- isolate(recodeTable$val)
    
    cat("add outcome", newOutcome, "to item",focusedOutcomeItem,
        "or set",focusedOutcomeSet, fill=T)
    
    if (focusedOutcomeItem != '-') {
      recodeTable$val <-
        rbind(origRecodeTable,
              data.frame(type="item",
                         name=focusedOutcomeItem,
                         nameHash='',
                         action="add",
                         from="",
                         to=newOutcome))
    } else if (focusedOutcomeSet != '-') {
      outcomes <- isolate(recodeOutcomes(input, rawData, recodeTable, permuteTable))
      oSetHash <- sapply(outcomes, function(col) digest(col, ascii=TRUE))
      nameHash <- oSetHash[focusedOutcomeSet == names(oSetHash)][[1]]
      recodeTable$val <-
        rbind(origRecodeTable,
              data.frame(type="outcomeSet",
                         name=focusedOutcomeSet,
                         nameHash=nameHash,
                         action="add",
                         from="",
                         to=newOutcome))
    } else {
      feedback[["newOutcomeAction"]] <- "Select an outcome set or item first."
      return()
    }
    feedback[["newOutcomeAction"]] <- ""
  })
  
  observe({
    if (is.null(recodeTable$val) || nrow(recodeTable$val) == 0) {
      updateNumericInput(session, "focusedRecodeRule", min=1, max=1)
    } else {
      updateNumericInput(session, "focusedRecodeRule", min=1, max=nrow(recodeTable$val))
    }
  })
  
  output$resetRecodeActionFeedback <- renderText(feedback[["resetRecodeAction"]])
  
  observe({
    if (input$resetRecodeAction == 0) return()
    
    origRecodeTable <- isolate(recodeTable$val)
    if (is.null(origRecodeTable) || nrow(origRecodeTable) == 0) {
      feedback[["resetRecodeAction"]] <- "The recode table is empty."
      return()
    }
    
    tbl <- origRecodeTable[-isolate(input$focusedRecodeRule),]
    rownames(tbl) <- NULL
    recodeTable$val <- tbl
    feedback[["resetRecodeAction"]] <- ""
  })
  
  observe({
    if (input$focusedOutcomeMapAction == 0) return()
    
    from <- isolate(input$focusedOutcomeMapFrom)
    to <- isolate(input$focusedOutcomeMapTo)
    if (to == "<Rename>") {
      to <- isolate(input$focusedOutcomeRenameTo)
    }
    if (from == to) {
      feedback[["focusedOutcomeMapAction"]] <-
        paste0("Outcome ",omxQuotes(from)," is already mapped to ",omxQuotes(to),".")
      return()
    }
    if (to == "<NA>") to <- ''
    
    focusedOutcomeItem <- isolate(input$focusedOutcomeItem)
    focusedOutcomeSet <- isolate(input$focusedOutcomeSet)
    origRecodeTable <- isolate(recodeTable$val)

    if (focusedOutcomeItem != '-') {
      recodeTable$val <-
        rbind(origRecodeTable,
              data.frame(type="item",
                         name=focusedOutcomeItem,
                         nameHash='',
                         action="recode",
                         from=from, to=to))
    } else if (focusedOutcomeSet != '-') {
      outcomes <- isolate(recodeOutcomes(input, rawData, recodeTable, permuteTable))
      oSetHash <- sapply(outcomes, function(col) digest(col, ascii=TRUE))
      nameHash <- oSetHash[focusedOutcomeSet == names(oSetHash)][[1]]
      recodeTable$val <-
        rbind(origRecodeTable,
              data.frame(type="outcomeSet",
                         name=focusedOutcomeSet,
                         nameHash=nameHash,
                         action="recode",
                         from=from, to=to))
    } else {
      feedback[["focusedOutcomeMapAction"]] <- "Select an outcome set or item first."
      return()
    }
    feedback[["focusedOutcomeMapAction"]] <- ""
  })

  output$focusedOutcomeMapActionFeedback <- renderText(feedback[["focusedOutcomeMapAction"]])

  output$focusedOutcomeTable <- renderTable(
    data.frame(outcome=getFocusedOutcomeDetail(input, rawData, recodeTable, permuteTable)))
  
  output$reorderOutcomesSorterUI <- renderUI({
    outcomes <- getFocusedOutcomeDetail(input, rawData, recodeTable, permuteTable)
    if (length(outcomes) == 0) {
      return(returnOrder("reorderOutcomesSorter", "Select an item"))
    }
    returnOrder("reorderOutcomesSorter", outcomes)
  })
  
  observe({
    newOrder <- input$reorderOutcomesSorter
    outcomes <- isolate(getFocusedOutcomeDetailUnordered(input, rawData, recodeTable, permuteTable))
    perm <- match(newOrder, outcomes)
    
    focusedOutcomeItem <- isolate(input$focusedOutcomeItem)
    focusedOutcomeSet <- isolate(input$focusedOutcomeSet)
  
    if (focusedOutcomeItem != '-') {
      fi <- focusedOutcomeItem
    } else if (focusedOutcomeSet != '-') {
      fi <- digest(outcomes, ascii=TRUE)
    } else {
      return()
    }
    
    oldPerm <- isolate(permuteTable$val[[fi]])
    if (!is.null(oldPerm) && length(perm) == length(oldPerm) &&
          all(perm == oldPerm)) return()

    if (all(perm == 1:length(perm)) && !is.null(oldPerm)) {
      permuteTable$val[[fi]] <- NULL
    } else if (any(perm != 1:length(perm))) {
      permuteTable$val[[fi]] <- perm
    }
  })
  
  output$permuteTable <- renderTable({
    pt <- permuteTable$val
    if (length(pt) == 0) return()
    t(mxSimplify2Array(pt))
  })
  
  output$reversePicker <- renderUI({
    cols <- dataColumnNames(input, rawData, permuteTable)
    revNames <- isolate(permuteTable$reversed)
    rmask <- !is.na(match(cols, revNames))
    chooserInput("reverseChooser", "Unreversed", "Reversed",
                 cols[!rmask], cols[rmask], size = 16, multiple = TRUE)
  })
  
  observe({
    got <- input$reverseChooser
    if (is.null(got)) return()
    
    oldRev <- isolate(permuteTable$reversed)
    newRev <- setdiff(union(oldRev, got$right), got$left)
    if (length(oldRev) == length(newRev) && all(oldRev == newRev)) return()

    if (verbose) cat("reverse chooser", paste(newRev, collapse=","), fill=TRUE)
    permuteTable$reversed <- newRev
  })
  
  codingToScript <- function() {
    out <- list("savedCodingVersion = 1")
    for (tbl in AppStateTables) {
      dat <- paste(deparse(reactiveValuesToList(get(tbl))), collapse="\n")
      out <- c(out, paste(tbl, "<-", dat))
    }
    inputState <- reactiveValuesToList(input)[c(AppStateSliderInput, AppStateTextInput)]
    out <- c(out, paste("input <-", paste(deparse(inputState), collapse="\n")))
    paste(out, collapse="\n")
  }
  
  output$downloadCoding <- downloadHandler(
    filename = function() {
      paste(rawData$stem, '-config.R', sep='')
    },
    content = function(file) {
      write(isolate(codingToScript()), file=file)
    }
  )
  
  output$codingFileFeedback <- renderText(feedback[["codingFile"]])
  
  output$debugSettingsOutput <- renderText({
    hit <- input$refreshSettingsAction
    isolate(codingToScript())
  })

  observe({
    inFile <- input$codingFile
    if (is.null(inFile)) return()
    bubble <- new.env()
    got <- try(source(inFile$datapath, local=bubble, echo=FALSE), silent=TRUE)
    if (inherits(got, "try-error")) {
      feedback[["codingFile"]] <- paste("Parse error", got, sep="\n")
      return()
    }
    if (!exists(envir=bubble, inherits=FALSE, "savedCodingVersion") ||
          get(envir=bubble, inherits=FALSE, "savedCodingVersion") != 1) {
      feedback[["codingFile"]] <- "This does not seem to be a saved coding file."
      return()
    }
    found <- list()
    e1 <- as.environment(-1L)
    for (tbl in AppStateTables) {
      if (!exists(envir=bubble, inherits=FALSE, tbl)) next
      rt <- mget(tbl, bubble, ifnotfound=list(NULL))[[tbl]]
      dest <- get(tbl)
      for (k in names(rt)) {
        dest[[k]] <- rt[[k]]
      }
      found <- c(found, paste0(tbl,': ', paste0(names(rt), collapse=" ")))
    }
    inputState <- mget("input", bubble, ifnotfound=list(NULL))[["input"]]
    if (!is.null(inputState)) {
      for (elem in AppStateSliderInput) {
        updateSliderInput(session, elem, value=inputState[[elem]])
      }
      for (elem in AppStateTextInput) {
        updateTextInput(session, elem, value=inputState[[elem]])
      }
      found <- c(found, paste0("input",': ',paste0(c(AppStateSliderInput, AppStateTextInput), collapse=" ")))
    }
    feedback[["codingFile"]] <- paste0("Found the following saved settings:\n\n",
                                       paste(found, collapse="\n"))
  })

  # ------------------------------------------------------------------ Item Model & Parameters
  
  observe({
    outcomes <- recodeOutcomes(input, rawData, recodeTable, permuteTable)
    if (is.null(outcomes)) return()
    
    numFactors <- input$numFactors
    if (numFactors) for (fx in 1:numFactors) {
      name <- paste0("nameOfFactor", fx)
      fn <- isolate(input[[name]])
      if (is.null(fn) || fn == "") {
        updateTextInput(session, name, value=sillyFactorName[fx])
      }
    }
    setupItemModels(input, rawData, itemModel, outcomes)
  })
  
  output$reorderItemsSorterUI <- renderUI({
    items <- dataColumnNamesUnsorted(input, rawData)
    if (length(items) == 0) {
      return(returnOrder("reorderItemsSorter", "No data loaded"))
    }
    perm <- isolate(permuteTable$items)
    if (!is.null(perm) && length(items) == length(perm)) items <- items[perm]
    returnOrder("reorderItemsSorter", items)
  })
  
  observe({
    newOrder <- input$reorderItemsSorter
    items <- isolate(dataColumnNamesUnsorted(input, rawData))
    perm <- match(newOrder, items)
    oldPerm <- isolate(permuteTable$items)
    if (!is.null(oldPerm) && length(perm) == length(oldPerm) &&
          all(perm == oldPerm)) return()
    if (all(perm == 1:length(perm)) && !is.null(oldPerm)) {
      permuteTable$items <- NULL
    } else if (any(perm != 1:length(perm))) {
      permuteTable$items <- perm
    }
  })
  
  observe({
    newModel <- input$focusedItemModel
    fi <- isolate(getFocusedItems(input, rawData, permuteTable))
    if (is.null(newModel) || is.null(fi) || newModel == 'as is' || fi == "No data loaded") return()
    if (verbose) cat("change", fi, "item model to", newModel, fill=T)
    outcomes <- isolate(recodeOutcomes(input, rawData, recodeTable, permuteTable))
    for (col in fi) {
      changeItemModel(itemModel, input, length(outcomes[[col]]), col, newModel)
    }
  })
  
  output$focusedItemModelTcFeedback <- renderText(feedback[['focusedItemModelTc']])
  
  observe({
    newTc <- input$focusedItemModelTc
    fi <- getFocusedItems(input, rawData, permuteTable)
    if (is.null(newTc) || is.null(fi) || newTc == 'as is' || fi == "No data loaded") return()
    if (verbose) cat("change", fi, "item model Tc to", newTc, fill=T)
    found <- FALSE
    numFactors <- isolate(input$numFactors)
    for (col in fi) {
      im <- isolate(itemModel[[col]])
      if (is.null(im$model) || im$model != "nrm") next
      found <- TRUE
      if (newTc != im$Tc) {
        if (verbose) cat("change",col,"from", im$Tc,"to", newTc,fill=T)
        im$Tc <- newTc
        itemModel[[col]] <- im
      }
    }
    if (!found) {
      feedback[['focusedItemModelTc']] <- "Intercept basis options are only available for the nominal model."
    } else {
      feedback[['focusedItemModelTc']] <- ''
    }
  })
  
  observe({
    fi <- getFocusedItems(input, rawData, permuteTable)
    if (length(fi) == 0) return()
    Tsel <- 'as is'
    Tchoices <- c('as is', 'trend', 'id', 'partial credit')
    if (length(fi) > 1) {
      sel <- 'as is'
      choices <- c('as is', 'drm', 'grm', 'nrm')
    } else {
      fi <- fi[[1]]
      outcomeMap <- recodeOutcomes(input, rawData, recodeTable, permuteTable)
      outcomes <- length(outcomeMap[[fi]])
      choices <- c('grm')
      if (outcomes == 2) choices <- c('drm', choices)
      if (outcomes > 2) choices <- c(choices, 'nrm')
      im <- isolate(itemModel[[fi]])
      sel <- im$model
      if (length(sel) && sel == 'nrm') {
        Tsel <- im$Tc
      } else {
        Tchoices <- 'as is'
      }
    }
    updateSelectInput(session, "focusedItemModel", choices=choices, selected=sel)
    updateSelectInput(session, "focusedItemModelTc", choices=Tchoices, selected=Tsel)
  })
  
  output$itemStartingValuesTable <- renderTable({
    if (length(names(itemModel)) == 0) return(NULL)
    sv <- buildParameterTable(input, rawData, permuteTable, itemModel, "starting")
    
    # adjust for equality constraints
    labels <- buildParameterTable(input, rawData, permuteTable, itemModel, "labels")
    ltbl <- table(labels)
    if (any(ltbl > 1)) {
      hasLabel <- !is.na(labels)
      for (lab in names(ltbl[ltbl > 1])) {
        sv[hasLabel & labels== lab] <-
          sample(sv[hasLabel & labels==lab], 1)
      }
    }

    sv
  })
  
  output$itemFreeTable <- renderTable({
    if (length(names(itemModel)) == 0) return(NULL)
    buildParameterTable(input, rawData, permuteTable, itemModel, "free")
  })
  
  output$itemLabelTable <- renderTable({
    if (length(names(itemModel)) == 0) return(NULL)
    buildParameterTable(input, rawData, permuteTable, itemModel, "labels")
  })
  
  output$itemPriorTable <- renderTable({
    if (length(names(itemModel)) == 0) return(NULL)
    tbl <- buildParameterTable(input, rawData, permuteTable, itemModel, "labels")
    map <- bayesianPrior$map
    mapped <- match(tbl[,], names(map))
    mapped[!is.na(mapped)] <- unlist(map[ mapped[!is.na(mapped)] ])
    tbl[,] <- mapped
    for (rx in rownames(tbl)) {
      if (rx == 'g' || rx == 'u') {
        # OK
      } else {
        tbl[rx,] <- sapply(tbl[rx,], function(pr)  ifelse(is.na(pr), NA, "??"))
      }
    }
    tbl
  })
  
  observe({
    im <- getFocusedItem(input, rawData, permuteTable, itemModel)
    if (is.null(im)) return()
    pname <- getFocusedParameterNames(input, im)
    prevSelect <- match(isolate(input$focusedItemParameter), pname)
    if (is.na(prevSelect)) prevSelect <- 1
    updateSelectInput(session, "focusedItemParameter",
                      choices=pname, selected=pname[prevSelect])
  })
  
  observe({
    im <- getFocusedItem(input, rawData, permuteTable, itemModel)
    if (is.null(im)) return()
    fx <- match(input$focusedItemParameter,
                isolate(getFocusedParameterNames(input, im)))
    if (is.na(fx)) return()
    
    spec <- itemToSpec(im)
    if (fx > rpf.numParam(spec)) return()

    choices <- c("Free", 0, 1)
    pi <- rpf.paramInfo(spec, fx)
    if (pi$type == "bound") choices <- c(choices, "inf")
    sel <- computeFreeSelected(im, fx)
    
    allFocused <- length(getFocusedItems(input, rawData, permuteTable)) > 1
    if (allFocused) {
      choices <- c("as is", choices)
      sel <- "as is"
    }
    updateSelectInput(session, "focusedParameterFree", choices=choices,
                      selected=sel)

    label <- im$labels[fx]
    if (is.na(label)) label <- "none"
    if (allFocused) {
      label <- "as is"
    }
    updateTextInput(session, "focusedParameterLabel", value=label)
  })

  observe({
    if (input$focusedParameterFree == "as is" ||
          input$focusedParameterFree == "No parameter selected") return()
    
    im <- isolate(getFocusedItem(input, rawData, permuteTable, itemModel))
    if (is.null(im)) return()
    origFx <- isolate(match(input$focusedItemParameter,
                            getFocusedParameterNames(input, im)))
    if (is.na(origFx)) return()
    
    pname <- names(im$starting)[origFx]
    
    dcols <- isolate(getFocusedItems(input, rawData, permuteTable))
    for (col in dcols) {
      maybeUpdateFree(input, itemModel, isolate(itemModel[[col]]), pname)
    }
  })
  
  observe({
    hit <- input$changeLabelAction
    
    if (hit == 0 || isolate(input$focusedParameterLabel) == "as is") return()

    im <- isolate(getFocusedItem(input, rawData, permuteTable, itemModel))
    if (is.null(im)) return()
    origFx <- isolate(match(input$focusedItemParameter,
                            getFocusedParameterNames(input, im)))
    if (is.na(origFx)) return()
    pname <- names(im$starting)[origFx]
    
    dcols <- isolate(getFocusedItems(input, rawData, permuteTable))
    for (col in dcols) {
      maybeUpdateLabel(input, itemModel, isolate(itemModel[[col]]), pname)
    }
  })
  
  output$focusedParameterPriorFeedback <- renderText(feedback[["focusedParameterPrior"]])
  
  observe({
    hit <- input$focusedParameterPriorSetAction
    
    if (hit == 0) return()
    
    im <- isolate(getFocusedItem(input, rawData, permuteTable, itemModel))
    if (is.null(im)) return()
    origFx <- isolate(match(input$focusedItemParameter,
                            getFocusedParameterNames(input, im)))
    if (is.na(origFx)) return()
    pname <- names(im$starting)[origFx]

    dcols <- isolate(getFocusedItems(input, rawData, permuteTable))
    for (col in dcols) {
      feedback[["focusedParameterPrior"]] <- 
        updatePriorMode(input, bayesianPrior, itemModel, isolate(itemModel[[col]]), pname, set=TRUE)
    }
  })
  
  observe({
    hit <- input$focusedParameterPriorClearAction
    
    if (hit == 0) return()
    
    im <- isolate(getFocusedItem(input, rawData, permuteTable, itemModel))
    if (is.null(im)) return()
    origFx <- isolate(match(input$focusedItemParameter,
                            getFocusedParameterNames(input, im)))
    if (is.na(origFx)) return()
    pname <- names(im$starting)[origFx]
    
    dcols <- isolate(getFocusedItems(input, rawData, permuteTable))
    for (col in dcols) {
      updatePriorMode(input, bayesianPrior, itemModel, itemModel[[col]], pname, set=FALSE)
    }
    feedback[["focusedParameterPrior"]] <- ''
  })
  
  output$excludePicker <- renderUI({
    cols <- dataColumnNames(input, rawData, permuteTable)
    exNames <- isolate(rawData$exclude)
    emask <- !is.na(match(cols, exNames))
    chooserInput("excludeChooser", "Included", "Excluded",
                 cols[!emask], cols[emask], size = 16, multiple = TRUE)
  })

  observe({
    got <- input$excludeChooser
    if (is.null(got)) return()
    
    oldRev <- isolate(rawData$exclude)
    newRev <- setdiff(union(oldRev, got$right), got$left)
    if (length(oldRev) == length(newRev) && all(oldRev == newRev)) return()
    
    if (verbose) cat("exclude chooser", paste(newRev, collapse=","), fill=TRUE)
    rawData$exclude <- newRev
  })
  
  output$itemModelAssignment <- renderTable({
    if (length(names(itemModel)) == 0) return(NULL)
    
    outcomes <- recodeOutcomes(input, rawData, recodeTable, permuteTable)
    massign <- mapply(function (col, cname) {
      im <- itemModel[[cname]]
      c(Outcomes=length(col),
        Model=im$model, T.a=im$Ta, T.c=im$Tc,
        Reversed=cname %in% permuteTable$reversed,
        Excluded=cname %in% rawData$exclude)
    }, outcomes, names(outcomes), SIMPLIFY=FALSE)
    tbl <- t(mxSimplify2Array(massign))
    for (bcol in c('Reversed','Excluded')) {
      tbl[tbl[,bcol] == "FALSE", bcol] <- NA
    }
    tbl
  })
  
  # ------------------------------------------------------------------ Preview & Download
  
  output$debugScriptOutput <- renderText({
    input$debugScriptAction
    validate(need(nrow(rawData$val), "No data is loaded"))
    str <- try(isolate(toScript(input, rawData, recodeTable, permuteTable, itemModel, bayesianPrior)), silent=TRUE)
    validate(need(!inherits(str, "try-error"),
                  paste("An error occurred. Please report this to the development team.",
                        str, sep="\n")))
    str
  })
  
  output$downloadScript <- downloadHandler(
    filename = function() {
      validate(need(!is.null(rawData$stem), "No data is loaded"))  # TODO, fix
      paste(rawData$stem, '.Rmd', sep='')
    },
    content = function(file) {
      write(isolate(toScript(input, rawData, recodeTable, permuteTable, itemModel, bayesianPrior)), file=file)
    }
  )
})
