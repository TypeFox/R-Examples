## this file is a redesign with extensions of function in ae.dotplot
## in file ae.dotplot.R

AEdotplot <- function(xr, ...)
  UseMethod("AEdotplot")

AEdotplot.data.frame <- function(xr, ...,
                                 conditionVariable=NULL,
                                 conditionName=deparse(substitute(xr)),
                                 useCondition=!is.null(conditionVariable),
                                 sub=list(conditionName, cex=.7)) {
  if (!all(c("RAND", "PREF", "SAE", "SN") %in% names(xr)))
    stop('In the AEdotplot.data.frame method,\n variables named c("RAND", "PREF", "SAE", "SN") must be used.', call.=FALSE)
  if (length(conditionName) != 1 || !is.character(conditionName))
    stop("conditionName must be a character vector of length one", call.=FALSE)
  if ("subset" %in% names(list(...)))
    stop("'subset=' doesn't work here.  Use explicit subscript.", call.=FALSE)
  if (is.null(conditionVariable) || !useCondition) {
    KEEP <- levels(xr$PREF) %in% unique(xr$PREF)
    xr$PREF <- ordered(xr$PREF, levels=levels(xr$PREF)[KEEP])
    xaer <- AElogrelrisk(xr, ...)
    result <- AEdotplot(xaer, ..., conditionName=conditionName, sub=sub)
    xaer <- list(xaer)
    class(xaer) <- c("AEtable", class(xaer))
    names(xaer) <- conditionName
  }
  else {
    xaer <- lapply(split(xr, conditionVariable, drop=TRUE),
                   function(xri) {
                     KEEP <- levels(xri$PREF) %in% unique(xri$PREF)
                     xri$PREF <- ordered(xri$PREF, levels=levels(xri$PREF)[KEEP])
                     AElogrelrisk(xri, ...)
                   })
    class(xaer) <- c("AEtable", class(xaer))
    tmp <- lapply(xaer, AEdotplot, ..., sub=sub)
    result <- do.call("c", tmp)
    for (i in seq(length(result)))
      result[[i]]$par.strip.text <- tmp[[1]][[1]]$par.strip.text
  }
  attr(result, "AEtable") <- xaer
  result
}

AEdotplot.AElogrelrisk <-
  function(xr,
           A.name=paste(levels(xr$RAND)[1], " (n=", xr$SN[1], ")", sep=""),
           B.name=paste(levels(xr$RAND)[2], " (n=", xr$SN[2], ")", sep=""),
           col.AB=c("red","blue"), pch.AB=c(16,17),
           main=if (sortbyRelativeRisk)
              list("Most Frequent On-Therapy Adverse Events Sorted by Relative Risk", cex=1)
           else
              list("Most Frequent On-Therapy Adverse Events", cex=1),
           cex.AB.points=NULL, cex.AB.y.scale=.6, cex.x.scale=.6,
           panel.widths=c(.55, .22, .23),
           key.y=-.2, CI.percent=95, conditionName=deparse(substitute(xr)),
           sortbyRelativeRisk=TRUE,
           ...,
           sub=list(conditionName, cex=.7), par.strip.text=list(cex=.7)) {

    ae.key <- list(y = key.y, x=.15,
                   points = list(col=col.AB, pch=pch.AB),
                   text = list(c(A.name, B.name), col=col.AB, cex=.8),
                   columns = 2,
                   between=.5,
                   space="bottom") ## r: bottom, s: top

    xr <- cbind(xr,
                condition=conditionName,
                left.x="Percent",
                right.x=paste("Relative Risk with ", CI.percent, "% CI", sep=""),
                text.x="nAE, pct, Relative Risk",
                stringsAsFactors=FALSE)
    if (any(table(xr$PREF) != 2))
       stop("At least one Adverse Effect has other than two Treatment levels.",
            call.=FALSE)
    ## construct left panel
    left.plot <- dotplot(PREF ~ PCT | left.x + condition,
                         data=xr, groups = xr$RAND,
                         col=col.AB, pch=pch.AB,
                         panel = panel.ae.leftplot,
                         xlab=NULL,
                         scales=list(
                           x=list(
                             cex=cex.x.scale,
                             limits=range(c(-.8, xr$PCT+1))),
                           y=list(cex=cex.AB.y.scale)),
                         par.strip.text=par.strip.text
                         )

    ## construct right panel
    if (is.null(xr$logrelrisk))
      stop("Variable 'logrelrisk' missing.\nPlease use the logrelisk() function before using ae.dotplot().")
    A.rows <- seq(1, nrow(xr), 2)
    B.rows <- A.rows + 1
    right.plot <- dotplot(PREF ~ logrelrisk | right.x + condition,
                          data=xr[A.rows,], pch=16,
                          lower=xr[A.rows, "logrelriskCI.lower"],
                          upper=xr[A.rows, "logrelriskCI.upper"],
                          panel=panel.ae.rightplot,
                          xlab=NULL,
                          scales=list(
                            x=list(
                              cex=cex.x.scale,
                              at=log(c(.125, .25, .5, 1, 2, 4, 8, 16, 32)),
                              labels=c(".125", "", ".5", 1, 2, 4, 8, 16, 32),
                              limits=range(log(c(.0625, 64)),
                                c(xr$logrelriskCI.lower, xr$logrelriskCI.upper),
                                finite=TRUE)),
                            y=list(cex=cex.AB.y.scale)
                            ),
                          par.strip.text=par.strip.text,
                          par.settings=list(layout.widths=list(
                                              left.padding=0,
                                              ylab.axis.padding=0))
                          )

    ## construct text panel
    xrA <- xr[A.rows, c("PREF", "SAE", "PCT"           )]
    names(xrA)[2:3] <- c("A:nn", "A:pctn")
    xrB <- xr[B.rows, c(        "SAE", "PCT", "relrisk")]
    names(xrB)[1:3] <- c("B:nn", "B:pctn", "relriskn")
    xrwide <- cbind(xrA, xrB)
    xrwidechar <- data.frame(xrwide,
                             "A:n"    =             xrwide[,"A:nn"],
                             "A:pct"  =format(round(xrwide[,"A:pctn"], 2)),
                             "B:n"    =             xrwide[,"B:nn"],
                             "B:pct"  =format(round(xrwide[,"B:pctn"], 2)),
                             "relrisk"=format(round(xrwide[,"relriskn"], 2)),
                             stringsAsFactors=FALSE,
                             check.names=FALSE
                             )
    `nAE, pct, Relative Risk`  <- factor(names(xrwidechar)[7:11], levels=names(xrwidechar)[7:11])
    xrcharlong <- data.frame(PREF=rep(xrwide$PREF, 5),
                             `nAE, pct, Relative Risk` =rep(`nAE, pct, Relative Risk` , each=length(levels(xrwide$PREF))),
                             text.x=xr$text.x[1],
                             value=unlist(xrwidechar[7:11]),
                             condition=xr$condition[1],
                             stringsAsFactors=FALSE,
                             check.names=FALSE
                             )
    text.plot <- with(xrcharlong,
                      xyplot(PREF ~ `nAE, pct, Relative Risk` | text.x + condition,
                             labels=value,
                             panel=panel.text,
                             xlab=NULL,
                             ylab=NULL,
                             cex=.5, adj=1,
                             scales=list(
                               x=list(
                                 cex=cex.x.scale,
                                 at=(1:5)-.1, labels=levels(`nAE, pct, Relative Risk` )),
                               y=list(cex=cex.AB.y.scale)
                               ),
                             par.strip.text=par.strip.text,
                             par.settings=list(layout.widths=list(
                                                 left.padding=0,
                                                 ylab.axis.padding=0))
                             ))
    result <- list(left.plot= useOuterStrips(left.plot),
                   right.plot=useOuterStrips(right.plot),
                   text.plot= useOuterStrips(text.plot))
    attr.list <- list(ae.key=ae.key,
                      main=main,
                      sub=sub,
                      n.events=length(levels(xr$PREF)),
                      panel.widths=panel.widths)
    attributes(result)[names(attr.list)] <- attr.list
    attr(result, "AEtable") <- list(xr)
    names(attr(result, "AEtable")) <- conditionName
    class(result) <- "AEdotplot"
    result
  }




print.AEdotplot <- function(x, ...,
                            main=attr(x, "main"),
                            sub=attr(x,"sub"),
                            ae.key=attr(x, "ae.key"),
                            panel.widths=attr(x,"panel.widths"),
                            AEtable=TRUE) {

  x.in <- x
  if (!AEtable && missing(panel.widths)) panel.widths <- c(.7, .3, 0)
  title.centers <- ifelse (panel.widths[3] > 0, .67, .90)

  title.adjust <- function(x, just) {
    if (is.null(x) ||
        is.na(x) ||
        (is.logical(x) && (x == FALSE)) ||
        (!is.list(x) && nchar(x) == 0) ||
        (is.list(x) && nchar(x[[1]]) == 0)) {
      x <- NULL
      x.blank <- x
    }
    else {
      if (!is.list(x)) x <- list(x, cex=1)
      if (!any(names(x) == "label")) {
        empty <- which(names(x) == "")
        if (length(empty) > 0) names(x)[empty[1]] <- "label"
      }
      x.blank <- x
      x.blank[["label"]] <- " "
      if (length(grep("\n", x[["label"]])) > 0)
        x.blank[["label"]] <-
          paste(" ", rep("\n", length(gregexpr("\n", x[["label"]])[[1]])),
                collapse="", sep="")
      x$just <- just
    }
    list(x, x.blank)
  }

  main.both <- title.adjust(main, just=title.centers)
  x$left.plot$main    <- main.both[[2]]
  x$right.plot$main   <- main.both[[1]]
  x$text.plot$main    <- main.both[[2]]

  sub.both <- title.adjust(sub, just=title.centers)
  x$left.plot$sub     <- sub.both[[2]]
  x$right.plot$sub    <- sub.both[[1]]
  x$text.plot$sub     <- sub.both[[2]]

  key.blank <- ae.key
  key.blank$points$col <- 0
  key.blank$text[[1]] <- c(" ", " ")

  x$left.plot$legend  <- list(bottom=list(
                                fun="draw.key",
                                args=list(key=ae.key)))
  x$right.plot$legend <- list(bottom=list(
                                fun="draw.key",
                                args=list(key=key.blank)))
  x$text.plot$legend  <- list(bottom=list(
                                fun="draw.key",
                                args=list(key=key.blank)))

  pw <- cumsum(c(0, panel.widths))
  pos1 <- c(pw[1], 0, pw[2], 1)
  pos2 <- c(pw[2], 0, pw[3], 1)
  pos3 <- c(pw[3], 0, pw[4], 1)
  if ((pos1[3] - pos1[1]) > 0)
    print(       x$left.plot                                 ,
          position=pos1, more=TRUE)

  if ((pos3[3] - pos3[1]) > 0)
    print(update(x$text.plot, scales=list(y=list(draw=FALSE)), strip.left=FALSE),
          position=pos3, more=TRUE)

  ## print right.plot (in the middle) last because it holds main, sub, and key
  if ((pos2[3] - pos2[1]) > 0)
    print(update(x$right.plot, scales=list(y=list(draw=FALSE)), strip.left=FALSE),
          position=pos2, more=TRUE)

  ##  lattice.setStatus() ## needed because all three panels are set to TRUE
  ## lattice:::lattice.setStatus(print.more = FALSE)
  lattice.lattice.setStatus(print.more = FALSE)

  invisible(x.in)
}


printOld.AEdotplot <- function(x, ...,
                           main.x.center=.6,
                           panel.widths=attr(x,"panel.widths")) {
### this is the minimalist version.  it prints three coordinated "trellis" objects.
### The main and sub and key are independently printed.
  pw <- cumsum(c(0, panel.widths))
  pos1 <- c(pw[1], 0.02, pw[2], 1)
  pos2 <- c(pw[2], 0.02, pw[3], 1)
  pos3 <- c(pw[3], 0.02, pw[4], 1)
  if ((pos1[3] - pos1[1]) > 0)
    print(       x[[1]]                                  ,
          position=pos1, more=TRUE )
  if ((pos2[3] - pos2[1]) > 0)
    print(update(x[[2]], scales=list(y=list(draw=FALSE)), strip.left=FALSE),
          position=pos2, more=TRUE )
  if ((pos3[3] - pos3[1]) > 0)
    print(update(x[[3]], scales=list(y=list(draw=FALSE)), strip.left=FALSE),
          position=pos3, more=TRUE)

  draw.key(list(text=attr(x, "main"), font=par("font.main")), draw=TRUE,
           vp=viewport(x=main.x.center - .1, y=.97))
  draw.key(list(text=attr(x, "sub"), font=par("font.main")), draw=TRUE,
           vp=viewport(x=main.x.center, y=.02))
  draw.key(attr(x, "ae.key"), draw=TRUE,
           vp=viewport(x=main.x.center, y=.045))

  ##  lattice.setStatus() ## needed because all three panels are set to TRUE
  ## lattice:::lattice.setStatus(print.more = FALSE)
  lattice.lattice.setStatus(print.more = FALSE)

  invisible(x)
}


AEdotplot.AEtable <- function(xr, ..., useCondition=TRUE,
                              sub="sub for AEsecond") {
  class(xr) <- class(xr)[-1]
  result <- list()
  for (i in names(xr)) {
    result[[i]] <- AEdotplot(xr[[i]], ..., conditionName=i, sub=sub)
    attr(result[[i]], "AEtable") <- xr[[i]]
  }
  do.call("c", result)
}

c.AEdotplot <- function(...,
                        panel.widths=attr(aedp[[1]], "panel.widths"),
                        par.strip.text=list(cex=.7)) {

  aedp <- list(...) ## (named) list of AEdotplot objects

  n.events <- unlist(sapply(aedp, attr, "n.events"))
  condlevels.2 <- unlist(sapply(aedp, function(x) x[[1]]$condlevels[[2]] ))
  if (!is.null(names(condlevels.2))) {
    condlevels.2 <-
      ifelse(nchar(names(condlevels.2))>0,
             names(condlevels.2),
             condlevels.2)
  }
  aedps <- vector("list", 3)
  names(aedps) <- names(aedp[[1]])
  for (i in names(aedps)) {
    aei <- do.call("c",  lapply(aedp, `[[`, i))
    aei <- update(aei,  layout=c(1, length(aedp)))
    aei$condlevels <- list(aedp[[1]][[i]]$condlevels[[1]], condlevels.2)
    aei$index.cond <- list(1, 1:length(condlevels.2))
    aei$perm.cond <- 1:2
    aei <- update(aei, xlab=NULL, between=list(y=1))
    aei <- useOuterStrips(aei)
    aei <- update(aei, par.strip.text=par.strip.text)
    aei <- resizePanels(aei, h=n.events)
    aedps[[i]] <- aei
  }
  aedps$left.plot <- update(aedps$left.plot,
                            scales=list(x=list(relation="same")))   ## required
  aedps$right.plot <- update(aedps$right.plot,
                            scales=list(x=list(relation="same")),
                             par.settings=list(layout.widths=list(
                                                 left.padding=0,
                                                 ylab.axis.padding=0)))
  aedps$right.plot$x.limits <- aedp[[1]][[2]]$x.limits
    ##rep(list(aedp[[1]][[2]]$x.limits),
    ##    length(aedps$right.plot$packet.sizes))
  aedps$right.plot$x.used.at <- NULL
  aedps$right.plot$x.num.limit <- NULL
  aedps$text.plot <- update(aedps$text.plot,
                            scales=list(x=list(relation="same")),   ## not original, not sure
                             par.settings=list(layout.widths=list(
                                                 left.padding=0,
                                                 ylab.axis.padding=0)))
  attributes(aedps) <- attributes(aedp[[1]])
  names(n.events) <- condlevels.2
  attr(aedps, "n.events") <- n.events
  attr(aedps, "panel.widths") <- panel.widths
  xaer <- unlist(lapply(aedp, function(x) attr(x, "AEtable")), recursive=FALSE)
  class(xaer) <- c("AEtable", class(xaer))
  attr(aedps, "AEtable") <- xaer
  aedps
}



## Calculate the percent, the relative risk, the log relative risk,
## and the confidence intervals.  if (sortbyRelativeRisk) {make PREF
## an ordered factor, sorted by the relative risk} else {use the order
## implied by levels(PREF) as the order.  This will normally be
## alphabetical unless the user has taken control}.  The user may
## instead specify the variable for sorting.  The names are the names
## used inside this function.
AElogrelrisk <- function(ae,
                         A.name=levels(ae$RAND)[1],
                         B.name=levels(ae$RAND)[2],
                         crit.value=1.96,
                         sortbyRelativeRisk=TRUE, ...,
                         sortbyVar=c("PREF", ## Event name
                           "PCT",            ## Percent
                           "SN",             ## Number of Patients
                           "SAE",            ## Number of Observed Events
                           "relrisk",        ## Relative Risk (RR)
                           "ase.logrelrisk", ## Asymptotic Standard Error(log(RR))
                           "relriskCI.lower","relriskCI.upper"), ## Confidence Interval Bounds
                         sortbyVarBegin=1) {

  if (any(ae$SN <= 0))
    stop("At least one AE has 0 patients at risk.", call.=FALSE)
  sortbyVar <- match.arg(sortbyVar)
  ae$PCT <- 100 * ae$SAE / ae$SN ## percent of patients

  ## Calculation of relative risk assumes both A and B PREF are in the
  ## same order.  Assignment of relative risk to the data.frame assumes
  ## that A and B for each PREF are in consecutive positions.
  ## This ordering enforces that assumption.
  ae <- ae[with(ae, order(PREF, RAND)),]
  AA <- ae$PREF[ae$RAND==A.name]
  BB <- ae$PREF[ae$RAND==B.name]
  if (length(AA) != length(BB) || !all(AA == BB))
    stop("Data problem: at least one AE is missing an A treatment or a B treatment.", call.=FALSE)

  tmp <-  ## sample relative risk
    ae$PCT[ae$RAND==B.name] /
      ae$PCT[ae$RAND==A.name]
  names(tmp) <- ae$PREF[seq(1,nrow(ae),2)]
  if (sortbyRelativeRisk) {
    ## sort by relrisk
    ae$PREF <- ordered(ae$PREF, levels=unique(names(sort(tmp, na.last=FALSE)))) ## allow multiple pairs
  }
  ae$relrisk <- as.vector(rbind(tmp,tmp))
  ae$logrelrisk <- log(ae$relrisk)

  ase.logrelrisk <-
    ## sample asymptotic standard error of relative risk
    ##  (Agresti, Equation 3.18)
    sqrt((1-ae$PCT[ae$RAND==B.name]/100) /
         (ae$PCT[ae$RAND==B.name]/100 *
          ae$SN[ae$RAND==B.name])
         + (1-ae$PCT[ae$RAND==A.name]/100) /
         (ae$PCT[ae$RAND==A.name]/100   *
          ae$SN[ae$RAND==A.name]  ))
  ae$ase.logrelrisk <- as.vector(rbind(ase.logrelrisk, ase.logrelrisk))

  ae$logrelriskCI.lower <- ae$logrelrisk - crit.value*ae$ase.logrelrisk
  ae$logrelriskCI.upper <- ae$logrelrisk + crit.value*ae$ase.logrelrisk
  ae$relriskCI.lower <- exp(ae$logrelriskCI.lower)
  ae$relriskCI.upper <- exp(ae$logrelriskCI.upper)

  if (!is.null(sortbyVar)) {
    sortrows <- seq(sortbyVarBegin, nrow(ae), 2)
    sortby <- ae[sortrows, c("PREF", sortbyVar)]
    sortby[order(sortby[[sortbyVar]]), "PREF"]
    ae$PREF <- ordered(ae$PREF, levels=sortby[order(sortby[[sortbyVar]]), "PREF"])
  }

  class(ae) <- c("AElogrelrisk", class(ae))
  ae
}

AEmatchSortorder <-
  function(AEstandard,
           AEsecond,
           AEsecond.AEtable=attr(AEsecond, "AEtable"),
           levels.order=lapply(attr(AEstandard,"AEtable"),
             function(AEsubtable)
             levels(AEsubtable$PREF)),
                    main.second=list(paste("Most Frequent On-Therapy Adverse Events",
                                          "Sorted to Match First Table"),
                                     cex=1)) {
    for (subgroup in 1:length(AEsecond.AEtable))
      AEsecond.AEtable[[subgroup]]$PREF <-
        ordered(AEsecond.AEtable[[subgroup]]$PREF,
                levels=levels.order[[subgroup]])
    newsub <- if (missing(AEsecond))
      "sub for AEsecond"
    else
      attr(AEsecond, "sub")
    AEsecond.matched <- AEdotplot(AEsecond.AEtable, useCondition=TRUE,
                                  main=main.second,
                                  sub=newsub)
    ## order of AE is now the same in AEstandard and AEsecond.matched.
    ## xlim may be different.
    ## xlim may now be manually adjusted with
    xlim.12 <- range(AEstandard$left.plot$x.limits,
                     AEsecond.matched$left.plot$x.limits)
    AEstandard$left.plot <- update(AEstandard$left.plot, xlim=xlim.12)
    AEsecond.matched$left.plot <- update(AEsecond.matched$left.plot, xlim=xlim.12)
    list(AEstandard=AEstandard, AEsecond.matched=AEsecond.matched)
  }

update.AEdotplot <- function(object, ... ) {
  AEnew <- lapply(object, update, ...)
  attributes(AEnew) <- attributes(object)
  ldotdotdot <- list(...)
  nl <- names(ldotdotdot)
  naAE <- names(attributes(AEnew))
  nl.match <- nl[nl %in% naAE]
  attributes(AEnew)[nl.match] <- ldotdotdot[nl.match]
  AEnew
}


AEdotplot.formula <- function(xr, groups=NULL, data=NULL,
                              sortbyRelativeRisk=TRUE,
                              ...,
                              ## sub=list(deparse(this.call[c(1,2,
                              ##   match(c("groups","data"),
                              ##         names(this.call)))],
                              sub=list(deparse(this.call[1:4],
                                width.cutoff=500), cex=.7)) {
  if (is.null(data))
    stop("'data=' must be specified.", call.=FALSE)

  dsg <- deparse(substitute(groups))
  if (dsg == "NULL")
    stop("'groups=' must be specified.", call.=FALSE)

  sdsg <- strsplit(dsg, "$", fixed=TRUE)[[1]]
  if (length(sdsg == 1)) {
    groups <- as.formula(paste("~", sdsg))
    TRT.name <- sdsg
  } else {
    TRT.name <- sdsg[2]
  }
  if (length(groups[[2]]) > 1)
    stop("'groups=' must be only the Treatment factor.", call.=FALSE)

  ## lPF <- latticeParseFormula(PREF ~ SAE/SN | OrgSys, groups=RAND, data=xae)
  lPF <- latticeParseFormula(xr, groups=groups, data=data)
  AE.name <- lPF$left.name

  strsplit.right.name <- strsplit(lPF$right.name, "/")[[1]]
  if (length(strsplit.right.name) != 2)
    stop("formula must have the form 'AE ~ nAE / nTRT'.", call.=FALSE)
  nAE.name  <- strsplit.right.name[1]
  nTRT.name <- strsplit.right.name[2]

  condition.name <- names(lPF$condition)
  if (length(condition.name) > 1)
    stop("At most one condition may be specified.", call.=FALSE)

  data.names <- names(data)

  formula.names <- c(PREF=AE.name,
                     SAE=nAE.name,
                     SN=nTRT.name,
                     RAND=TRT.name,
                     condition=condition.name)
  matched.names <- match(formula.names, data.names, 0)
  newdata <- data[matched.names]
  names(newdata) <- names(formula.names)

##  this.call <- sys.call(1)
  this.call <- match.call()
  names(this.call)[2] <- ""
  this.call[[1]] <- as.name("AEdotplot")  ## drop the ".formula"

  if (is.null(condition.name)) {             ## formula
    conditionName <- list(...)$conditionName ## argument
    if (is.null(conditionName)) {
      conditionName <- deparse(substitute(data))
      AEdotplot(newdata, ...,
                conditionName=conditionName,
                sortbyRelativeRisk=sortbyRelativeRisk,
                sub=sub)
    }
    else
      AEdotplot(newdata, ...,
                sortbyRelativeRisk=sortbyRelativeRisk,
                sub=sub)
  }
  else
    AEdotplot(newdata, ...,
              conditionVariable=newdata$condition,
              sortbyRelativeRisk=sortbyRelativeRisk,
              sub=sub)
}
