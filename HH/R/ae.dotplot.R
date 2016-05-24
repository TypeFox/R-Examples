## variable names in the input data.frame to ae.dotplot.long
##
## RAND   treatment as randomized
## PREF   adverse event symptom name
## SN     number of patients in treatment group
## SAE    number of patients  in each group for whom the event PREF was observed
##
## Input sort order should be PREF/RAND.
## logrelrisk orders the data to assure itself of that order.

## Calculate the percent, the relative risk, the log relative risk,
## and the confidence intervals.
## Make PREF an ordered factor, sorted by the relative risk.
logrelrisk <- function(ae,
                       A.name=levels(ae$RAND)[1],
                       B.name=levels(ae$RAND)[2],
                       crit.value=1.96) {
  ae$PCT <- 100 * ae$SAE / ae$SN ## percent of patients

  ## Calculation of relative risk assumes both A and B PREF are in the
  ## same order.  Assignment of relative risk to the data.frame assumes
  ## that A and B for each PREF are in consecutive positions.
  ae <- ae[with(ae, order(PREF, RAND)),]
  AA <- ae$PREF[ae$RAND==A.name]
  BB <- ae$PREF[ae$RAND==B.name]
  if (length(AA) != length(BB) || !all(AA == BB))
    stop("Data problem: at least one AE is missing an A treatment or a B treatment.", call.=FALSE)

  tmp <-  ## sample relative risk
    ae$PCT[ae$RAND==B.name] /
      ae$PCT[ae$RAND==A.name]

  ## sort by relrisk
  tmp.order <- order(tmp) ## sort by logrelrisk
  ##  ordered(ae$PREF) <- as.character(unique(ae$PREF))[tmp.order]  ## S-Plus only
  ae$PREF <- ordered(ae$PREF, levels=(ae$PREF[ae$RAND==A.name])[tmp.order]) ## R and S-Plus

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

  ae
}


panel.ae.leftplot <- function(x, y, groups, col.AB, ...) {
  panel.abline(h=y, lty=2, lwd=.3, col=1)
  panel.abline(v=0, lty=3, lwd=.3, col=1)
  panel.superpose(x, y, groups=groups, col=col.AB, ...)
}

panel.ae.rightplot <- function(x, y, ..., lwd=6, lower, upper, cex=.7) {
  if.R(r={}, s={panel.segments <- segments; panel.points <- points})
  panel.abline(v=0, lty=3, lwd=.3)
  panel.abline(h=y, lty=2, lwd=.3, col=1)
  panel.segments(lower, y, upper, y, lwd=2)
  panel.points(lower, y, pch=3, col=1, cex=.4)
  panel.xyplot(x, y, ..., col=1, cex=cex)
  panel.points(upper, y, pch=3, col=1, cex=.4)
}
## assignInNamespace("panel.ae.leftplot",  panel.ae.leftplot,  "HH")
## assignInNamespace("panel.ae.rightplot", panel.ae.rightplot, "HH")

panel.ae.dotplot <- function(x, y, groups, ..., col.AB, pch.AB, lower, upper) {
  panel.num <-
    ## if.R(s=get("cell", frame=sys.parent()),
    ##                 r=
    panel.number()
  ##          )
  if (panel.num==1)
    panel.ae.leftplot(x, y, groups=groups, col.AB=col.AB, pch=pch.AB, ...)
  if (panel.num==2)
    panel.ae.rightplot(x, y, ..., lwd=6, pch=16,
                       lower=lower, upper=upper)
}


ae.dotplot.long <- if.R(r={
  function(xr,
           A.name=levels(xr$RAND)[1],
           B.name=levels(xr$RAND)[2],
           col.AB=c("red","blue"), pch.AB=c(16,17),
           main.title = paste("Most Frequent On-Therapy Adverse Events",
                              "Sorted by Relative Risk"),
           main.cex=1,
           cex.AB.points=NULL, cex.AB.y.scale=.6,
           position.left= c(0,   0, .70, 1.), ## ignored in R
           position.right=c(.61, 0, .98, 1.), ## ignored in R
           key.y=-.2, CI.percent=95) {

    if (is.null(xr$logrelrisk))
      stop("Variable 'logrelrisk' missing.\nPlease use the logrelisk() function before using ae.dotplot().")
    result <-
      dotplot(PREF ~ PCT + logrelrisk,
              groups=xr$RAND, data=xr, outer=TRUE,
              lower=xr$logrelriskCI.lower,
              upper=xr$logrelriskCI.upper,
              panel=panel.ae.dotplot,
              scales=list(
                x=list(
                  relation="free",
                  at=list(
                    seq(0, 100, 10),
                    log(c(.125, .25, .5, 1, 2, 4, 8, 16, 32))
                    ),
                  labels=list(
                    seq(0, 100, 10),
                    c(".125", "", ".5", 1, 2, 4, 8, 16, 32)
                    ),
                  limits=list(
                    range(xr$PCT),
                    range(xr$logrelriskCI.lower, xr$logrelriskCI.upper))
                  ),
                y=list(cex=cex.AB.y.scale)),
              A.name=A.name, B.name=B.name,
              col.AB=col.AB, pch.AB=pch.AB,
              cex.AB.points=cex.AB.points,
              cex.AB.y.scale=cex.AB.y.scale,
              main=list(main.title, cex=main.cex),
              xlab=NULL,
              between=list(x=1),
              key=list(y = key.y, x=.15,
                points = list(col=col.AB, pch=pch.AB),
                text = list(c(A.name, B.name), col=col.AB, cex=.9),
                columns = 2,
                between=.5,
                space="bottom")
              )
    if.R(r=result$condlevels[[1]] <-
         c("Percent",
           paste("Relative Risk with ", CI.percent, "% CI", sep=""))
         )
    result
  }
},s={

  function(xr,
           A.name=levels(xr$RAND)[1],
           B.name=levels(xr$RAND)[2],
           col.AB=c("red","blue"), pch.AB=c(16,17),
           main.title="Most Frequent On-Therapy Adverse Events Sorted by Relative Risk",
           main.cex=1,
           cex.AB.points=NULL, cex.AB.y.scale=.6,
           position.left= c(0,   0, .70, 1.),
           position.right=c(.61, 0, .98, 1.),
           key.y=-.2, CI.percent=95) {

    ae.key <- list(y = -.3,
                   points = list(col=0),
                   text = list(""),
                   columns = 2,
                   space="bottom")

    ae.main <- list(" ", cex=1.5)

    ## construct left panel
    left.plot <- dotplot(PREF ~ PCT, data=xr, groups = RAND,
                         col=col.AB, pch=pch.AB,
                         panel = panel.ae.leftplot,
                         xlab = "Percent",
                         scales=list(y=list(cex=cex.AB.y.scale)),
                         main = ae.main,
                         cex=cex.AB.points,
                         key = ae.key
                         )

    ## construct right panel
    if (is.null(xr$logrelrisk))
      stop("Variable 'logrelrisk' missing.\nPlease use the logrelisk() function before using ae.dotplot().")
    right.plot <- dotplot(PREF ~ logrelrisk, data=xr, pch=16,
                          lower=xr$logrelriskCI.lower,
                          upper=xr$logrelriskCI.upper,
                          panel=panel.ae.rightplot,
                          xlab = paste("Relative Risk with ",
                            CI.percent, "% CI", sep=""),
                          scales=list(
                            x=list(at=log(c(.125, .25, .5, 1, 2, 4, 8, 16, 32)),
                              labels=c(".125", "", ".5", 1, 2, 4, 8, 16, 32)),
                            y=list(cex=cex.AB.y.scale)
                            ),
                          xlim=c(-2.3,4),
                          main = ae.main,
                          key = ae.key
                          )
    ##  right.plot
    if.R(s=
         right.plot$scales$y$labels[] <- " " ## suppress the AE names on the right panel
         ,r=
         right.plot$y.limits[] <- ""
         )
    ##  right.plot

    ## print both plots on current device
    print(left.plot, position=position.left,  more=TRUE)
    print(right.plot, position=position.right, more=FALSE)
    title(main.title, cex=main.cex)
    key(y = key.y, x=.15,
        points = list(col=col.AB, pch=pch.AB),
        text = list(c(A.name, B.name), col=col.AB, cex=.9),
        columns = 2,
        between=.5,
        space="top")
    invisible(list(left.plot, right.plot))
  }
}
                   )

aeReshapeToLong <- function(aewide) {
  ## this code works in both R and S-plus
  aewide$PREF<- ordered(aewide$Event, rev(aewide$Event))
  aewide$Event <- NULL
  if (! "PCT.A" %in% names(aewide))
    aewide$PCT.A <- aewide$AE.A / aewide$N.A
  if (! "PCT.B" %in% names(aewide))
    aewide$PCT.B <- aewide$AE.B / aewide$N.B

  aewide.names <- names(aewide)
  aewide.namesA <- match(c("N.A", "AE.A", "PCT.A"), aewide.names)
  aewide.namesB <- match(c("N.B", "AE.B", "PCT.B"), aewide.names)
  aewideA <- cbind(aewide, RAND="A")
  aewideB <- cbind(aewide, RAND="B")
  names(aewideA)[aewide.namesA] <- c("SN","SAE","PCT")
  names(aewideB)[aewide.namesB] <- c("SN","SAE","PCT")
  aelong <- rbind(aewideA[-aewide.namesB], aewideB[-aewide.namesA])

  n <- nrow(aewide)
  aelong[as.vector(outer(c(0,n), 1:n, "+")), ]
}

## if (FALSE)  ## aeReshapeToLongR uses reshape and works only in R.
##   if.R(r=
##        aeReshapeToLongR <- function(aewide) {
##          aewide$PREF<- ordered(aewide$Event, rev(aewide$Event))
##          aelong$Event <- NULL
##          if (! "PCT.A" %in% names(aewide))
##            aewide$PCT.A <- aewide$AE.A / aewide$N.A
##          if (! "PCT.B" %in% names(aewide))
##            aewide$PCT.B <- aewide$AE.B / aewide$N.B
##          aelong <- reshape(aewide, v.names=c("SN","SAE","PCT"),
##                            varying=list(
##                              c("N.A", "N.B"),
##                              c("AE.A", "AE.B"),
##                              c("PCT.A", "PCT.B")),
##                            idvar="PREF",
##                            timevar="RAND",
##                            times=factor(c("A","B")),
##                            direction="long")
##          n <- nrow(aewide)
##          aelong[as.vector(outer(c(0,n), 1:n, "+")), ]
##        }
##        ,s={})
##

ae.dotplot <- function(ae, ...) {
  if (!("PREF" %in% names(ae)))
    ae <- aeReshapeToLong(ae)
  result <- ae.dotplot.long(ae, ...)
  if.R(r=result,
       s=invisible(result))
}

## ae.dotplot(aeanonymr,
##            A.name="TREATMENT A (N=216)",
##            B.name="TREATMENT B (N=431)")
