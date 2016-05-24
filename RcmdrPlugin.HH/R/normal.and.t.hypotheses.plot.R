normal.and.t.hypotheses.plot <- function() {
  initializeDialog(title=gettextRcmdr("Normal and t Distributions"))
  muVar <- tclVar("")
  muEntry <- tkentry(top, width="6", textvariable=muVar)
  sigmaVar <- tclVar("")
  sigmaEntry <- tkentry(top, width="6", textvariable=sigmaVar)
  degfreeVar <- tclVar("")
  degfreeEntry <- tkentry(top, width="6", textvariable=degfreeVar)
  nVar <- tclVar("")
  nEntry <- tkentry(top, width="6", textvariable=nVar)
  leftAlphaVar <- tclVar("")
  leftAlphaEntry <- tkentry(top, width="6", textvariable=leftAlphaVar)
  rightAlphaVar <- tclVar("")
  rightAlphaEntry <- tkentry(top, width="6", textvariable=rightAlphaVar)
  obsValueVar <- tclVar("")
  obsValueEntry <- tkentry(top, width="6", textvariable=obsValueVar)
  muAltVar <- tclVar("")
  muAltEntry <- tkentry(top, width="6", textvariable=muAltVar)
  ymaxVar <- tclVar("")
  ymaxEntry <- tkentry(top, width="6", textvariable=ymaxVar)

  radioButtons(name="hypoth",
               buttons=c("Hypoth", "Conf"),
               labels=gettextRcmdr(c("Hypothesis Test", "Confidence Interval")),
               title=gettextRcmdr("Hypothesis or Confidence"))

  onOK <- function() {
    closeDialog()
    mu <- as.numeric(tclvalue(muVar))
    if (is.na(mu)) mu <- 0
    sigma <- as.numeric(tclvalue(sigmaVar))
    if (!is.na(sigma) && sigma <= 0) {
      errorCondition(recall=normal.and.t.hypotheses.plot,
                     message=gettextRcmdr("Standard deviation must be positive."))
            return()
            }
    degfree <- as.numeric(tclvalue(degfreeVar))
    if (!is.na(degfree) && degfree <= 0) {
      errorCondition(recall=normal.and.t.hypotheses.plot,
                     message=gettextRcmdr("Degrees of freedom must be positive."))
            return()
            }
    n <- as.numeric(tclvalue(nVar))
    if (!is.na(n) && n <= 0) {
      errorCondition(recall=normal.and.t.hypotheses.plot,
                     message=gettextRcmdr("Sample size must be positive."))
            return()
            }
    leftAlpha <- as.numeric(tclvalue(leftAlphaVar))
    if (!is.na(leftAlpha) && (leftAlpha < 0 || 1 < leftAlpha)) {
      errorCondition(recall=normal.and.t.hypotheses.plot,
                     message=gettextRcmdr("leftAlpha must be between 0 and 1."))
            return()
            }
    rightAlpha <- as.numeric(tclvalue(rightAlphaVar))
    if (!is.na(rightAlpha) && (rightAlpha < 0 || 1 < rightAlpha)) {
      errorCondition(recall=normal.and.t.hypotheses.plot,
                     message=gettextRcmdr("rightAlpha must be between 0 and 1."))
            return()
            }
    if (!(is.na(leftAlpha) || is.na(rightAlpha))&& leftAlpha + rightAlpha > 1) {
      errorCondition(recall=normal.and.t.hypotheses.plot,
                     message=gettextRcmdr("left Alpha + rightAlpha must be less than 1."))
            return()
            }

    hypoth <- as.character(tclvalue(hypothVariable))
    obsValue <- as.numeric(tclvalue(obsValueVar))
    if (hypoth=="Conf" && is.na(obsValue)) {
      errorCondition(recall=normal.and.t.hypotheses.plot,
                     message=gettextRcmdr("Confidence Intervals require an observed mean."))
            return()
            }
    
    muAlt <- as.numeric(tclvalue(muAltVar))
    ymax <- as.numeric(tclvalue(ymaxVar))
    if (!is.na(ymax) && ymax <= 0) {
      errorCondition(recall=normal.and.t.hypotheses.plot,
                     message=gettextRcmdr("ymax must be positive."))
            return()
            }
    
    command <- paste('normal.and.t.dist.result <- ',
                     'normal.and.t.dist(',
                     'std.dev=', sigma,
                     ', n=', n,
                     ', deg.freedom=', degfree,
                     ', mu.H0=', mu,
                     ', xmin=', NA,
                     ', xmax=', NA,
                     ', gxbar.min=', NA,
                     ', gxbar.max=', ymax,
                     ', Use.alpha.right=', !is.na(rightAlpha),
                     ', alpha.right=', rightAlpha,
                     ', Use.alpha.left=', !is.na(leftAlpha),
                     ', alpha.left=', leftAlpha,
                     ', Use.mu.H1=', !is.na(muAlt),
                     ', mu.H1=', muAlt,
                     ', Use.obs.mean=', !is.na(obsValue),
                     ', obs.mean=', obsValue,
                     ', hypoth.or.conf="', hypoth, '"',
                     ', col.mean=', '"lime green"',
                     ', col.alpha=', '"blue"',
                     ', col.beta=', '"red"',
                     ', col.conf=', '"pale green"',
                     ', col.conf.arrow=', '"dark green"',
                     ', col.mean.label=', '"lime green"',
                     ', col.alpha.label=', '"blue"',
                     ', col.beta.label=', '"red"',
                     ', col.conf.label=', '"dark green"',
                     ', cex.crit=', 1.2,
                     ')', sep='') 
        doItAndPrint(command)
      
      tkfocus(CommanderWindow())
    }
  OKCancelHelp(helpSubject=if (is.na(as.numeric(tclvalue(degfreeVar)))) "dnorm" else "dt")
  tkgrid(tklabel(top, text=gettextRcmdr("mu (mean)")), muEntry, sticky="e")
  tkgrid(tklabel(top, text=gettextRcmdr("standard deviation (sigma or s)")), sigmaEntry, sticky="e")
  tkgrid(tklabel(top, text=gettextRcmdr("df (degrees of freedom)")), degfreeEntry, sticky="e")
  tkgrid(tklabel(top, text=gettextRcmdr("n (sample size)")), nEntry, sticky="e")
  tkgrid(tklabel(top, text=gettextRcmdr("left alpha")), leftAlphaEntry, sticky="e")
  tkgrid(tklabel(top, text=gettextRcmdr("right alpha")), rightAlphaEntry, sticky="e")
  tkgrid(tklabel(top, text=gettextRcmdr("Observed Value")), obsValueEntry, sticky="e")
  tkgrid(tklabel(top, text=gettextRcmdr("mu (Alternate Hypothesis)")), muAltEntry, sticky="e")
  tkgrid(tklabel(top, text=gettextRcmdr("ymax (right-hand side)")), ymaxEntry, sticky="e")
  tkgrid.configure(muEntry, sticky="w")
  tkgrid.configure(sigmaEntry, sticky="w")
  tkgrid.configure(degfreeEntry, sticky="w")
  tkgrid.configure(nEntry, sticky="w")
  tkgrid.configure(leftAlphaEntry, sticky="w")
  tkgrid.configure(rightAlphaEntry, sticky="w")
  tkgrid.configure(obsValueEntry, sticky="w")
  tkgrid.configure(muAltEntry, sticky="w")
  tkgrid.configure(ymaxEntry, sticky="w")
  tkgrid(hypothFrame, sticky="w")
  tkgrid(buttonsFrame, columnspan=2, sticky="w")
  dialogSuffix(rows=10, columns=2, focus=muEntry)
}

## source("~/HH-R.package/RcmdrPlugin.HH/R/normal.and.t.hypotheses.plot.R")
