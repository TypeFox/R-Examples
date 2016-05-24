#' achievable.se.sp(): target sensitivity and specificity for a retrospective study.
#'
#' For a retrospective study design, given a prevalence value,
#' produce a plot displaying the achievable contours of
#' either predictive values or NNT values.
#' The calculation uses the "contra-Bayes" theorem, sesp.from.pv.
#'
#' @param the.prev Prevalence (prior probability)
#' @param axes Should the axes be predictive values ("pv") or NNT values? Default is "pv".
#' @param sesp.seq Sequence of values at which the sensitivity and specificity will be explored.
#' @param drawNNTaxes (default=TRUE) Option for tweaking the plot.
#' @param drawPVaxes (default=FALSE) Option for tweaking the plot.
#' @param drawArrows (default=TRUE) Arrow option; deprecated.
#' @param drawTable (default=TRUE) Option for tweaking the plot.
#' @param latexTable (default=TRUE) Option for tweaking the plot.
#' @param placePointLabels (default=TRUE)  Write A, B, C, ... in circles where sensitivity=specificity.
#' @param cexText (default=0.5) Option for tweaking the plot.
#' @param cexSubtitle (default=0.5) Option for tweaking the plot.
#' @param cexTitle (default=0.7) Option for tweaking the plot.
#' @param y0arrow (default=0.25) Arrow option; deprecated.
#' @param lwdArrow (default=1) Arrow option; deprecated.
#' @param ltyArrow (default=2) Arrow option; deprecated.
#' @param title (default=FALSE) Option for tweaking the plot.
#' @param  mtext (default=FALSE) Option for tweaking the plot.
#' @param  contours (default=TRUE) Option for tweaking the plot.
#' @param ... Options to pass to plot.default()
#' @return The predictive values when sensitivity equals specificity:  sesp.to.pv(cbind(sesp.seq,sesp.seq), prev=the.prev))

achievable.se.sp = function(the.prev = 0.5,
                            axes=c("pv", "NNT"),
                            sesp.seq = seq(0.5, 1, 0.10),
                            drawNNTaxes=TRUE,
                            drawPVaxes=FALSE,
                            drawArrows=TRUE,
                            drawTable=TRUE,
                            latexTable=TRUE,
                            placePointLabels=TRUE, ### A, B, C, ...
                            cexText=0.5,
                            cexSubtitle=0.5,
                            cexTitle=0.7,
                            y0arrow=0.25,
                            lwdArrow=1,
                            ltyArrow=2,
                            title=FALSE, mtext=FALSE, contours=TRUE,
                            ...) {
  if(axes[1] == "NNT")
    diagFun = function(se, sp, prev){
      pv = sesp.to.pv(se, sp, prev=the.prev)
      print(pv)
      NNT = pv.to.NNT(pv=pv)
      print(NNT)
      return(NNT)
    }
  else
    diagFun = function(se,sp,prev)
      sesp.to.pv(se, sp, prev = the.prev)
  diagonal.values = diagFun(se=sesp.seq, sp=sesp.seq,
                            prev=the.prev)
  print(diagonal.values)
  #browser()
  theCall = as.list(match.call())
  ylim = unlist(as.list(theCall[["ylim"]])[-1])
  if(is.null(ylim))   ylim = c(max(0, 1 - (1-min(diagonal.values[ , 2]))), 1)
  xlim = unlist(as.list(theCall[["xlim"]])[-1])
  if(is.null(xlim))  xlim = 0:1
  print(xlim)
  print(ylim)
  plot(diagonal.values, type="l", lwd=2,
       xlab=if(axes[1]=="pv") "positive predictive value"
       else expression(bolditalic(NNT[Pos])),
       ylab=if(axes[1]=="pv") "negative predictive value"
       else expression(bolditalic(NNT[Neg])),
       xaxs="i", yaxs="i",
       ...)

  if(contours) {
    for(se in sesp.seq)
      lines((diagFun(se=rep(se, length(sesp.seq)), sp=sesp.seq, prev=the.prev)),
            lwd=2, col="blue")
    for(sp in sesp.seq)
      lines((diagFun(se=sesp.seq, sp=rep(sp, length(sesp.seq)), prev=the.prev)),
            lwd=2, col="red")
  }
  if(placePointLabels) {
    solidcircle = 16
    points(xpd=NA,
           diagonal.values, col="black", pch=solidcircle, cex=2.7)
    graphics::text(xpd=NA, diagonal.values, labels=LETTERS[1:nrow(diagonal.values)],
#        vfont=c("serif", "bold"),  # Beware of shinyjs::text
         font=4, family="serif",
         col="white", cex=0.9)
  }
  # This commented-out section was written for a publication diagram.
  #  if(drawArrows){
  #     x0positions = seq(xlim[1]+0.1*(xlim[2]-xlim[1]),
  #                       xlim[2]-0.1*(xlim[2]-xlim[1]), along=sesp.seq)
  #     arrows(x0=x0positions,
  #            y0=rep( par()$usr[3] + (par()$usr[4]-par()$usr[3])*y0arrow,	length(sesp.seq)),
  #            x1=diagonal.values[1, ],
  #            y1=diagonal.values[2, ],
  #            lwd=lwdArrow, lty=ltyArrow, length=0.1, angle=10)
  #   }
  #
  #     text(pos=1, x0positions, cex=cexText,
  #          rep( par()$usr[3] + (par()$usr[4]-par()$usr[3])*(y0arrow*0.8),
  #               length(sesp.seq)),
  #          sesp.seq %&% "\n" %&%
  #            sesp.seq %&% "\n" %&%
  #            signif(digits=3, diagonal.values[1, ] ) %&%
  #            "\n" %&% signif(digits=3, diagonal.values[2, ])
  #          %&% "\n" %&% NNTpos
  #          %&% "\n" %&% NNTneg
  #     )
  #   }

if(axes[1]=="NNT") {
      NNTpos = signif(digits=3, diagonal.values[ ,1 ])
      NNTneg = signif(digits=3, diagonal.values[ ,2 ])
      PV = NNT.to.pv(NNTpos, NNTneg)
    }
    else {
      NNTpos = signif(digits=3, NNT.from.pv(pv=diagonal.values)$NNTpos)
      NNTneg = signif(digits=3, NNT.from.pv(pv=diagonal.values)$NNTneg)
      PV = diagonal.values
    }
  cat("PV\n"); print(PV)
  if(drawTable){
    graphics::text(pos=2, -0.03, cex=cexText,
         -0.018 +rep( par()$usr[3] + (par()$usr[4]-par()$usr[3])*(y0arrow*0.8),
                      length(sesp.seq)),
         "se=\n" %&%
           "sp=\nppv=\nnpv=\nNNTpos=\nNNTneg="
    )
  }
  if(latexTable)
    xtable::xtable(digits=3, t(data.frame(
      sensitivity=sesp.seq, specificity=sesp.seq,
      PPV=PV[ , "ppv"],
      NPV=PV[ , "npv"],
      NNTpos=NNTpos,
      NNTneg=NNTneg
    )))
  if(title)
    title(paste("prevalence = ", the.prev), cex=cexTitle)
  if(mtext)
    mtext(paste(
      "orange (steeper) lines: fixed sensitivity", "\n",
      "    blue (gentler) lines: fixed specificity", collapse=" "
    ), side=3, cex=cexSubtitle)
  if(drawNNTaxes)  {
      #1=below, 2=left, 3=above and 4=right.
      yaxp = par()$yaxp
      verticalTickMarks = seq(yaxp[1], yaxp[2], length=yaxp[3] + 1) [-(yaxp[3]+1)][-1]
      xaxp = par()$xaxp
      horizontalTickMarks = seq(xaxp[1], xaxp[2], length=xaxp[3] + 1) [-(xaxp[3]+1)][-1]
      NNT_ = NNT.from.pv(ppv = horizontalTickMarks, npv = verticalTickMarks)
      NNTpos = NNT_$"NNTpos"
      NNTneg = NNT_$"NNTneg"
      mtext(expression(NNTpos %->% phantom(0)), at=xaxp[1], line=1, side=3)
      mtext(expression(NNTneg %->% phantom(0)), at=yaxp[1], line=1, side=4)
      axis(3, at=horizontalTickMarks, labels=round(NNTpos, 2))
      axis(4, at=verticalTickMarks, labels=round(NNTneg, 2))
  }
  return(invisible(sesp.to.pv(cbind(sesp.seq,sesp.seq), prev=the.prev)))
}
# END  achievable.se.sp()


