lmplot <- function(lm.object, ..., main=NULL,
                   col=trellis.par.get("superpose.symbol")$col[1:2],
                   ylim=NULL) {
  if (!missing(col)) {
    old.col <- trellis.par.get("superpose.symbol")$col
    col <- rep(col, length=7)
    trellis.par.set(superpose.symbol=list(col=col))
    on.exit(trellis.par.set(superpose.symbol=list(col=old.col)))
  }

  A <- residVSfitted(lm.object, pch=c(25,24), fill=col[1:2])
  B <- scaleLocation(lm.object, pch=c(25,24), fill=col[1:2])
  BA <- c("Scale-Location"=B,
          "Residuals vs Fitted"=
          update(A, scales=list(y=list(at=-100, alternating=3))),
          layout=c(1,2))
  BAu <-
    update(BA,
           ylab=c(B$ylab, A$ylab),
           ylab.right=c(B$ylab.right, A$ylab.right),
           xlab.top=NULL,
           main=main,
           between=list(y=1),
           scales=list(y=list(at=-100)),
           par.settings=list(layout.widths=list(ylab.axis.padding=3, axis.key.padding=2, ylab.right=6))
           )

  BAu <-
  update(BAu,
         strip=strip.custom(factor.levels=c(
                              expression(sqrt(""~abs(""~Residuals~"")~"") ~~~ widetilde("  ") ~~~ "Fitted Values"),
                              expression(Residuals ~~~ widetilde("  ") ~~~ "Fitted Values"))
           ),
         par.strip.text=list(lines=1.3))
  if (!is.null(ylim))
    BAu <- update(BAu,
                  scales=list(y=list(limits=list(c(0,max(abs(ylim))), ylim))))

  print(BAu,
        split=c(1,1,2,1), more=TRUE)


  C <- diagQQ(lm.object, col=col[1], main=main)

  Cu <-
  update(c("Normal Q-Q"=C), xlab.top=NULL, strip=TRUE,
               par.strip.text=list(lines=1.3))
  if (!is.null(ylim))
    Cu <- update(Cu, scales=list(y=list(limits=ylim)))

  print(Cu,
        ## split=c(2,1,2,2),
        position=c(.53, .54, 1, 1),  ## .54 is function of device and size
        more=TRUE)


  D <- diagplot5new(lm.object, col=col[1], main=main)
  Du <-
  update(D, xlab.top=NULL,
         strip=strip.custom(factor.levels=c(
                              expression(hat(y)-bar(y)),
                              expression( y -hat(y)))),
         ylab="Centered Fitted Values",
         ylab.right="Residuals",
         par.strip.text=list(lines=1.3))
  if (!is.null(ylim))
    Du <- update(Du, scales=list(y=list(limits=ylim)))

  print(Du,
        ## split=c(2,2,2,2),
        position=c(.53, 0, 1, .57),  ## .57 is function of device and size
        more=FALSE)
  ## the .54 and .57 work nicely with the default quartz window on Mac OS X.

  invisible(list(BAu, Cu, Du))
}
