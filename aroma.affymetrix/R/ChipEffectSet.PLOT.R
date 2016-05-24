# @author "HB, EP"
setMethodS3("plotBoxplot", "ChipEffectSet", function(this, type, transform=NULL, ..., main=NULL, ylab=NULL) {
  # Argument 'type':
  # Further validation in 'boxplotStats()':
  type <- Arguments$getCharacter(type, length=c(1,1));



  if (is.null(main))
    main <- type;

  if (is.null(ylab)) {
    if (type == "theta") {
      if (is.function(transform) && identical(transform, log2)) {
        ylab <- expression(log[2](theta));
      } else {
        ylab <- expression(theta);
      }
    } else if (type == "RLE") {
      ylab <- expression(M == log[2](theta/bar(theta)));
    } else if (type == "NUSE") {
      ylab <- expression(paste(rho/bar(rho), "  [", rho==log[2](sigma), "]"));
    } else {
      ylab <- type;
    }
  }


  stats <- boxplotStats(this, type=type, transform=transform, ...);

  plotBoxplotStats(stats, main=main, ylab=ylab, ...);

  invisible(stats);
})


##########################################################################
# HISTORY:
# 2008-02-25
# o Created from EPs code.
##########################################################################
