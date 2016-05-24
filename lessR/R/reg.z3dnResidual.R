.reg3dnResidual <-
function(lm.out, pdf=FALSE, pdf.width=5, pdf.height=5, manage.gr=FALSE, ...) {

  nm <- all.vars(lm.out$terms)  # names of vars in the model
  n.vars <- length(nm)
  n.keep <- nrow(lm.out$model)


  # pdf graphics option
  if (pdf) {
    pdf.file <- "RegResiduals.pdf"
    pdf(file=pdf.file, width=pdf.width, height=pdf.height)
  }

  # keep track of the plot in this routine
  plt.i <- 0L
  plt.title  <- character(length=0)

  plt.i <- plt.i + 1L
  plt.title[plt.i] <- "Distribution of Residuals"

  # frequency distribution of residuals
  .dn.main(lm.out$residuals, 
    bw="nrd0", type="both",
    bin.start=NULL, bin.width=NULL,
    col.fill=getOption("col.fill.pt"),
    col.bg=getOption("col.bg"), col.grid=getOption("col.grid"),
    col.box=getOption("col.box"),
    col.nrm="gray40", col.gen="gray40",
    col.fill.nrm="transparent", col.fill.gen="transparent",
    cex.axis=.85, col.axis="gray30", rotate.values=0, offset=0.5, 
    x.pt=NULL, xlab="Residuals",
    main="", y.axis=FALSE, 
    x.min=NULL, x.max=NULL, band=FALSE, quiet=TRUE)

  if (pdf) {
    dev.off()
    .showfile(pdf.file, "residuals plot")
  }

  invisible(list(i=plt.i, ttl=plt.title))

}
