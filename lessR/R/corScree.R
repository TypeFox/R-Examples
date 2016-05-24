corScree <- 
function (x=mycor, 
          main=NULL, pdf=FALSE, pdf.width=5, pdf.height=5, ...) {


  # cor matrix:  mycor as class out_all, mycor$cors, or stand-alone matrix
  cor.nm <- deparse(substitute(x))
  .cor.exists(cor.nm)  # see if matrix exists in one of the 3 locations
  if (class(x) == "out_all")
    x <- eval(parse(text=paste(cor.nm, "$cors", sep="")))  # go to $cors 

  
  # extract eigenvectors
  eig <- eigen(x, symmetric=TRUE, only.values=TRUE)
  ev <- eig$values

  # see if graphics are to be managed
  manage.gr <- .graphman()

  # if manage, set up graphics system for 2 windows default
  if (!pdf) {
    if (manage.gr) {
      .graphwin(2)
      dev.set(which=3)
    }
  }
  else { 
    pdf.file <- "Scree.pdf"
    pdf(file=pdf.file, width=pdf.width, height=pdf.height)
  }

  if (getOption("colors") == "gray" || getOption("colors") == "gray.black")
    col.ln <- getOption("col.fill.bar")
  else
    col.ln <- getOption("col.stroke.bar")

  # keep track of generated graphics
  plot.i <- 0
  plot.title  <- character(length=0)

  plot.i <- plot.i + 1
  plot.title[plot.i] <- "Eigenvalues"

  # scree plot
  .lc.main(ev, type=NULL, 
         col.line=col.ln,
         col.area=NULL, col.box="black",
         col.stroke=getOption("col.stroke.pt"),
         col.fill=getOption("col.fill.pt"),
         shape.pts=21, col.grid=getOption("col.grid"),
         col.bg=getOption("col.bg"),
         cex.axis=.85, col.axis="gray30", rotate.values=0, offset=0.5,
         xy.ticks=TRUE, line.width=1.1,
         xlab=NULL, ylab="Eigenvalues", main=main, sub=NULL, cex=NULL,
         x.start=NULL, x.end=NULL, y.start=NULL, y.end=NULL,
         time.start=NULL, time.by=NULL, time.reverse=FALSE,
         center.line="off", quiet=TRUE, ...)

  if (pdf) {
    dev.off()
    .showfile(pdf.file, "scree chart")
  }

  if (!pdf) {
    if (manage.gr) {
      dev.set(which=4) 
    }
  }
  else { 
    pdf.file <- "ScreeDiff.pdf"
    pdf(file=pdf.file, width=pdf.width, height=pdf.height)
  }

  # differences scree plot
  ev.diff <- -diff(ev)

  plot.i <- plot.i + 1
  plot.title[plot.i] <- "Differences of Successive Eigenvalues"

  .lc.main(ev.diff, type=NULL, 
         col.line=col.ln,
         col.area=NULL, col.box="black",
         col.stroke=getOption("col.stroke.pt"),
         col.fill=getOption("col.fill.pt"),
         shape.pts=21, col.grid=getOption("col.grid"),
         col.bg=getOption("col.bg"),
         cex.axis=.85, col.axis="gray30", rotate.values=0, offset=0.5,
         xy.ticks=TRUE, line.width=1.1,
         xlab=NULL, ylab="Differences of Successive Eigenvalues",
         main=main, sub=NULL, cex=NULL,
         x.start=NULL, x.end=NULL, y.start=NULL, y.end=NULL,
         time.start=NULL, time.by=NULL, time.reverse=FALSE,
         center.line="off", quiet=TRUE, ...)
   n.dregs <- ceiling(length(ev.diff)/1.35)  # get bottom sequence of ev differences
   dregs <- numeric(length=n.dregs)
   for (i in 1:n.dregs) dregs[i] <- ev.diff[length(ev.diff)-(i-1)] 
   abline(h=mean(dregs), col="gray50", lwd=2)

  if (pdf) {
    dev.off()
    .showfile(pdf.file, "scree difference chart")
    cat("\n\n")
  }

  cat("\n")
  cat("Eigenvalues of", cor.nm, "\n")
  .dash(15+nchar(cor.nm))
  for (i in 1:length(ev)) cat(round(ev[i],3), " ")
  cat("\n")

  cat("\n")
  cat("Differences of Successive Eigenvalues of", cor.nm, "\n")
  .dash(41+nchar(cor.nm))
  for (i in 1:length(ev.diff)) cat(round(ev.diff[i],3), " ")
  cat("\n")

  if (is.null(options()$knitr.in.progress))
    .plotList(plot.i, plot.title)

  cat("\n")

}
