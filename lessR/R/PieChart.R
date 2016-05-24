PieChart <-
function(x, data=mydata, 
         color.fill=NULL, color.low=NULL, color.hi=NULL,
         colors=c("rainbow", "terrain", "heat"),
         color.random=FALSE, main=NULL, cex=1, cex.main=1,
         quiet=getOption("quiet"),
         pdf.file=NULL, pdf.width=5, pdf.height=5, ...) {


  if (missing(colors)) 
    colors <- getOption("colors")
  else
    colors <- match.arg(colors)

  dots <- list(...)  # check for deprecated parameters
  if (length(dots) > 0) {
    for (i in 1:length(dots)) {
      old.nm <- c("col.fill", "col.low", "col.hi")
      if (names(dots)[i] %in% old.nm) {
        cat("\n"); stop(call.=FALSE, "\n","------\n",
          "options that began with the abbreviation  col  now begin with  ",
          "color \n\n")
      }
    }
  }

  if (!is.null(pdf.file))
    if (!grepl(".pdf", pdf.file)) pdf.file <- paste(pdf.file, ".pdf", sep="")

  x.name <- deparse(substitute(x)) 
  options(xname = x.name)

  # get data frame name
  dname <- deparse(substitute(data))
  options(dname = dname)

  # get conditions and check for data existing
  xs <- .xstatus(x.name, dname, quiet)
  is.frml <- xs$ifr
  from.data <- xs$fd
  in.global <- xs$ig 

  # see if variable exists in the data frame, if x not in Global Env or function call 
  if (!missing(x) && !in.global)  .xcheck(x.name, dname, data)

  if (!in.global) x.call <- eval(substitute(data$x))
  else {  # vars that are function names get assigned to global
    x.call <- x
    if (is.function(x.call)) x.call <- eval(substitute(data$x))
  }

  # set up graphics system
  .opendev(pdf.file, pdf.width, pdf.height)

  #orig.params <- par(no.readonly=TRUE)
  #on.exit(par(orig.params))

  .pc.main(x.call, 
       color.random, color.fill, color.low, color.hi,
       colors, cex, cex.main, quiet, main, 
       pdf.file, pdf.width, pdf.height, ...)

  # terminate pdf graphics system
  if (!is.null(pdf.file)) {
    dev.off()
    .showfile(pdf.file, "pie chart")
  }

}
