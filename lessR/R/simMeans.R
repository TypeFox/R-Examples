simMeans <- 
function(ns, n, mu=0, sigma=1, ylim.bound=NULL, 
         show.title=TRUE, show.data=TRUE, max.data=10, 
         color.grid="grey90", pause=FALSE,
         sort=NULL, set.mu=FALSE, digits.d=2,
         main=NULL, pdf.file=NULL, pdf.width=5, pdf.height=5, ...) {


  if (missing(ns)) {
    cat("\n"); stop(call.=FALSE, "\n","------\n",
      "Specify the number of samples, each of a given size, with:  ns\n\n")
  }

  if (missing(n)) {
    cat("\n"); stop(call.=FALSE, "\n","------\n",
      "Specify the size of each sample, the number of data values, with:  n\n\n")
  }

  dots <- list(...)  # check for deprecated parameters
  if (length(dots) > 0) {
    for (i in 1:length(dots)) {
      if (substr(names(dots)[i], 1, 4) == "col.") {
        cat("\n"); stop(call.=FALSE, "\n","------\n",
          "options that began with the abbreviation  col  now begin with  ",
          "color \n\n")
      }
    }
  }

  if (sigma < 0) { 
    cat("\n"); stop(call.=FALSE, "\n","------\n",
      "Standard deviation, sigma, cannot be negative.\n\n")
  }

  if (!is.null(pdf.file))
    if (!grepl(".pdf", pdf.file)) pdf.file <- paste(pdf.file, ".pdf", sep="")

  if (is.null(sort)) if (pause) sort <- FALSE else sort <- TRUE

  if (set.mu) {
    mu <- sample(0:100, size=1)
    sigma <- sample(1:25, size=1)
  }
   
  if (!set.mu && !pause) {
    cat("\nPopulation mean, mu:", mu, "\n")
    cat("Pop std dev, sigma :", sigma, "\n")
  }
  cat("\nNumber of samples  :", ns, "\n")
  cat("Size of each sample:", n, "\n")
   
  # data generation
  data.raw <- rnorm(ns*n, mu, sigma)
  data.rep <- matrix(data.raw, ns, n)

  # summary stats
  Ymean <- apply(data.rep, 1, mean)
  Ysd <- apply(data.rep, 1, sd)
  max.Y <- max(Ymean)
  min.Y <- min(Ymean)

  # sort
  if (sort) {
    o <- order(Ymean)
    Ymean <- Ymean[o]
    Ysd <- Ysd[o]
    data.rep <- data.rep[o,]
  }
  else o <- 1:ns

  # keep mu centered with symmetric upper and lower limits on y-axis
  if (is.null(ylim.bound)) {
    l <- min.Y
    u <- max.Y
    max.dev <- max( abs(mu-l), abs(u-mu) )
    l <- mu - max.dev
    u <- mu + max.dev
  }
  else {
    l <- mu - ylim.bound
    u <- mu + ylim.bound
  }


  # plot 

  # set up graphics system
  .opendev(pdf.file, pdf.width, pdf.height)

  orig.params <- par(no.readonly=TRUE)
  par(mar=c(3,3,1.75,2), mgp=c(1.75,.5,0))

  plot(0, type = "n", xlim=c(1,ns), ylim = c(l,u), xlab = "", ylab = "Sample Mean", 
       cex.main=.95, cex.axis=.8)

  # colored plotting area
  usr <- par("usr")
  rect(usr[1], usr[3], usr[2], usr[4], col="ghostwhite", border="black")

  # grid lines
  vx <- pretty(c(usr[1],usr[2]))
  vy <- pretty(c(usr[3],usr[4]))
  abline(v=seq(vx[1],vx[length(vx)],vx[2]-vx[1]), col=color.grid, lwd=.5)
  abline(h=seq(vy[1],vy[length(vy)],vy[2]-vy[1]), col=color.grid, lwd=.5)

  if (!set.mu && !pause) {
    if (show.title) title(main = bquote(paste(mu, "=", .(mu), "  ", sigma, "=", .(sigma), 
        "  ", "  n=", .(n))), cex.main=1)
    mtext(bquote(paste(" ", mu)), side=4, cex=1.5, col="darkslateblue", las=2)
    abline(h=mu, col="darkslateblue", lwd=1.5) # horizontal centerline at mu
  }

  # directions
  if (pause) cat("\n>>> Press Enter to obtain the next sample <<< \n\n")

  # get maximum length of each field 
  max.ln <- 0
  for (i in 1:ns) {
    nc <- nchar(as.character(round(Ymean[i],3)))
    if (nc > max.ln) max.ln <- nc
  }
  max.ln <- max.ln + 2

  # output
  if (n <= max.data) maxd <- n else maxd <- max.data
  cat("\nSample", "  Mean   ", "  SD   ")
  if (show.data) for (i in 1:maxd) cat(format(toString(i), width=max.ln,
                                       justify="right", sep=""))
  cat("\n")
  for (i in 1:ns) {
    if (pause) invisible(readline())
    cat(format(o[i], width=5, justify="right", sep=""))
    cat(format(sprintf("%.*f", digits.d, Ymean[i]), width=max.ln,
               justify="right", sep=""))
    cat(format(sprintf("%.*f", digits.d, Ysd[i]), width=max.ln,
               justify="right", sep=""))
    cat("   ")
    if (show.data) for (j in 1:maxd) 
      cat(format(sprintf("%.*f", digits.d, data.rep[i,j]), width=max.ln,
          justify="right", sep=""))
    if (n > max.data) cat(" ...")
    cat("\n")
    points(i, Ymean[i], pch=21, col="darkblue", bg="plum")
  }

  if (!pause) {
    if (show.title) title(main = bquote(paste(mu, "=", .(mu), "  ", sigma, "=", .(sigma), 
        "  ", "  n=", .(n))), cex.main=1)
    mtext(bquote(paste(" ", mu)), side=4, cex=1.5, col="darkslateblue", las=2)
    abline(h=mu, col="darkslateblue", lwd=1.5) # horizontal centerline at mu
  }

  if (pause) cat("\n")

  cat("\nAnalysis of Sample Means\n")
  cat("   Mean:", format(sprintf("%.*f", digits.d, mean(Ymean)), width=max.ln, justify="right", sep=""), "\n")
  cat("Std Dev:", format(sprintf("%.*f", digits.d, sd(Ymean)), width=max.ln, justify="right", sep=""))
  cat("    Direct estimate of the standard error of the sample mean\n")

  if (pause) {
    mu <<- mu
    sigma <<- sigma 
  }

  cat("\n")

  par(orig.params)

  # terminate pdf graphics system
  if (!is.null(pdf.file)) {
    dev.off()
    .showfile(pdf.file, "mean plot")
  }

}
