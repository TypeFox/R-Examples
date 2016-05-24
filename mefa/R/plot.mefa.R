`plot.mefa` <-
function(x, stat=1:4, type=c("hist", "rank"), trafo=c("none", "log",
"ratio"), show=TRUE, ylab, xlab, ...)
{
    if (missing(ylab))
        ylab <- NULL
    if (missing(xlab))
        xlab <- NULL
   if (!all(stat %in% 1:4))
       stop("'stat' must be in 1:4")
   if (!length(stat) == 1) stat <- 1
   if (!length(type) == 1) type <- type[1]
   if (!length(trafo) == 1) trafo <- trafo[1]
   trafo <- match.arg(trafo, c("none", "log", "ratio"))
   ## maintaining back compatibility
   type <- match.arg(type, c("hist", "rank", "bar"))
   if (type == "hist")
       type <- "bar"

   if (!is.null(ylab))
       ylab2 <- ylab
   if (!is.null(xlab))
       xlab2 <- xlab
   if (is.null(ylab) && type=="bar")
       ylab2 <- "Frequency"
   if (is.null(ylab) && type=="rank")
       ylab2 <- "Rank"

   if (stat == 1) {
       if (is.null(ylab)) ylab2 <- paste(ylab2, "(samples)")
       if (is.null(xlab)) xlab2 <- "Number of taxa"
       }
   if (stat == 2) {
       if (is.null(ylab)) ylab2 <- paste(ylab2, "(samples)")
       if (is.null(xlab)) xlab2 <- "Number of individuals"
       }
   if (stat == 3) {
       if (is.null(ylab)) ylab2 <- paste(ylab2, "(taxa)")
       if (is.null(xlab)) xlab2 <- "Frequency of occurrence"
       }
   if (stat == 4) {
       if (is.null(ylab)) ylab2 <- paste(ylab2, "(taxa)")
       if (is.null(xlab)) xlab2 <- "Abundance"
       }
    yvar <- summary(x)[[stat]]

   if (trafo=="log") {
       yvar <- log10(yvar)
       if (is.null(ylab) && type=="bar")
           ylab2 <- paste("log10", ylab2)
       if (is.null(ylab) && type=="rank")
           xlab2 <- paste("log10", xlab2)}
   if (trafo=="ratio") {
       yvar <- yvar / max(yvar)
       if (is.null(ylab) && type=="bar")
           ylab2 <- paste("Relative", tolower(ylab2))
       if (is.null(ylab) && type=="rank")
           xlab2 <- paste("Relative", tolower(xlab2))}


       if (type=="bar") {
           if (all(as.integer(yvar) == as.numeric(yvar))) {
               yvar <- table(yvar)
               if (show)
                    plot(yvar, xlab=xlab2, ylab=ylab2, ...)
               } else {
                    hist(yvar, xlab=xlab2, ylab=ylab2, ...)
               }
       }
       if (type=="rank") {
           yvar <- yvar[order(yvar, decreasing=TRUE)]
           if (show)
                plot(yvar, type="l", xlab=ylab2, ylab=xlab2, ...)}
    if (show)
        invisible(yvar) else return(yvar)
}
