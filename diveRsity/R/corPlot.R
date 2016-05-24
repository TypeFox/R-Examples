# corPlot, plot the relationship between fastDivPart stats and number of 
# alleles
#' @export
corPlot <-function(infile = NULL, write = FALSE, plot.format = NULL){
  # make sure mfrow is reset on exit
  orig <- par()$mfrow
  on.exit(par(mfrow = orig))
  # calculate the number of alleles per locus
  na <- sapply(diveRsity::rgp(infile)$af, function(x){
    mean(colSums(x > 0), na.rm = TRUE)
  })
  # generate the differentiation stats
  colselect <- c("Fst", "gst", "Gst", "D")
  dat <- diveRsity::diffCalc(infile, fst = TRUE)$std_stats[,colselect]
  dat <- dat[-nrow(dat),]
  colselect <- c("Fst", "gst", "Gprimest", "D")
  colnames(dat) <- colselect
  dat$na <- na
  ynms <- c(expression("F"[ST]), expression("G"[ST]), expression("G'"[ST]), 
            expression("D"[JOST]))
  if(write){
    plts <- mapply(corPlotter, y = colselect, yname = ynms, 
                   MoreArgs = list(x = "na", dat = dat,
                                   write = TRUE, plot.format = plot.format), 
                   SIMPLIFY = FALSE)
    multiplot(plts[[1]], plts[[2]], plts[[3]], plts[[4]], cols = 2)
  } else {
    plts <- mapply(corPlotter, y = colselect, yname = ynms, 
                   MoreArgs = list(x = "na", dat = dat), SIMPLIFY = FALSE)
    multiplot(plts[[1]], plts[[2]], plts[[3]], plts[[4]], cols = 2)
  }
}