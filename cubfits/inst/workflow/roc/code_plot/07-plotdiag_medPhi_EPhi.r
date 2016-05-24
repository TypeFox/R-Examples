### This script plot X-Y protein production rates.

rm(list = ls())

suppressMessages(library(cubfits, quietly = TRUE))

source("00-set_env.r")
source(paste(prefix$code.plot, "u0-get_case_main.r", sep = ""))

### Plot.
for(i.case in case.names){
  ### All mcmc outputs.
  fn.in <- paste(prefix$subset, i.case, ".rda", sep = "")
  if(!file.exists(fn.in)){
    cat(paste("File not found: ", fn.in, "\n", sep = ""))
    next
  }
  load(fn.in)

  medPhi <- apply(phi.mcmc, 2, median)
  meanPhi <- apply(phi.mcmc, 2, mean)

  ### Plot.
  fn.out <- paste(prefix$plot.diag, "medPhi_EPhi_", i.case, ".pdf", sep = "")
  pdf(fn.out, width = 5, height = 5)
    plotprxy(meanPhi, medPhi,
             log10.x = FALSE, log10.y = FALSE,
             xlab = "Mean of Phi", ylab = "Median of Phi",
             main = "By Iteration")
    mtext(paste(workflow.name, ", ", get.case.main(i.case, model), sep = ""),
          line = 3, cex = 0.6)
    mtext(date(), line = 2.5, cex = 0.4)
  dev.off()
}

