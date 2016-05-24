### This script plot X-Y protein production rates.

rm(list = ls())

suppressMessages(library(cubfits, quietly = TRUE))

source("00-set_env.r")
source(paste(prefix$code.plot.ps, "u0-get_case_main.r", sep = ""))

### Plot.
for(i.case in case.names){
  ### Subset of mcmc output.
  fn.in <- paste(prefix$subset, i.case, "_PM.rda", sep = "")
  if(!file.exists(fn.in)){
    cat("File not found: ", fn.in, "\n", sep = "")
    next
  }
  load(fn.in)

  ### Subset of mcmc output with scaling.
  fn.in <- paste(prefix$subset, i.case, "_PM_scaling.rda", sep = "")
  if(!file.exists(fn.in)){
    cat("File not found: ", fn.in, "\n", sep = "")
    next
  }
  load(fn.in)

  ### Plot.
  fn.out <- paste(prefix$plot.ps.diag, "EPhi_", i.case, ".pdf", sep = "")
  pdf(fn.out, width = 5, height = 5)
    ret <- hist(phi.PM, nclass = 50, xlab = "EPhi", main = i.case)
    mtext(paste(workflow.name, ", ", get.case.main(i.case, model), sep = ""),
          line = 3, cex = 0.6)
    mtext(date(), line = 2.5, cex = 0.4)
    abline(v = mean(phi.PM), col = 2)
    text(mean(phi.PM) + 0.05 * diff(range(ret$breaks)),
         max(ret$counts) * 1.01, sprintf("%.4f", mean(phi.PM)),
         cex = 0.5)
  dev.off()
}

