### This script plot X-Y protein production rates.

rm(list = ls())

suppressMessages(library(cubfits, quietly = TRUE))

source("00-set_env.r")
source(paste(prefix$code.plot, "u0-get_case_main.r", sep = ""))

### Pre processed phi.Obs.
fn.in <- paste(prefix$data, "pre_process.rda", sep = "")
load(fn.in)

### Keep _wphi_ cases only.
case.names <- case.names[grep("_wphi_", case.names)]

### Plot.
for(i.case in case.names){
  ### Subset of mcmc output.
  fn.in <- paste(prefix$subset, i.case, ".rda", sep = "")
  if(!file.exists(fn.in)){
    cat("File not found: ", fn.in, "\n", sep = "")
    next
  }
  load(fn.in)
  fn.in <- paste(prefix$subset, i.case, "_PM.rda", sep = "")
  if(!file.exists(fn.in)){
    cat("File not found: ", fn.in, "\n", sep = "")
    next
  }
  load(fn.in)

  ### Plot.
  fn.out <- paste(prefix$plot.diag, "sigmaW_", i.case, ".pdf", sep = "")
  pdf(fn.out, width = 5, height = 5)
    ret <- hist(p.mcmc[1,], nclass = 40,
                xlab = "std. err. of measurement errors", main = i.case)
    mtext(paste(workflow.name, ", ", get.case.main(i.case, model), sep = ""),
          line = 3, cex = 0.6)
    mtext(date(), line = 2.5, cex = 0.4)
    abline(v = p.PM[1], col = 2)
    text(p.PM[1] + 0.05 * diff(range(ret$breaks)),
         max(ret$counts) * 1.01, sprintf("%.4f", p.PM[1]),
         cex = 0.5)
  dev.off()
}

