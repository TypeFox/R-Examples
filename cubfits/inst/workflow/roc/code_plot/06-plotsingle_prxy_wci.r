### This script plot X-Y protein production rates with labeled outliers.

rm(list = ls())

suppressMessages(library(cubfits, quietly = TRUE))

source("00-set_env.r")
source(paste(prefix$code.plot, "u0-get_case_main.r", sep = ""))

### Pre processed phi.Obs.
fn.in <- paste(prefix$data, "pre_process.rda", sep = "")
load(fn.in)

for(i.case in case.names){
  ### Subset of mcmc output.
  fn.in <- paste(prefix$subset, i.case, "_PM.rda", sep = "")
  if(!file.exists(fn.in)){
    cat("File not found: ", fn.in, "\n", sep = "")
    next
  }
  load(fn.in)
  # fn.in <- paste(prefix$subset, i.case, "_PM_scaling.rda", sep = "")
  # if(!file.exists(fn.in)){
  #   cat("File not found: ", fn.in, "\n", sep = "")
  #   next
  # }
  # load(fn.in)

  ### Plot posterior mean.
  fn.out <- paste(prefix$plot.single,
                  "prxy_wci_", i.case, "_nps.pdf", sep = "")
  pdf(fn.out, width = 5, height = 5)
    ### x-axis: predicted, y-axis: observed.
    plotprxy(phi.Obs, phi.PM,
             y.ci = phi.CI,
             xlab = "Observed Production Rate (log10)",
             ylab = "Predicted Production Rate (log10)",
             weights = 1 / phi.STD.log10,
             main = paste(i.case, " posterior mean", sep = ""))
    mtext(paste(workflow.name, ", ", get.case.main(i.case, model), sep = ""),
          line = 3, cex = 0.6)
    mtext(date(), line = 2.5, cex = 0.4)
  dev.off()

  ### Plot posterior median.
  fn.out <- paste(prefix$plot.single,
                  "prxy_wci_med_", i.case, "_nps.pdf", sep = "")
  pdf(fn.out, width = 5, height = 5)
    ### x-axis: predicted, y-axis: observed.
    plotprxy(phi.Obs, phi.MED,
             y.ci = phi.CI,
             xlab = "Observed Production Rate (log10)",
             ylab = "Predicted Production Rate (log10)",
             weights = 1 / phi.STD.log10,
             main = paste(i.case, " posterior median", sep = ""))
    mtext(paste(workflow.name, ", ", get.case.main(i.case, model), sep = ""),
          line = 3, cex = 0.6)
    mtext(date(), line = 2.5, cex = 0.4)
  dev.off()

  ### Plot posterior log10 mean.
  fn.out <- paste(prefix$plot.single,
                  "prxy_wci_log10_", i.case, "_nps.pdf", sep = "")
  pdf(fn.out, width = 5, height = 5)
    ### x-axis: predicted, y-axis: observed.
    plotprxy(phi.Obs, 10^(phi.PM.log10),
             xlab = "Observed Production Rate (log10)",
             ylab = "Predicted Production Rate (log10)",
             y.ci = 10^(phi.CI.log10),
             weights = 1 / phi.STD.log10,
             main = paste(i.case, " posterior log10 mean", sep = ""))
    mtext(paste(workflow.name, ", ", get.case.main(i.case, model), sep = ""),
          line = 3, cex = 0.6)
    mtext(date(), line = 2.5, cex = 0.4)
  dev.off()
}

