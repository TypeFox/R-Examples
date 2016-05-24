### This script plot binning results and predicted curves from multinomial
### logistic regression without measurement errors.

rm(list = ls())

suppressMessages(library(cubfits, quietly = TRUE))

### Load environment and set data.
source("00-set_env.r")
fn.in <- paste(prefix$data, "pre_process.rda", sep = "")
load(fn.in)
fn.in <- paste(prefix$data, "init_", model, ".rda", sep = "")
load(fn.in)

### Plot bin and model for measurements.
fn.out <- paste(prefix$plot.diag, "prxy_init.pdf", sep = "")
pdf(fn.out, width = 8, height = 9)
  nf <- layout(matrix(c(1, 1, 2:5), nrow = 3, ncol = 2, byrow = TRUE),
               c(1, 1), c(2, 8, 8), respect = FALSE)
  ### Plot title.
  par(mar = c(0, 0, 0, 0))
  plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
  text(0.5, 0.6, workflow.name)
  text(0.5, 0.4, date(), cex = 0.6)
  par(mar = c(5.1, 4.1, 4.1, 2.1))

  ### Plot rests.
  plotprxy(phi.Obs, SCUO, xlab = "phi.Obs", ylab = "SCUO")
  plotprxy(phi.init.PM, phi.init.SCUO,
           xlab = "phi.init.PM", ylab = "phi.init.SCUO")
  hist(log(phi.Obs), nclass = 40, main = "log phi.Obs")
  hist(log(SCUO), nclass = 40, main = "log SCUO")
dev.off()

