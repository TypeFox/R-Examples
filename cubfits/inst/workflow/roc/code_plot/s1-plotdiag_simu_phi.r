### Plot some X-Y plots to make sure simulations are OK.

rm(list = ls())

suppressMessages(library(cubfits, quietly = TRUE))

### Set environment and read in data.
source("00-set_env.r")
fn.in <- paste(prefix$data, "simu_phi.tsv", sep = "")
if(file.exists(fn.in)){
  phi <- read.phi.df(fn.in)
} else{
  stop(paste(fn.in, " is not found.", sep = ""))
}

### Convert.
phi.Obs <- gen.phi.Obs(phi)
EPhi <- phi$true.phi

### Plot.
fn.out <- paste(prefix$plot.diag, "phiObs_vs_true.pdf", sep = "")
pdf(fn.out, width = 5, height = 5)
  ### x-axis: true, y-axis: observed.
  plotprxy(EPhi, phi.Obs, main = "True vs Observed",
           xlab = "Phi True", ylab = "Observed")
  mtext(workflow.name, line = 3, cex = 0.6)
  mtext(date(), line = 2.5, cex = 0.4)
dev.off()


