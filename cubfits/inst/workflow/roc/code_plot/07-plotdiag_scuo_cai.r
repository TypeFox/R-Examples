### This script plot X-Y protein production rates.

rm(list = ls())

suppressMessages(library(cubfits, quietly = TRUE))

source("00-set_env.r")

### Pre processed phi.Obs.
fn.in <- paste(prefix$data, "pre_process.rda", sep = "")
load(fn.in)

### Plot.
### x-axis: predicted, y-axis: observed.
fn.in <- file.data$tsv
if(file.exists(fn.in)){
  tmp <- log10(phi.Obs / mean(phi.Obs))
  fn.out <- paste(prefix$plot.diag, "prxy_scuo_cai.pdf", sep = "")
  pdf(fn.out, width = 12, height = 5)
    nf <- layout(matrix(c(1, 1, 1, 2, 3, 4),
                        nrow = 2, ncol = 3, byrow = TRUE),
                 c(1, 1, 1), c(2, 8), respect = FALSE)
    ### Plot title.
    par(mar = c(0, 0, 0, 0))
    plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
    text(0.5, 0.6, workflow.name)
    text(0.5, 0.4, date(), cex = 0.6)
    par(mar = c(5.1, 4.1, 4.1, 2.1))

    ### Plot SCUO vs phi.Obs.
    plotprxy(SCUO, tmp, log10.x = FALSE, log10.y = FALSE,
             xlab = "SCUO", ylab = "Observed Production Rate (log10)",
             main = "SCUO vs phi.Obs")

    ### Plot CAI vs phi.Obs.
    plotprxy(CAI$CAI, tmp, log10.x = FALSE, log10.y = FALSE,
             xlab = "CAI", ylab = "Observed Production Rate (log10)",
             main = "CAI vs phi.Obs")

    ### Plot SCUO vs CAI.
    plotprxy(SCUO, CAI$CAI, log10.x = FALSE, log10.y = FALSE,
             xlab = "SCUO", ylab = "CAI", main = "SCUO vs CAI")
  dev.off()
} else{
  fn.out <- paste(prefix$plot.diag, "prxy_scuo_cai.pdf", sep = "")
  pdf(fn.out, width = 4, height = 4)
    plotprxy(SCUO, CAI$CAI, log10.x = FALSE, log10.y = FALSE,
             xlab = "SCUO", ylab = "CAI", main = "SCUO vs CAI")
    mtext(workflow.name, line = 3, cex = 0.6)
    mtext(date(), line = 2.5, cex = 0.4)
  dev.off()
}

