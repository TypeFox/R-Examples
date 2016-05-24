rm(list = ls())

suppressMessages(library(cubfits, quietly = TRUE))

source("00-set_env.r")
source(paste(prefix$code.plot.ps, "u0-get_case_main.r", sep = ""))

for(i.case in case.names){
  ### Subset of mcmc output.
  fn.in <- paste(prefix$subset, i.case, ".rda", sep = "")
  if(!file.exists(fn.in)){
    cat("File not found: ", fn.in, "\n", sep = "")
    next
  }
  load(fn.in)

  ### Subset of mcmc output with scaling.
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

  phi.mcmc.1 <- rowMeans(phi.mcmc[, -1] != phi.mcmc[, -ncol(phi.mcmc)])

  ### Plot Phi_g.
  fn.out <- paste(prefix$plot.ps.diag, "acceptvsPhi_", i.case, ".pdf", sep = "")
  pdf(fn.out, width = 5, height = 5)
    plot(log10(phi.PM / mean(phi.PM)), phi.mcmc.1,
         xlab = "Predicted Expression (log10)", ylab = "acceptance",
         main = paste(i.case, " (posterior mean)", sep = ""),
         pch = 20, cex = 0.6)
    mtext(paste(workflow.name, ", ", get.case.main(i.case, model), sep = ""),
          line = 3, cex = 0.6)
    mtext(date(), line = 2.5, cex = 0.4)
  dev.off()
}

