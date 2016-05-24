rm(list = ls())

suppressMessages(library(cubfits, quietly = TRUE))

source("00-set_env.r")
source(paste(prefix$code.plot, "u0-get_case_main.r", sep = ""))

for(i.case in case.names){
  ### Subset of mcmc output.
  fn.in <- paste(prefix$subset, i.case, ".rda", sep = "")
  if(!file.exists(fn.in)){
    cat("File not found: ", fn.in, "\n", sep = "")
    next
  }
  load(fn.in)

  b.mcmc.1 <- rowMeans(b.mcmc[, -1] != b.mcmc[, -ncol(b.mcmc)])
  phi.mcmc.1 <- rowMeans(phi.mcmc[, -1] != phi.mcmc[, -ncol(phi.mcmc)])

  ### Plot Delta.t/M.
  fn.out <- paste(prefix$plot.diag, "accept_logmu_deltat_", i.case, ".pdf", sep = "")
  pdf(fn.out, width = 5, height = 5)
    hist(b.mcmc.1, main = paste(i.case, "_logmu_deltat", sep = ""),
         xlim = c(0, 1), nclass = 25)
    mtext(paste(workflow.name, ", ", get.case.main(i.case, model), sep = ""),
          line = 3, cex = 0.6)
    mtext(date(), line = 2.5, cex = 0.4)
  dev.off()

  ### Plot Phi_g.
  fn.out <- paste(prefix$plot.diag, "accept_Phi_", i.case, ".pdf", sep = "")
  pdf(fn.out, width = 5, height = 5)
    hist(phi.mcmc.1, main = paste(i.case, "_Phi", sep = ""),
         xlim = c(0, 1), nclass = 25)
    mtext(paste(workflow.name, ", ", get.case.main(i.case, model), sep = ""),
          line = 3, cex = 0.6)
    mtext(date(), line = 2.5, cex = 0.4)
  dev.off()
}

