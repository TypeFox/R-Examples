### This script plots comparisions across cases.

rm(list = ls())

suppressMessages(library(cubfits, quietly = TRUE))

source("00-set_env.r")

if(length(case.names) < 4){
  stop("Need 4 cases to match with.")
}

### Ordered by "wophi_pm", "wophi_scuo", "wphi_pm", and "wphi_scuo".
phi.mean <- list(NULL, NULL, NULL, NULL)
phi.median <- list(NULL, NULL, NULL, NULL)
phi.std <- list(NULL, NULL, NULL, NULL)
phi.ci <- list(NULL, NULL, NULL, NULL)
for(i.case in 1:4){
  ### Subset of mcmc output.
  fn.in <- paste(prefix$subset, case.names[i.case], "_PM.rda", sep = "")
  if(!file.exists(fn.in)){
    cat("File not found: ", fn.in, "\n", sep = "")
    next
  }
  load(fn.in)

  ### Subset of mcmc output with scaling.
  fn.in <- paste(prefix$subset, case.names[i.case], "_PM_scaling.rda", sep = "")
  if(!file.exists(fn.in)){
    cat("File not found: ", fn.in, "\n", sep = "")
    next
  }
  load(fn.in)

  phi.mean[[i.case]] <- phi.PM
  phi.median[[i.case]] <- phi.MED
  phi.std[[i.case]] <- phi.STD.log10
  phi.ci[[i.case]] <- phi.CI
}

### Plot posterior mean.
if(!is.null(phi.mean[[3]]) && !is.null(phi.mean[[1]])){
  fn.out <- paste(prefix$plot.ps.match, "prxy_wci_pm.pdf", sep = "")
  pdf(fn.out, width = 5, height = 5)
    ### x-axis: with phi, y-axis: without phi.
    plotprxy(phi.mean[[3]], phi.mean[[1]], weights = 1 / phi.std[[1]],
             y.ci = phi.ci[[1]],
             xlab = "Production Rate with phi (log10)",
             ylab = "Production Rate without phi (log10)",
             main = "fits vs appr (pm, posterior mean)")
    mtext(workflow.name, line = 3, cex = 0.6)
    mtext(date(), line = 2.5, cex = 0.4)
  dev.off()
}

if(!is.null(phi.mean[[4]]) && !is.null(phi.mean[[2]])){
  fn.out <- paste(prefix$plot.ps.match, "prxy_wci_scuo.pdf", sep = "")
  pdf(fn.out, width = 5, height = 5)
    ### x-axis: with phi, y-axis: without phi.
    plotprxy(phi.mean[[4]], phi.mean[[2]], weights = 1 / phi.std[[2]],
             y.ci = phi.ci[[2]],
             xlab = "Production Rate with phi (log10)", 
             ylab = "Production Rate without phi (log10)",
             main = "fits vs appr (scuo, posterior mean)")
    mtext(workflow.name, line = 3, cex = 0.6)
    mtext(date(), line = 2.5, cex = 0.4)
  dev.off()
}

### Plot posterior median.
if(!is.null(phi.mean[[3]]) && !is.null(phi.mean[[1]])){
  fn.out <- paste(prefix$plot.ps.match, "prxy_wci_med_pm.pdf", sep = "")
  pdf(fn.out, width = 5, height = 5)
    ### x-axis: with phi, y-axis: without phi.
    plotprxy(phi.median[[3]], phi.median[[1]], weights = 1 / phi.std[[1]],
             y.ci = phi.ci[[1]],
             xlab = "Production Rate with phi (log10)",
             ylab = "Production Rate without phi (log10)",
             main = "fits vs appr (pm, posterior median)")
    mtext(workflow.name, line = 3, cex = 0.6)
    mtext(date(), line = 2.5, cex = 0.4)
  dev.off()
}

if(!is.null(phi.mean[[4]]) && !is.null(phi.mean[[2]])){
  fn.out <- paste(prefix$plot.ps.match, "prxy_wci_med_scuo.pdf", sep = "")
  pdf(fn.out, width = 5, height = 5)
    ### x-axis: with phi, y-axis: without phi.
    plotprxy(phi.median[[4]], phi.median[[2]], weights = 1 / phi.std[[2]],
             y.ci = phi.ci[[2]],
             xlab = "Production Rate with phi (log10)", 
             ylab = "Production Rate without phi (log10)",
             main = "fits vs appr (scuo, posterior median)")
    mtext(workflow.name, line = 3, cex = 0.6)
    mtext(date(), line = 2.5, cex = 0.4)
  dev.off()
}

