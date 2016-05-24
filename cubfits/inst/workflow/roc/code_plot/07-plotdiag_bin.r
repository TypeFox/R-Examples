rm(list = ls())

suppressMessages(library(cubfits, quietly = TRUE))

### Load environment and set data.
source("00-set_env.r")
source(paste(prefix$code.plot, "u0-get_case_main.r", sep = ""))
fn.in <- paste(prefix$data, "pre_process.rda", sep = "")
load(fn.in)

### Arrange data.
phi.Obs.lim <- range(phi.Obs)
aa.names <- names(reu13.df.obs)
ret.phi.Obs <- prop.bin.roc(reu13.df.obs, phi.Obs,
                            bin.class = run.info$bin.class)

for(i.case in case.names){
  ### Subset of mcmc output.
  fn.in <- paste(prefix$subset, i.case, "_PM.rda", sep = "")
  if(!file.exists(fn.in)){
    cat("File not found: ", fn.in, "\n", sep = "")
    next
  }
  load(fn.in)

  ### Subset of mcmc output with scaling.
  # fn.in <- paste(prefix$subset, i.case, "_PM_scaling.rda", sep = "")
  # if(!file.exists(fn.in)){
  #   cat("File not found: ", fn.in, "\n", sep = "")
  #   next
  # }
  # load(fn.in)

  ### To adjust to similar range of phi.Obs.
  ret.EPhi <- prop.bin.roc(reu13.df.obs, phi.PM,
                           bin.class = run.info$bin.class)
  b.PM <- convert.bVec.to.b(b.PM, aa.names)
  EPhi.lim <- range(c(phi.Obs.lim, phi.PM))
  predict.roc <- prop.model.roc(b.PM, EPhi.lim)

### For phiObs.
  ### Fix xlim at log10 scale.
  lim.bin <- range(log10(ret.phi.Obs[[1]]$center))
  xlim <- c(lim.bin[1] - diff(lim.bin) / 4, lim.bin[2] + diff(lim.bin) / 4)

  ### Plot bin and model for measurements.
  fn.out <- paste(prefix$plot.diag, "bin_pred_phiObs_", i.case, ".pdf", sep = "")
  pdf(fn.out, width = 16, height = 11)
    mat <- matrix(c(rep(1, 5), 2:21, rep(22, 5)),
                  nrow = 6, ncol = 5, byrow = TRUE)
    mat <- cbind(rep(23, 6), mat, rep(24, 6))
    nf <- layout(mat, c(3, rep(8, 5), 2), c(3, 8, 8, 8, 8, 3), respect = FALSE)
    ### Plot title.
    par(mar = c(0, 0, 0, 0))
    plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
    text(0.5, 0.6,
         paste(workflow.name, ", ", get.case.main(i.case, model),
               ", bin: observed phi", sep = ""))
    text(0.5, 0.4, date(), cex = 0.6)

    ### Plot results.
    for(i.aa in 1:length(aa.names)){
      tmp.obs <- ret.phi.Obs[[i.aa]]
      tmp.roc <- predict.roc[[i.aa]]
      plotbin(tmp.obs, tmp.roc, main = "", xlab = "", ylab = "",
              lty = 1, axes = FALSE, xlim = xlim)
      box()
      text(mean(xlim), 1, aa.names[i.aa], cex = 1.5)
      if(i.aa %in% c(1, 6, 11, 16)){
        axis(2)
      }
      if(i.aa %in% 16:19){
        axis(1)
      }
      if(i.aa %in% 1:5){
        axis(3)
      }
      if(i.aa %in% c(5, 10, 15)){
        axis(4)
      }
      axis(1, tck = 0.02, labels = FALSE)
      axis(2, tck = 0.02, labels = FALSE)
      axis(3, tck = 0.02, labels = FALSE)
      axis(4, tck = 0.02, labels = FALSE)
    }

    ### For cases with less aa.
    i.aa <- 19 - i.aa
    if(i.aa > 0){
      for(i.plot in 1:i.aa){
        plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1),
             xlab = "", ylab = "", main = "", axes = FALSE)
      }
    }

    ### Add histogram.
    p.1 <- hist(log10(phi.Obs), nclass = 40, plot = FALSE)
    hist.ylim <- range(p.1$counts)
    hist.ylim[2] <- hist.ylim[2] + 0.2 * diff(hist.ylim)
    plot(p.1, xlim = xlim, ylim = hist.ylim, main = "", xlab = "", ylab = "",
         axes = FALSE)
    axis(1)
    axis(4)

    ### Add label.
    model.label <- c("MCMC Posterior")
    model.lty <- 1
    legend(xlim[1], hist.ylim[2],
           model.label, lty = model.lty, box.lty = 0, cex = 0.8)

    ### Plot xlab.
    plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
    text(0.5, 0.5,
         expression(paste(log[10], "(Observed Production Rate)", sep = "")))

    ### Plot ylab.
    plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
    text(0.5, 0.5, "Codon Frequency", srt = 90)
  dev.off()

### For EPhi.
  ### Fix xlim at log10 scale.
  lim.bin <- range(log10(ret.EPhi[[1]]$center))
  xlim <- c(lim.bin[1] - diff(lim.bin) / 4, lim.bin[2] + diff(lim.bin) / 4)

  ### Plot bin and model for predictions.
  fn.out <- paste(prefix$plot.diag, "bin_pred_EPhi_", i.case, ".pdf", sep = "")
  pdf(fn.out, width = 16, height = 11)
    mat <- matrix(c(rep(1, 5), 2:21, rep(22, 5)),
                  nrow = 6, ncol = 5, byrow = TRUE)
    mat <- cbind(rep(23, 6), mat, rep(24, 6))
    nf <- layout(mat, c(3, rep(8, 5), 2), c(3, 8, 8, 8, 8, 3), respect = FALSE)
    ### Plot title.
    par(mar = c(0, 0, 0, 0))
    plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
    text(0.5, 0.6,
         paste(workflow.name, ", ", get.case.main(i.case, model),
               ", bin: posterior mean of Phi", sep = ""))
    text(0.5, 0.4, date(), cex = 0.6)

    ### Plot results.
    for(i.aa in 1:length(aa.names)){
      tmp.obs <- ret.EPhi[[i.aa]]
      tmp.roc <- predict.roc[[i.aa]]
      plotbin(tmp.obs, tmp.roc, main = "", xlab = "", ylab = "",
              lty = 1, axes = FALSE, xlim = xlim)
      box()
      text(mean(xlim), 1, aa.names[i.aa], cex = 1.5)
      if(i.aa %in% c(1, 6, 11, 16)){
        axis(2)
      }
      if(i.aa %in% 16:19){
        axis(1)
      }
      if(i.aa %in% 1:5){
        axis(3)
      }
      if(i.aa %in% c(5, 10, 15)){
        axis(4)
      }
      axis(1, tck = 0.02, labels = FALSE)
      axis(2, tck = 0.02, labels = FALSE)
      axis(3, tck = 0.02, labels = FALSE)
      axis(4, tck = 0.02, labels = FALSE)
    }

    ### For cases with less aa.
    i.aa <- 19 - i.aa
    if(i.aa > 0){
      for(i.plot in 1:i.aa){
        plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1),
             xlab = "", ylab = "", main = "", axes = FALSE)
      }
    }

    ### Add histogram.
    p.1 <- hist(log10(phi.PM), nclass = 40, plot = FALSE)
    hist.ylim <- range(p.1$counts)
    hist.ylim[2] <- hist.ylim[2] + 0.2 * diff(hist.ylim)
    plot(p.1, xlim = xlim, ylim = hist.ylim, main = "", xlab = "", ylab = "",
         axes = FALSE)
    axis(1)
    axis(4)

    ### Add label.
    model.label <- c("MCMC Posterior")
    model.lty <- 1
    legend(xlim[1], hist.ylim[2],
           model.label, lty = model.lty, box.lty = 0, cex = 0.8)

    ### Plot xlab.
    plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
    text(0.5, 0.5,
         expression(paste(log[10], "(Estimated Production Rate)", sep = "")))

    ### Plot ylab.
    plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
    text(0.5, 0.5, "Codon Frequency", srt = 90)
  dev.off()
}
