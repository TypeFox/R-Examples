### This is for simulation only, ploting correlation aganst true values.

rm(list = ls())

suppressMessages(library(cubfits, quietly = TRUE))

source("00-set_env.r")
source(paste(prefix$code.plot, "u0-get_case_main.r", sep = ""))
source(paste(prefix$code, "u1-get_negsel.r", sep = ""))
source(paste(prefix$code.plot, "u2-plot_b_corr.r", sep = ""))
source(paste(prefix$code.plot, "u5-new_page.r", sep = ""))
source(paste(prefix$code.plot, "u6-adjust_focal_codon.r", sep = ""))

### Load true Phi.
fn.in <- paste(prefix$data, "simu_true_", model, ".rda", sep = "")
if(file.exists(fn.in)){
  load(fn.in)
} else{
  stop(paste(fn.in, " is not found.", sep = ""))
}
b.Init <- convert.b.to.bVec(Eb)


### Load results from collapsed runs for yassour with phi since
### wphi_wophi does not need to simulate new sequences.
# fn.in <- paste(prefix$param, "small_bInit.rda", sep = "")
# load(fn.in)
# b.Init <- convert.b.to.bVec(b.InitList.roc)
# fn.in <- paste(prefix$param, "small_train.rda", sep = "")
# load(fn.in)
# EPhi <- phi.Obs


### Load data.
fn.in <- paste(prefix$data, "pre_process.rda", sep = "")
load(fn.in)

### Get AA and synonymous codons.
aa.names <- names(reu13.df.obs)
coef.names <- cubfits:::get.my.coefnames(model)
label <- NULL
b.names <- NULL
for(i.aa in aa.names){
  tmp <- sort(unique(reu13.df.obs[[i.aa]]$Codon))
  tmp <- tmp[-length(tmp)]
  label <- c(label, paste(i.aa, tmp, sep = "."))
  b.names <- c(b.names, rep(coef.names, each = length(tmp)))
}

### Get true values.
all.names <- b.names 
id.intercept <- grep("log.mu", all.names)
id.slop <- grep("Delta.t", all.names)

scale.EPhi <- mean(EPhi)
b.Init[id.slop] <- b.Init[id.slop] * scale.EPhi
b.Init.negsel <- get.negsel(b.Init, id.intercept, id.slop, aa.names, label)

### True SCU.
EPhi <- EPhi / scale.EPhi
Eb <- convert.bVec.to.b(b.Init, aa.names)
SCU.true <- calc_scu_values(Eb, y.list, EPhi)


### Plot by case.
for(i.case in case.names){
  ### Check files.
  fn.in <- paste(prefix$subset, i.case, ".rda", sep = "")
  if(!file.exists(fn.in)){
    cat("File not found: ", fn.in, "\n", sep = "")
    next
  }
  fn.in <- paste(prefix$subset, i.case, "_PM.rda", sep = "")
  if(!file.exists(fn.in)){
    cat("File not found: ", fn.in, "\n", sep = "")
    next
  }
  fn.in <- paste(prefix$output, i.case, "/output_mcmc.rda", sep = "")
  if(!file.exists(fn.in)){
    cat("File not found: ", fn.in, "\n", sep = "")
    next
  }


  ### Subset of mcmc output.
  fn.in <- paste(prefix$subset, i.case, ".rda", sep = "")
  load(fn.in)

  ### Subset of mcmc output.
  fn.in <- paste(prefix$subset, i.case, "_PM.rda", sep = "")
  load(fn.in)

  ### Subset of mcmc output with scaling.
  # fn.in <- paste(prefix$subset, i.case, "_PM_scaling.rda", sep = "")
  # load(fn.in)

  ### All mcmc outputs.
  fn.in <- paste(prefix$output, i.case, "/output_mcmc.rda", sep = "")
  load(fn.in)

  ### Since my.appr() doesn't have phi.Mat, but have phi.pred.Mat
  if(is.null(ret[["phi.Mat"]])){
    ret$phi.Mat <- ret$phi.pred.Mat
  }


  ### Set layout.
  fn.out <- paste(prefix$plot.multi, "true_", i.case, "_nps.pdf", sep = "")
  pdf(fn.out, width = 6, height = 10)

### New page.
    new.page(workflow.name, i.case, model)

    ### Plot Delta.t.
    x <- b.Init.negsel$b.negsel.PM
    y <- b.negsel.PM
    y.ci <- b.negsel.ci.PM
    x.label <- b.Init.negsel$b.negsel.label
    plot.b.corr(x, y, x.label, y.ci = y.ci,
                # xlab = "True", ylab = "Estimated", main = "Delta.t",
                add.lm = TRUE, add.ci = TRUE)
    mtext(expression(paste(italic(Delta[t]), " , true")),
          side = 1, line = 2.5, cex = 0.8)
    if(length(grep("wophi", i.case)) > 0){
      mtext(expression(paste(italic(Delta[t]), " , without ", italic(X[obs]))),
            side = 2, line = 2.5, cex = 0.8)
    } else{
      mtext(expression(paste(italic(Delta[t]), " , with ", italic(X[obs]))),
            side = 2, line = 2.5, cex = 0.8)
    }
    x.label.focal <- x.label

    ### Plot log(mu).
    x <- b.Init.negsel$b.logmu.PM
    y <- b.logmu.PM
    y.ci <- b.logmu.ci.PM
    x.label <- b.logmu.label
    ### Adjust focal codons if needed and dispatch after adjusting.
    new.order <- adjust.focal.codon(y, x.label, x.label.focal, y.ci = y.ci)
    y <- new.order$y
    y.ci <- new.order$y.ci
    plot.b.corr(x, y, x.label.focal, y.ci = y.ci,
                # xlab = "True", ylab = "Estimated", main = "log(mu)",
                add.lm = TRUE, add.ci = TRUE)
    mtext(expression(paste(italic(M), " , true")),
          side = 1, line = 2.5, cex = 0.8)
    if(length(grep("wophi", i.case)) > 0){
      mtext(expression(paste(italic(M), " , without ", italic(X[obs]))),
            side = 2, line = 2.5, cex = 0.8)
    } else{
      mtext(expression(paste(italic(M), " , with ", italic(X[obs]))),
            side = 2, line = 2.5, cex = 0.8)
    }

    ### Overlap two histograms.
    p.1 <- hist(log10(EPhi / mean(EPhi)), nclass = 50, plot = FALSE)
    p.2 <- hist(log10(phi.PM / mean(phi.PM)), nclass = 50, plot = FALSE)
    xlim <- range(p.1$breaks, p.2$breaks)
    ylim <- range(p.1$counts, p.2$counts)
    plot(p.1, col = "#0000FF50", xlim = xlim, ylim = ylim,
         xlab = "Production Rate (log10)",
         main = "Expression (Posterior Mean)")
    plot(p.2, col = "#FF000050", xlim = xlim, ylim = ylim, add = TRUE)
    legend(xlim[1], ylim[2], c("True", "Predicted"),
           pch = c(15, 15), col = c("#0000FF50", "#FF000050"), cex = 0.8)

    ### Overlap two histograms.
    p.1 <- hist(log10(EPhi / mean(EPhi)), nclass = 50, plot = FALSE)
    p.2 <- hist(log10(phi.MED / mean(phi.MED)), nclass = 50, plot = FALSE)
    xlim <- range(p.1$breaks, p.2$breaks)
    ylim <- range(p.1$counts, p.2$counts)
    plot(p.1, col = "#0000FF50", xlim = xlim, ylim = ylim,
         xlab = "Production Rate (log10)",
         main = "Expression (Posterior Median)")
    plot(p.2, col = "#FF000050", xlim = xlim, ylim = ylim, add = TRUE)
    legend(xlim[1], ylim[2], c("True", "Predicted"),
           pch = c(15, 15), col = c("#0000FF50", "#FF000050"), cex = 0.8)

    ### Plot SCU.
    b <- convert.bVec.to.b(b.PM, aa.names)
    SCU <- calc_scu_values(b, y.list, phi.PM)

    plotprxy(SCU.true$SCU, SCU$SCU,
             xlab = "True SCU (log10)",
             ylab = "Predicted SCU (log10)",
             main = "SCU (Posterior Mean)")

    ### Plot mSCU.
    plotprxy(SCU.true$mSCU, SCU$mSCU,
             log10.x = FALSE, log10.y = FALSE,
             xlab = "True mSCU",
             ylab = "Predicted mSCU",
             main = "mSCU (Posterior Mean)")


### New page.
    new.page(workflow.name, i.case, model)

    ### For x-y plots.
    x <- log10(EPhi / mean(EPhi))
    y <- log10(phi.PM / mean(phi.PM))
    xlim <- ylim <- range(c(x, y))

    ### Plot prxy.
    plotprxy(EPhi, 10^(phi.PM.log10),
             weights = 1 / phi.STD.log10,
             xlim = xlim, ylim = ylim,
             xlab = "True Production Rate (log10)",
             ylab = "Predicted Production Rate (log10)",
             main = "Expression (Posterior log10 Mean)")

    ### Plot prxy with outliers labeled.
    plotprxy(EPhi, 10^(phi.PM.log10),
             y.ci = 10^(phi.CI.log10),
             weights = 1 / phi.STD.log10,
             xlim = xlim, ylim = ylim,
             xlab = "True Production Rate (log10)",
             ylab = "Predicted Production Rate (log10)",
             main = "Expression (Posterior log10 Mean)")

    ### Plot prxy.
    plotprxy(EPhi, phi.MED,
             weights = 1 / phi.STD.log10,
             xlim = xlim, ylim = ylim,
             xlab = "True Production Rate (log10)",
             ylab = "Predicted Production Rate (log10)",
             main = "Expression (Posterior Median)")

    ### Plot prxy with outliers labeled.
    plotprxy(EPhi, phi.MED,
             y.ci = phi.CI,
             weights = 1 / phi.STD.log10,
             xlim = xlim, ylim = ylim,
             xlab = "True Production Rate (log10)",
             ylab = "Predicted Production Rate (log10)",
             main = "Expression (Posterior Median)")

    ### Plot prxy.
    plotprxy(EPhi, phi.PM,
             weights = 1 / phi.STD.log10,
             xlim = xlim, ylim = ylim,
             xlab = "True Production Rate (log10)",
             ylab = "Predicted Production Rate (log10)",
             main = "Expression (Posterior Mean)")

    ### Plot prxy with outliers labeled.
    plotprxy(EPhi, phi.PM,
             y.ci = phi.CI,
             weights = 1 / phi.STD.log10,
             xlim = xlim, ylim = ylim,
             xlab = "True Production Rate (log10)",
             ylab = "Predicted Production Rate (log10)",
             main = "Expression (Posterior Mean)")


### New page.
    new.page(workflow.name, i.case, model)

    ### Add qqplot.
    x <- log10(EPhi / mean(EPhi))
    y <- log10(phi.PM / mean(phi.PM))
    xlim <- ylim <- range(c(x, y))
    qqplot(x, y,
           xlim = xlim, ylim = ylim,
           xlab = "True Production Rate (log10)",
           ylab = "Predicted Production Rate (log10)",
           main = "Q-Q (Posterior Mean)",
           pch = 20, cex = 0.6)
    abline(a = 0, b = 1, col = 4, lty = 2)

    ### Add qqplot for non-log.
    x <- EPhi / mean(EPhi)
    y <- phi.PM / mean(phi.PM)
    xlim <- ylim <- range(c(x, y))
    qqplot(x, y,
           xlim = xlim, ylim = ylim,
           xlab = "True Production Rate",
           ylab = "Predicted Production Rate",
           main = "Q-Q (Posterior Mean)",
           pch = 20, cex = 0.6)
    abline(a = 0, b = 1, col = 4, lty = 2)

    ### Add qqplot.
    x <- log10(EPhi / mean(EPhi))
    y <- log10(phi.MED / mean(phi.MED))
    xlim <- ylim <- range(c(x, y))
    qqplot(x, y,
           xlim = xlim, ylim = ylim,
           xlab = "True Production Rate (log10)",
           ylab = "Predicted Production Rate (log10)",
           main = "Q-Q (Posterior Median)",
           pch = 20, cex = 0.6)
    abline(a = 0, b = 1, col = 4, lty = 2)

    ### Add qqplot for non-log.
    x <- EPhi / mean(EPhi)
    y <- phi.MED / mean(phi.MED)
    xlim <- ylim <- range(c(x, y))
    qqplot(x, y,
           xlim = xlim, ylim = ylim,
           xlab = "True Production Rate",
           ylab = "Predicted Production Rate",
           main = "Q-Q (Posterior Median)",
           pch = 20, cex = 0.6)
    abline(a = 0, b = 1, col = 4, lty = 2)


### New page.
    new.page(workflow.name, i.case, model)

    ### Plot mediadn Phi vs mean Phi.
    medPhi <- apply(phi.mcmc, 2, median)
    meanPhi <- apply(phi.mcmc, 2, mean)
    plotprxy(meanPhi, medPhi,
             log10.x = FALSE, log10.y = FALSE,
             xlab = "Mean of Phi", ylab = "Median of Phi",
             main = "By Iteration")

    ### Plot trace of EPhi.
    trace <- lapply(1:length(ret$phi.Mat), function(i){ mean(ret$phi.Mat[[i]]) })
    trace <- do.call("c", trace)

    xlim <- c(1, length(ret$phi.Mat))
    ylim <- range(c(range(trace), 1))
    plot(NULL, NULL, xlim = xlim, ylim = ylim,
         xlab = "Iterations", ylab = "Mean of EPhi",
         main = "Trace of mean expression")
    lines(x = 1:length(ret$phi.Mat), y = trace)
    abline(h = 1, col = 2)

    ### Get pacf.
    b.pacf <- NULL
    for(i in 1:nrow(b.mcmc)){
      tmp <- which(abs(pacf(b.mcmc[i,], plot = FALSE)$acf[, 1, 1]) < 0.1)
      if(length(tmp) == 0){
        tmp <- NA
      }
      b.pacf <- c(b.pacf, tmp[1])
    }

    ### Plot PACF Delta.t.
    x <- b.pacf[id.slop]
    hist(x, nclass = 10,
         xlab = "Min Lag Partial Corr. < 0.1",
         main = paste("PACF Delta.t\nmean=",
                      sprintf("%4.2f", mean(x, na.rm = TRUE)),
                      " med.=",
                      sprintf("%4.2f", median(x, na.rm = TRUE)), sep = ""))

    ### Plot PACF log(mu).
    x <- b.pacf[id.intercept]
    hist(x, nclass = 10,
         xlab = "Min Lag Partial Corr. < 0.1",
         main = paste("PACF log(mu)\nmean=",
                      sprintf("%4.2f", mean(x, na.rm = TRUE)),
                      " med.=",
                      sprintf("%4.2f", median(x, na.rm = TRUE)), sep = ""))

  ### Finish all plots.
  dev.off()
}


    ### Plot histogram of true expression and prior.
    # x <- log(EPhi / mean(EPhi))
    # mu.prior <- mean(p.mcmc[nrow(p.mcmc) - 1,])
    # sd.prior <- mean(p.mcmc[nrow(p.mcmc),])
    # xx <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length = 100)
    # d.prior <- dnorm(xx, mu.prior, sd.prior)
    # tmp <- hist(x, nclass = 50, plot = FALSE)
    # ylim <- c(0, max(max(d.prior), max(tmp$counts) / sum(tmp$counts))) * 1.2
    # ylim <- c(0, 0.7)
    # hist(x, nclass = 50, freq = FALSE, ylim = ylim,
    #      xlab = "True Production Rate (log)", ylab = "Density",
    #      main = "Posterior of Expression")
    # lines(x = xx, y = d.prior, col = 2, lwd = 1.5)
    # legend(min(x), ylim[2], c("log normal"), col = 2, lwd = 1.5)


    ### For trace of prior.
    # trace <- do.call("cbind", ret$p.Mat)
    # p.names <- rownames(trace)
    # p.names[p.names == "nu.Phi"] <- "m.Phi"
    # p.names[p.names == "bsig.Phi"] <- "s.Phi"
    # xlim <- c(1, length(ret$p.Mat))

    ### Plot trace of meanlog(m).
    # i.p <- which(p.names == "m.Phi")
    # ylim <- range(trace[i.p,])
    # plot(NULL, NULL, xlim = xlim, ylim = ylim,
    #      xlab = "Iterations", ylab = p.names[i.p],
    #      main = "Trace of m.Phi")
    # lines(x = 1:length(ret$p.Mat), y = trace[i.p,])

    # Plot trace of sdlog(s).
    # i.p <- which(p.names == "s.Phi")
    # ylim <- range(trace[i.p,])
    # plot(NULL, NULL, xlim = xlim, ylim = ylim,
    #      xlab = "Iterations", ylab = p.names[i.p],
    #      main = "Trace of s.Phi")
    # lines(x = 1:length(ret$p.Mat), y = trace[i.p,])

