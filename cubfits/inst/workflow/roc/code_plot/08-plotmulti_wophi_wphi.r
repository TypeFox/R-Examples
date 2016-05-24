### This script plots comparisions across cases.

rm(list = ls())

suppressMessages(library(cubfits, quietly = TRUE))

source("00-set_env.r")
source(paste(prefix$code.plot, "u0-get_case_main.r", sep = ""))
source(paste(prefix$code, "u1-get_negsel.r", sep = ""))
source(paste(prefix$code.plot, "u2-plot_b_corr.r", sep = ""))
source(paste(prefix$code.plot, "u5-new_page.r", sep = ""))
source(paste(prefix$code.plot, "u6-adjust_focal_codon.r", sep = ""))

### Load data.
fn.in <- paste(prefix$data, "pre_process.rda", sep = "")
load(fn.in)

### Get AA and synonymous codons.
aa.names <- names(reu13.df.obs)
label <- NULL
for(i.aa in aa.names){
  tmp <- sort(unique(reu13.df.obs[[i.aa]]$Codon))
  tmp <- tmp[-length(tmp)]
  label <- c(label, paste(i.aa, tmp, sep = "."))
}

### Get possible match cases. wophi fits vs wphi fits on wphi EPhi.
match.case <- rbind(
  c("wophi_scuo", "wphi_scuo"),
  c("wphi_wophi_scuo", "wphi_scuo"),
  c("wophi_scuo", "wphi_pm")
)
match.case <- matrix(paste(model, match.case, sep = "_"), ncol = 2)

### Plot matched cases if found.
for(i.match in 1:nrow(match.case)){
  ### Initial.
  phi.pm <- list()
  phi.median <- list()
  phi.ci <- list()
  phi.pm.log10 <- list()
  phi.std.log10 <- list()
  phi.ci.log10 <- list()
  b.pm <- list()
  b.ci <- list()
  b.negsel <- list()
  b.negsel.ci <- list()
  b.negsel.label.list <- list()
  b.logmu <- list()
  b.logmu.ci <- list()
  b.logmu.label.list <- list()
  SCU.mean <- list()
  SCU.median <- list()

  ### Load data.
  for(i.case in match.case[i.match,]){
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

    ### Reassign.
    phi.pm[[i.case]] <- phi.PM
    phi.median[[i.case]] <- phi.MED
    phi.ci[[i.case]] <- phi.CI
    phi.pm.log10[[i.case]] <- phi.PM.log10
    phi.std.log10[[i.case]] <- phi.STD.log10
    phi.ci.log10[[i.case]] <- phi.CI.log10

    b.pm[[i.case]] <- b.PM
    b.ci[[i.case]] <- b.ci.PM

    ### Delta.t
    b.negsel[[i.case]] <- b.negsel.PM
    b.negsel.ci[[i.case]] <- b.negsel.ci.PM
    b.negsel.label.list[[i.case]] <- b.negsel.label

    ### log.mu
    b.logmu[[i.case]] <- b.logmu.PM
    b.logmu.ci[[i.case]] <- b.logmu.ci.PM
    b.logmu.label.list[[i.case]] <- b.logmu.label

    b <- convert.bVec.to.b(b.PM, aa.names)
    SCU.mean[[i.case]] <- calc_scu_values(b, y.list, phi.PM)
    SCU.median[[i.case]] <- calc_scu_values(b, y.list, phi.MED)
  }

  ### Check if both cases are matched.
  if(length(b.pm) != 2){
    next
  }

  ### Get negsel.
  all.names <- names(b.pm[[1]])
  id.intercept <- grep("log.mu", all.names)
  id.slop <- grep("Delta.t", all.names)


  ### Set layout.
  case.main <- paste(match.case[i.match, 1], " ",
                     match.case[i.match, 2], sep = "")
  fn.out <- paste(prefix$plot.multi, match.case[i.match, 1], "_",
                  match.case[i.match, 2], ".pdf", sep = "")
  pdf(fn.out, width = 6, height = 10)

  ### New page.
    new.page(workflow.name, case.main = case.main)

    ### Plot Delta.t.
    x <- b.negsel[[2]]
    y <- b.negsel[[1]]
    x.ci <- b.negsel.ci[[2]]
    y.ci <- b.negsel.ci[[1]]
    x.label <- b.negsel.label.list[[1]]
    plot.b.corr(x, y, x.label, x.ci = x.ci, y.ci = y.ci,
                # xlab = "With Phi", ylab = "Without Phi", main = "Delta.t",
                add.lm = TRUE, add.ci = TRUE)
    mtext(expression(paste(italic(Delta[t]), " , with ", italic(X[obs]))),
          side = 1, line = 2.5, cex = 0.8)
    mtext(expression(paste(italic(Delta[t]), " , without ", italic(X[obs]))),
          side = 2, line = 2.5, cex = 0.8)
    x.label.focal <- x.label

    ### Plot log(mu).
    x <- b.logmu[[2]]
    y <- b.logmu[[1]]
    x.ci <- b.logmu.ci[[2]]
    y.ci <- b.logmu.ci[[1]]
    x.label <- b.logmu.label.list[[1]]
    ### Adjust focal codons if needed and dispatch after adjusting.
    new.order <- adjust.focal.codon(y, x.label, x.label.focal, y.ci = y.ci)
    y <- new.order$y
    y.ci <- new.order$y.ci
    plot.b.corr(x, y, x.label.focal, x.ci = x.ci, y.ci = y.ci,
                # xlab = "With Phi", ylab = "Without Phi", main = "log(mu)",
                add.lm = TRUE, add.ci = TRUE)
    mtext(expression(paste(italic(M), " , with ", italic(X[obs]))),
          side = 1, line = 2.5, cex = 0.8)
    mtext(expression(paste(italic(M), " , without ", italic(X[obs]))),
          side = 2, line = 2.5, cex = 0.8)

    ### Overlap two histograms.
    p.1 <- hist(log10(phi.pm[[2]] / mean(phi.pm[[2]])),
                nclass = 50, plot = FALSE)
    p.2 <- hist(log10(phi.pm[[1]] / mean(phi.pm[[1]])),
                nclass = 50, plot = FALSE)
    xlim <- range(p.1$breaks, p.2$breaks)
    ylim <- range(p.1$counts, p.2$counts)
    plot(p.1, col = "#0000FF50", xlim = xlim, ylim = ylim,
         xlab = "Production Rate (log10)",
         main = "Expression (Posterior Mean)")
    plot(p.2, col = "#FF000050", xlim = xlim, ylim = ylim, add = TRUE)
    legend(xlim[1], ylim[2], c("With Phi", "Without Phi"),
           pch = c(15, 15), col = c("#0000FF50", "#FF000050"), cex = 0.8)

    ### Overlap two histograms.
    p.1 <- hist(log10(phi.median[[2]] / mean(phi.median[[2]])),
                nclass = 50, plot = FALSE)
    p.2 <- hist(log10(phi.median[[1]] / mean(phi.median[[1]])),
                nclass = 50, plot = FALSE)
    xlim <- range(p.1$breaks, p.2$breaks)
    ylim <- range(p.1$counts, p.2$counts)
    plot(p.1, col = "#0000FF50", xlim = xlim, ylim = ylim,
         xlab = "Production Rate (log10)",
         main = "Expression (Posterior Median)")
    plot(p.2, col = "#FF000050", xlim = xlim, ylim = ylim, add = TRUE)
    legend(xlim[1], ylim[2], c("With Phi", "Without Phi"),
           pch = c(15, 15), col = c("#0000FF50", "#FF000050"), cex = 0.8)

  ### New page.
    new.page(workflow.name, case.main = case.main)

    ### Plot SCU
    plotprxy(SCU.mean[[2]]$SCU, SCU.mean[[1]]$SCU,
             xlab = "With Phi SCU (log10)",
             ylab = "Without Phi SCU (log10)",
             main = "SCU (Posterior Mean)")

    ### Plot mSCU
    plotprxy(SCU.mean[[2]]$mSCU, SCU.mean[[1]]$mSCU,
             log10.x = FALSE, log10.y = FALSE,
             xlab = "With Phi mSCU",
             ylab = "Without Phi mSCU",
             main = "mSCU (Posterior mean)")

    ### Plot SCU
    plotprxy(SCU.median[[2]]$SCU, SCU.median[[1]]$SCU,
             xlab = "With Phi SCU (log10)",
             ylab = "Without Phi SCU (log10)",
             main = "SCU (Posterior Mean)")

  ### New page.
    new.page(workflow.name, case.main = case.main)

    ### For x-y plots.
    x <- log10(phi.pm[[2]] / mean(phi.pm[[2]]))
    y <- log10(phi.pm[[1]] / mean(phi.pm[[1]]))
    xlim <- ylim <- range(c(x, y))

    ### Plot prxy.
    plotprxy(10^(phi.pm.log10[[2]]), 10^(phi.pm.log10[[1]]),
             weights = 1 / phi.std.log10[[1]],
             xlim = xlim, ylim = ylim,
             xlab = "With Phi Production Rate (log10)",
             ylab = "Without Phi Production Rate (log10)",
             main = "Expression (Posterior log10 Mean)")

    ### Plot prxy with outliers labeled.
    plotprxy(10^(phi.pm.log10[[2]]), 10^(phi.pm.log10[[1]]),
             y.ci = 10^(phi.ci.log10[[1]]),
             weights = 1 / phi.std.log10[[1]],
             xlim = xlim, ylim = ylim,
             xlab = "With Phi Production Rate (log10)",
             ylab = "Without Phi Production Rate (log10)",
             main = "Expression (Posterior log10 Mean)")

    ### Plot prxy.
    plotprxy(phi.median[[2]], phi.median[[1]],
             weights = 1 / phi.std.log10[[1]],
             xlim = xlim, ylim = ylim,
             xlab = "With Phi Production Rate (log10)",
             ylab = "Without Phi Production Rate (log10)",
             main = "Expression (Posterior Median)")

    ### Plot prxy with outliers labeled.
    plotprxy(phi.median[[2]], phi.median[[1]],
             y.ci = phi.ci[[1]],
             weights = 1 / phi.std.log10[[1]],
             xlim = xlim, ylim = ylim,
             xlab = "With Phi Production Rate (log10)",
             ylab = "Without Phi Production Rate (log10)",
             main = "Expression (Posterior Median)")

    ### Plot prxy.
    plotprxy(phi.pm[[2]], phi.pm[[1]],
             weights = 1 / phi.std.log10[[1]],
             xlim = xlim, ylim = ylim,
             xlab = "With Phi Production Rate (log10)",
             ylab = "Without Phi Production Rate (log10)",
             main = "Expression (Posterior Mean)")

    ### Plot prxy with outliers labeled.
    plotprxy(phi.pm[[2]], phi.pm[[1]],
             y.ci = phi.ci[[1]],
             weights = 1 / phi.std.log10[[1]],
             xlim = xlim, ylim = ylim,
             xlab = "With Phi Production Rate (log10)",
             ylab = "Without Phi Production Rate (log10)",
             main = "Expression (Posterior Mean)")

  ### New page.
    new.page(workflow.name, case.main = case.main)

    ### Add qqplot.
    x <- log10(phi.pm[[2]] / mean(phi.pm[[2]]))
    y <- log10(phi.pm[[1]] / mean(phi.pm[[1]]))
    xlim <- ylim <- range(c(x, y))
    qqplot(x, y,
           xlim = xlim, ylim = ylim,
           xlab = "With Phi Production Rate (log10)",
           ylab = "Without Phi Production Rate (log10)",
           main = "Q-Q (Posterior Mean)",
           pch = 20, cex = 0.6)
    abline(a = 0, b = 1, col = 4, lty = 2)

    ### Add qqplot for non-log.
    x <- phi.pm[[2]] / mean(phi.pm[[2]])
    y <- phi.pm[[1]] / mean(phi.pm[[1]])
    xlim <- ylim <- range(c(x, y))
    qqplot(x, y,
           xlim = xlim, ylim = ylim,
           xlab = "With Phi Production Rate",
           ylab = "Without Phi Production Rate",
           main = "Q-Q (Posterior Mean)",
           pch = 20, cex = 0.6)
    abline(a = 0, b = 1, col = 4, lty = 2)

    ### Add qqplot.
    x <- log10(phi.median[[2]] / mean(phi.median[[2]]))
    y <- log10(phi.median[[1]] / mean(phi.median[[1]]))
    xlim <- ylim <- range(c(x, y))
    qqplot(x, y,
           xlim = xlim, ylim = ylim,
           xlab = "With Phi Production Rate (log10)",
           ylab = "Without Phi Production Rate (log10)",
           main = "Q-Q (Posterior Median)",
           pch = 20, cex = 0.6)
    abline(a = 0, b = 1, col = 4, lty = 2)

    ### Add qqplot for non-log.
    x <- phi.median[[2]] / mean(phi.median[[2]])
    y <- phi.median[[1]] / mean(phi.median[[1]])
    xlim <- ylim <- range(c(x, y))
    qqplot(x, y,
           xlim = xlim, ylim = ylim,
           xlab = "With Phi Production Rate",
           ylab = "Without Phi Production Rate",
           main = "Q-Q (Posterior Median)",
           pch = 20, cex = 0.6)
    abline(a = 0, b = 1, col = 4, lty = 2)

  dev.off()
}
