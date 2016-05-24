### This script collect the poster means of MCMC runs.

rm(list = ls())

suppressMessages(library(cubfits, quietly = TRUE))
suppressMessages(library(pbdMPI, quietly = TRUE))
init(set.seed = FALSE)
source("00-set_env.r")
source(paste(prefix$code, "u1-get_negsel.r", sep = ""))

### Load data.
fn.in <- paste(prefix$data, "pre_process.rda", sep = "")
load(fn.in)

### Get AA and synonymous codons.
aa.names <- names(reu13.df.obs)
coef.names <- cubfits:::get.my.coefnames(model)
b.label <- NULL
b.names <- NULL
for(i.aa in aa.names){
  tmp <- sort(unique(reu13.df.obs[[i.aa]]$Codon))
  tmp <- tmp[-length(tmp)]
  b.label <- c(b.label, paste(i.aa, tmp, sep = "."))
  b.names <- c(b.names, rep(coef.names, each = length(tmp)))
}

### Get all cases.
# for(i.case in case.names){
all.jobs <- function(i.job){
  i.case <- case.names[i.job]

  ### Check first.
  fn.in <- paste(prefix$output, i.case, "/output_mcmc.rda", sep = "")
  if(!file.exists(fn.in)){
    # cat("File not found: ", fn.in, "\n", sep = "")
    # next
    stop(paste("File not found: ", fn.in, "\n", sep = ""))
  }

  ### Load MCMC output.
  cat(i.job, ": ", i.case, ", load: ", fn.in, "\n", sep = "")
  load(fn.in)

  ### Obtain the last 5000 iterations.
  b.mcmc <- do.call("cbind", ret$b.Mat[range$subset])
  rownames(b.mcmc) <- b.names

  p.mcmc <- do.call("cbind", ret$p.Mat[range$subset])
  rownames(p.mcmc) <- names(ret$p.Mat[[1]])

  ### This is a risky matching by names since R use $ for lazy matching,
  ### though it is fine for "phi.pred.Mat" when "phi.Mat" is missing.
  # phi.mcmc <- do.call("cbind", ret$phi.Mat[range$subset])
  # rownames(phi.mcmc) <- names(ret$phi.Mat[[1]])
  ### Since my.appr() cases don't have phi.Mat but have phi.pred.Mat
  if(is.null(ret[["phi.Mat"]])){
    ret$phi.Mat <- ret$phi.pred.Mat
  }
  phi.mcmc <- do.call("cbind", ret[["phi.Mat"]][range$subset])
  rownames(phi.mcmc) <- names(ret[["phi.Mat"]][[1]])

  ### Dump posterior distributions.
  fn.out <- paste(prefix$subset, i.case, ".rda", sep = "")
  cat(i.job, ": ", i.case, ", dump: ", fn.out, "\n", sep = "")
  save(b.mcmc, p.mcmc, phi.mcmc, file = fn.out)

  ### Find indices.
  all.names <- rownames(b.mcmc)
  id.slop <- grep("Delta.t", all.names)
  id.intercept <- grep("log.mu", all.names)

### Original scale. ###
  ### Obtain posterior means.
  b.PM <- rowMeans(b.mcmc)
  b.STD <- apply(b.mcmc, 1, sd)
  b.ci.PM <- t(apply(b.mcmc, 1, quantile, prob = ci.prob))
  b.MED <- apply(b.mcmc, 1, median)

  p.PM <- rowMeans(p.mcmc)
  p.STD <- apply(p.mcmc, 1, sd)
  p.CI <- t(apply(p.mcmc, 1, quantile, prob = ci.prob))
  p.MED <- apply(p.mcmc, 1, median)

  phi.PM <- rowMeans(phi.mcmc)
  phi.STD <- apply(phi.mcmc, 1, sd)
  phi.CI <- t(apply(phi.mcmc, 1, quantile, prob = ci.prob))
  phi.MED <- apply(phi.mcmc, 1, median)

  phi.mcmc.log10 <- log10(phi.mcmc)
  phi.PM.log10 <- rowMeans(phi.mcmc.log10)
  phi.STD.log10 <- apply(phi.mcmc.log10, 1, sd)
  phi.CI.log10 <- t(apply(phi.mcmc.log10, 1, quantile, prob = ci.prob))

  ### Negative selection.
  ret <- get.negsel(b.PM, id.intercept, id.slop, aa.names, b.label,
                    b.ci.PM = b.ci.PM)
  ### Delta.t
  b.negsel.PM <- ret$b.negsel.PM
  b.negsel.ci.PM <- ret$b.negsel.ci.PM
  b.negsel.label <- ret$b.negsel.label
  ### log.mu
  b.logmu.PM <- ret$b.logmu.PM
  b.logmu.ci.PM <- ret$b.logmu.ci.PM
  b.logmu.label <- ret$b.logmu.label

  ### Negative selection by median.
  ret <- get.negsel(b.MED, id.intercept, id.slop, aa.names, b.label,
                    b.ci.PM = b.ci.PM)
  ### Delta.t
  b.negsel.MED <- ret$b.negsel.PM
  ### log.mu
  b.logmu.MED <- ret$b.logmu.PM

  ### Dump summarized results.
  fn.out <- paste(prefix$subset, i.case, "_PM.rda", sep = "")
  cat(i.job, ": ", i.case, ", dump: ", fn.out, "\n", sep = "")
  save(b.PM, b.STD, b.ci.PM, b.MED, b.label,
       p.PM, p.STD, p.CI, p.MED,
       phi.PM, phi.STD, phi.CI, phi.MED,
       phi.PM.log10, phi.STD.log10, phi.CI.log10,
       b.negsel.PM, b.negsel.ci.PM, b.negsel.MED, b.negsel.label,
       b.logmu.PM, b.logmu.ci.PM, b.logmu.MED, b.logmu.label,
       file = fn.out)

  ### Thinning.
  # id.save <- (range$subset %% range$thinning) == 1
  # b.mcmc <- b.mcmc[, id.save]
  # p.mcmc <- p.mcmc[, id.save]
  # phi.mcmc <- phi.mcmc[, id.save]

  ### Dump posterior distributions.
  # fn.out <- paste(prefix$subset, i.case, "_thinning.rda", sep = "")
  # cat(i.job, ": ", i.case, ", dump: ", fn.out, "\n", sep = "")
  # save(b.mcmc, p.mcmc, phi.mcmc, file = fn.out)

  ### Obtain posterior means.
  # b.PM <- rowMeans(b.mcmc)
  # p.PM <- rowMeans(p.mcmc)
  # phi.PM <- rowMeans(phi.mcmc)

  ### Dump posterior means.
  # fn.out <- paste(prefix$subset, i.case, "_PM_thinning.rda", sep = "")
  # cat(i.job, ": ", i.case, ", dump: ", fn.out, "\n", sep = "")
  # save(b.PM, p.PM, phi.PM, file = fn.out)

### Post-scaling. ###
  ### Scaling.
  scale.EPhi <- colMeans(phi.mcmc)
  phi.mcmc <- t(t(phi.mcmc) / scale.EPhi)
  b.mcmc[id.slop,] <- t(t(b.mcmc[id.slop,]) * scale.EPhi)
  b.PM <- rowMeans(b.mcmc)
  b.STD <- apply(b.mcmc, 1, sd)
  b.ci.PM <- t(apply(b.mcmc, 1, quantile, prob = ci.prob))
  b.MED <- apply(b.mcmc, 1, median)

  phi.PM <- rowMeans(phi.mcmc)
  phi.STD <- apply(phi.mcmc, 1, sd)
  phi.CI <- t(apply(phi.mcmc, 1, quantile, prob = ci.prob))
  phi.MED <- apply(phi.mcmc, 1, median)

  phi.mcmc.log10 <- log10(phi.mcmc)
  phi.PM.log10 <- rowMeans(phi.mcmc.log10)
  phi.STD.log10 <- apply(phi.mcmc.log10, 1, sd)
  phi.CI.log10 <- t(apply(phi.mcmc.log10, 1, quantile, prob = ci.prob))

  ### Negative selection.
  ret <- get.negsel(b.PM, id.intercept, id.slop, aa.names, b.label,
                    b.ci.PM = b.ci.PM)
  ### Delta.t
  b.negsel.PM <- ret$b.negsel.PM
  b.negsel.ci.PM <- ret$b.negsel.ci.PM
  b.negsel.label <- ret$b.negsel.label
  ### log.mu
  b.logmu.PM <- ret$b.logmu.PM
  b.logmu.ci.PM <- ret$b.logmu.ci.PM
  b.logmu.label <- ret$b.logmu.label

  ### Negative selection by median.
  ret <- get.negsel(b.MED, id.intercept, id.slop, aa.names, b.label,
                    b.ci.PM = b.ci.PM)
  ### Delta.t
  b.negsel.MED <- ret$b.negsel.PM
  ### log.mu
  b.logmu.MED <- ret$b.logmu.PM

  ### Dump summarized results.
  fn.out <- paste(prefix$subset, i.case, "_PM_scaling.rda", sep = "")
  cat(i.job, ": ", i.case, ", dump: ", fn.out, "\n", sep = "")
  save(b.PM, b.STD, b.ci.PM, b.MED, b.label,
       phi.PM, phi.STD, phi.CI, phi.MED,
       phi.PM.log10, phi.STD.log10, phi.CI.log10,
       b.negsel.PM, b.negsel.ci.PM, b.negsel.MED, b.negsel.label,
       b.logmu.PM, b.logmu.ci.PM, b.logmu.MED, b.logmu.label,
       file = fn.out)

  return(c(comm.rank(), i.job))
} # End of all.jobs().

ret <- task.pull(1:length(case.names), all.jobs, try.silent = TRUE)
comm.print(ret)
finalize()
