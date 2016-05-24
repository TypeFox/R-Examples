suppressMessages(library(cubfits, quietly = TRUE))

source("00-set_env.r")
source(paste(prefix$code.plot, "u0-get_case_main.r", sep = ""))

### Load true Phi.
fn.in <- paste(prefix$data, "simu_true_", model, ".rda", sep = "")
if(file.exists(fn.in)){
  load(fn.in)
} else{
  stop(paste(fn.in, " is not found.", sep = ""))
}

### Pre processed phi.Obs.
fn.in <- paste(prefix$data, "pre_process.rda", sep = "")
load(fn.in)

### True SCU.
b.Init <- convert.b.to.bVec(Eb)
all.names <- names(b.Init)
id.slop <- grep("Delta.t", all.names)
scale.EPhi <- mean(EPhi)
b.Init[id.slop] <- b.Init[id.slop] * scale.EPhi
EPhi <- EPhi / scale.EPhi
Eb <- convert.bVec.to.b(b.Init, names(reu13.df.obs))
SCU.true <- calc_scu_values(Eb, y.list, EPhi)

for(i.case in case.names){
  ### Subset of mcmc output.
  fn.in <- paste(prefix$subset, i.case, "_PM.rda", sep = "")
  if(!file.exists(fn.in)){
    cat("File not found: ", fn.in, "\n", sep = "")
    next
  }
  load(fn.in)
  # fn.in <- paste(prefix$subset, i.case, "_PM_scaling.rda", sep = "")
  # if(!file.exists(fn.in)){
  #   cat("File not found: ", fn.in, "\n", sep = "")
  #   next
  # }
  # load(fn.in)

  b <- convert.bVec.to.b(b.PM, names(reu13.df.obs))
  SCU <- calc_scu_values(b, y.list, phi.PM)

  ### Plot SCU.
  fn.out <- paste(prefix$plot.single,
                  "scu_true_", i.case, "_nps.pdf", sep = "")
  pdf(fn.out, width = 5, height = 5)
    plotprxy(SCU.true$SCU, SCU$SCU,
             xlab = "True SCU (log10)",
             ylab = "Predicted SCU (log10)",
             main = "SCU (Posterior Mean)")
    mtext(paste(workflow.name, ", ", get.case.main(i.case, model), sep = ""),
          line = 3, cex = 0.6)
  dev.off()

  ### Plot mSCU.
  fn.out <- paste(prefix$plot.single,
                  "mscu_true_", i.case, "_nps.pdf", sep = "")
  pdf(fn.out, width = 5, height = 5)
    plotprxy(SCU.true$mSCU, SCU$mSCU,
             log10.x = FALSE, log10.y = FALSE,
             xlab = "True mSCU",
             ylab = "Predicted mSCU",
             main = "mSCU (Posterior Mean)")
    mtext(paste(workflow.name, ", ", get.case.main(i.case, model), sep = ""),
          line = 3, cex = 0.6)
  dev.off()
}
