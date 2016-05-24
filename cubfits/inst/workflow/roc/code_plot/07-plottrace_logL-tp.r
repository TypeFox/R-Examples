rm(list = ls())

suppressMessages(library(cubfits, quietly = TRUE))

### Load environment and set data.
suppressMessages(library(pbdMPI, quietly = TRUE))
init(set.seed = FALSE)
source("00-set_env.r")
source(paste(prefix$code.plot, "u0-get_case_main.r", sep = ""))
fn.in <- paste(prefix$data, "pre_process.rda", sep = "")
load(fn.in)

### Initial.
init.function(model = model)

### Trace each run.
# for(i.case in case.names){
all.jobs <- function(i.job){
  i.case <- case.names[i.job]

  ### All mcmc outputs.
  fn.in <- paste(prefix$output, i.case, "/output_mcmc.rda", sep = "")
  if(!file.exists(fn.in)){
    # cat("File not found: ", fn.in, "\n", sep = "")
    # next
    stop(paste("File not found: ", fn.in, "\n", sep = ""))
  }
  cat(i.job, ": ", i.case, ", load: ", fn.in, "\n", sep = "")
  load(fn.in)

  ### Since my.appr() doesn't have phi.Mat, but have phi.pred.Mat
  if(is.null(ret[["phi.Mat"]])){
    ret$phi.Mat <- ret$phi.pred.Mat
  }

  if("sigmaW" %in% rownames(ret$p.Mat[[1]])){
    logL <- lapply(1:length(ret$phi.Mat),
                   function(i.iter){
                     xx <- ret$phi.Mat[[i.iter]]
                     b.Init <- convert.bVec.to.b(ret$b.Mat[[i.iter]],
                                                names(reu13.df.obs),
                                                model = model)
                     b.Init <- lapply(b.Init, function(B) B$coefficients)
                     sigmaWsq <- ret$p.Mat[[i.iter]][1]^2
                     tmp <- .cubfitsEnv$my.logLAll(xx, phi.Obs, y, n, b.Init,
                                                   sigmaWsq,
                                                   reu13.df = reu13.df.obs)
                     sum(tmp)
                  })
  } else{
    logL <- lapply(1:length(ret$phi.Mat),
                   function(i.iter){
                     xx <- ret$phi.Mat[[i.iter]]
                     b.Init <- convert.bVec.to.b(ret$b.Mat[[i.iter]],
                                                names(reu13.df.obs),
                                                model = model)
                     b.Init <- lapply(b.Init, function(B) B$coefficients)
                     tmp <- .cubfitsEnv$my.logLAllPred(xx, y, n, b.Init,
                                                       reu13.df = reu13.df.obs)
                     sum(tmp)
                  })
  }

  ### Load logL mean results.
  fn.in <- paste(prefix$subset, i.case, "_PM.rda", sep = "")
  load(fn.in)
  # fn.in <- paste(prefix$subset, i.case, "_PM_scaling.rda", sep = "")
  # load(fn.in)

  b.Init <- convert.bVec.to.b(b.PM, names(reu13.df.obs), model = model)
  b.Init <- lapply(b.Init, function(B) B$coefficients)

  if("sigmaW" %in% rownames(p.PM)){
    sigmaWsq <- p.PM[1]^2
    tmp <- .cubfitsEnv$my.logLAll(phi.PM, phi.Obs, y, n, b.Init, sigmaWsq,
                                  reu13.df = reu13.df.obs)
  } else{
    tmp <- .cubfitsEnv$my.logLAllPred(phi.PM, y, n, b.Init,
                                      reu13.df = reu13.df.obs)
  }
  logL.PM <- sum(tmp)

  x <- 1:length(ret$phi.Mat)
  xlim <- range(x)
  ylim <- range(range(logL), logL.PM)

  ### Trace of logL.
  fn.out <- paste(prefix$plot.trace, "logL_", i.case, "_nps.pdf", sep = "")
  cat(i.job, ": ", i.case, ", plot: ", fn.out, "\n", sep = "")
  pdf(fn.out, width = 6, height = 4)
    plot(NULL, NULL, xlim = xlim, ylim = ylim,
         main = paste("PM of logL = ", sprintf("%.4f", logL.PM), sep = ""),
         xlab = "Iterations", ylab = "Prop. logL")
    mtext(paste(workflow.name, ", ", get.case.main(i.case, model), sep = ""),
          line = 3, cex = 0.6)
    mtext(date(), line = 2.5, cex = 0.4)
    lines(x = x, y = logL)
    abline(h = logL.PM, col = 2)
  dev.off()

  ### Dump logL.
  fn.out <- paste(prefix$subset, "trace_logL_", i.case, ".rda", sep = "")
  cat(i.job, ": ", i.case, ", dump: ", fn.out, "\n", sep = "")
  save(logL, logL.PM, file = fn.out)

  return(c(comm.rank(), i.job))
} # End of all.jobs().

ret <- task.pull(1:length(case.names), all.jobs, try.silent = TRUE)
comm.print(ret)
finalize()
