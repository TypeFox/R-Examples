rm(list = ls())

suppressMessages(library(cubfits, quietly = TRUE))

### Load environment and set data.
source("00-set_env.r")
source(paste(prefix$code.plot, "u0-get_case_main.r", sep = ""))

### Trace each run.
for(i.case in case.names){
  ### All mcmc outputs.
  fn.in <- paste(prefix$output, i.case, "/output_mcmc.rda", sep = "")
  if(!file.exists(fn.in)){
    cat("File not found: ", fn.in, "\n", sep = "")
    next
  }
  load(fn.in)

  x <- 1:length(ret$p.Mat)
  xlim <- range(x)

  trace <- do.call("cbind", ret$p.Mat)
  if(nrow(trace) == 4){
    p.names <- c("sigmaW", "m.Phi", "s.Phi", "bias.Phi")
  } else if(nrow(trace) == 3){
    p.names <- c("sigmaW", "m.Phi", "s.Phi")
  } else if(nrow(trace) == 2){
    p.names <- c("m.Phi", "s.Phi")
  } else{
    cat("Not log normal model?\n")
    next
  }

  ### Plot priors.
  for(i.p in 1:length(p.names)){
    fn.out <- paste(prefix$plot.trace, "prior_", p.names[i.p], "_", i.case,
                    ".pdf", sep = "")
    pdf(fn.out, width = 6, height = 4)
      ylim <- range(trace[i.p,])
      plot(NULL, NULL, xlim = xlim, ylim = ylim,
           xlab = "Iterations", ylab = p.names[i.p])
      mtext(paste(workflow.name, ", ", get.case.main(i.case, model), sep = ""),
            line = 3, cex = 0.6)
      mtext(date(), line = 2.5, cex = 0.4)
      lines(x = x, y = trace[i.p,])
    dev.off()
  }
}

