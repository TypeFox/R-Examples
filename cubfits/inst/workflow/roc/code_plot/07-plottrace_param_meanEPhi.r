rm(list = ls())

suppressMessages(library(cubfits, quietly = TRUE))

### Load environment and set data.
source("00-set_env.r")
source(paste(prefix$code.plot, "u0-get_case_main.r", sep = ""))
source(paste(prefix$code.plot, "u3-plot_trace.r", sep = ""))
fn.in <- paste(prefix$data, "init_", model, ".rda", sep = "")
load(fn.in)
coef.names <- cubfits:::get.my.coefnames(model)

### Load true Phi.
fn.in <- paste(prefix$data, "simu_true_", model, ".rda", sep = "")
if(file.exists(fn.in)){
  load(fn.in)
  b.Init <- lapply(Eb, function(x){
                        tmp <- x$coefficients
                        names(tmp) <- rep(coef.names,
                                          eacho = length(tmp) /
                                                  length(coef.names))
                        tmp
                      })
  b.Init <- do.call("c", b.Init)

  ### Get true values and scale accordingly.
  all.names <- names(b.Init)
  id.slop <- grep("Delta.t", all.names)
  scale.EPhi <- mean(EPhi)
  b.Init.scaled <- b.Init 
  b.Init.scaled[id.slop] <- b.Init.scaled[id.slop] * scale.EPhi
} else{
  b.Init <- NULL
  b.Init.scaled <- NULL
}

### Check split.S.
names.aa <- names(fitlist)
split.S <- FALSE
if("Delta.t" %in% fitlist){
  if("Z" %in% names.aa){
    split.S <- TRUE
  }
  if(ncol(fitlist[["S"]]$coef.mat) == 3){
    split.S <- TRUE
  }
}
names.b <- lapply(1:length(names.aa),
             function(i.aa){
               tmp <- fitlist[[i.aa]]$coefficients
               tmp <- rep(coef.names, each = length(tmp) / length(coef.names))
               cbind(names.aa[i.aa], tmp)
             })
names.b <- do.call("rbind", names.b)
id.intercept <- grep("log.mu", names.b[, 2])
id.slop <- grep("Delta.t", names.b[, 2])

### Trace each run.
for(i.case in case.names){
  ### All mcmc outputs.
  fn.in <- paste(prefix$output, i.case, "/output_mcmc.rda", sep = "")
  if(!file.exists(fn.in)){
    cat("File not found: ", fn.in, "\n", sep = "")
    next
  }
  load(fn.in)

  ### Since my.appr() doesn't have phi.Mat, but have phi.pred.Mat
  if(is.null(ret[["phi.Mat"]])){
    ret$phi.Mat <- ret$phi.pred.Mat
  }

  ### For original plots.
  ret.b.Mat <- ret$b.Mat
  ret.phi.Mat <- ret$phi.Mat

  ### For logmu.
  fn.out <- paste(prefix$plot.trace, "param_logmu_", i.case, ".pdf", sep = "")
  pdf(fn.out, width = 12, height = 11)
    plottrace.param(ret.b.Mat, names.b, names.aa, id.intercept,
                    workflow.name, i.case, model,
                    b.Init = b.Init, param = "logmu")
  dev.off()

  ### For deltat.
  fn.out <- paste(prefix$plot.trace, "param_deltat_", i.case, ".pdf", sep = "")
  pdf(fn.out, width = 12, height = 11)
    plottrace.param(ret.b.Mat, names.b, names.aa, id.slop,
                    workflow.name, i.case, model,
                    b.Init = b.Init, param = "deltat")
  dev.off()

  ### For mean EPhi.
  fn.out <- paste(prefix$plot.trace, "meanEPhi_", 
                  i.case, ".pdf", sep = "")
  pdf(fn.out, width = 6, height = 4)
    plottrace.meanEPhi(ret.phi.Mat,
                       workflow.name, i.case, model)
  dev.off()

  ### For scaled plots.
  ret.b.Mat.scaled <- lapply(1:length(ret.b.Mat),
                        function(i){
                          tmp <- ret.b.Mat[[i]]
                          tmp[id.slop] <- tmp[id.slop] * mean(ret.phi.Mat[[i]])
                          tmp
                        })

  ### For deltat.
  fn.out <- paste(prefix$plot.trace, "scaled_param_deltat_",
                  i.case, ".pdf", sep = "")
  pdf(fn.out, width = 12, height = 11)
    plottrace.param(ret.b.Mat.scaled, names.b, names.aa, id.slop,
                    paste(workflow.name, ", scaled", sep = ""),
                    i.case, model,
                    b.Init = b.Init.scaled, param = "deltat")
  dev.off()
}
