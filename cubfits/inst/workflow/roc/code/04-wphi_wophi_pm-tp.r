### This is similar to "04-wophi_bInit-tp.r", but the b.Init is provided by
### the poster means of "04-wphi_pm-tp.r".

rm(list = ls())

suppressMessages(library(pbdMPI, quietly = TRUE))
init(set.seed = FALSE)
suppressMessages(library(cubfits, quietly = TRUE))

### Set environment.
source("00-set_env.r")
set.seed(simulation$seed)
case.name <- "wphi_wophi_pm"
case.name <- paste(model, "_", case.name, sep = "")

### Check output directory.
fn.out <- paste(prefix$output, case.name, sep = "")
if(!file.exists(fn.out)){
  comm.stop(paste(fn.out, " is not found.", sep = ""))
}

### Load data.
fn.in <- paste(prefix$data, "pre_process.rda", sep = "")
load(fn.in)

### Load b.PM and phi.PM from previous summarized MCMC outputs.
fn.in <- paste(prefix$subset, model, "_wphi_pm_PM.rda", sep = "")
load(fn.in)
# b.Init <- convert.bVec.to.b(b.PM, names(reu13.df.obs), model = model)
phi.init.PM <- phi.PM

### Change initial for fitsappr.
nIter <- run.info$nIter

### For configuration.
.CF.DP$dump <- run.info$dump
.CF.DP$prefix.dump <- run.info$prefix.dump
.CF.CT$parallel <- run.info$parallel
.CF.CONF$estimate.bias.Phi <- FALSE
if(.CF.CT$type.p[1] == "lognormal_bias"){
  .CF.CT$type.p[1] <- "lognormal_RW"
}

### Run.
if(.CF.CONF$scale.phi.Obs || .CF.CONF$estimate.bias.Phi){
  phi.init.PM <- phi.init.PM / mean(phi.init.PM)
}
ret <- cubappr(reu13.df.obs, phi.init.PM, y, n,
               nIter = nIter,
               # b.Init = b.Init,
               p.nclass = p.nclass,
               model = model, verbose = TRUE, report = 10)

### Dump results.
if(comm.rank() == 0){
  ret.time <- proc.time()
  print(ret.time)

  fn.out <- paste(prefix$output, case.name, "/output_mcmc.rda", sep = "")
  save(list = c("nIter", "ret", "ret.time"),
       file = fn.out)

  fn.out <- paste(prefix$output, case.name, "/output_env.rda", sep = "")
  save(list = ls(envir = .cubfitsEnv),
       file = fn.out, envir = .cubfitsEnv)

  warnings()
}

finalize()
