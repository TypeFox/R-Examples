### This script is an approximation starting from true values since this is
### an simulation.

rm(list = ls())

suppressMessages(library(pbdMPI, quietly = TRUE))
init(set.seed = FALSE)
suppressMessages(library(cubfits, quietly = TRUE))

### Set environment.
source("00-set_env.r")
set.seed(simulation$seed)
case.name <- "wophi_true"
case.name <- paste(model, "_", case.name, sep = "")
fn.in <- paste(prefix$data, "pre_process.rda", sep = "")
load(fn.in)
fn.in <- paste(prefix$data, "simu_true_", model, ".rda", sep = "")
if(file.exists(fn.in)){
  load(fn.in)
} else{
  comm.stop(paste(fn.in, " is not found.", sep = ""))
}
phi.init.true <- EPhi
b.Init <- Eb

### Initial.
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
  phi.init.true <- phi.init.true / mean(phi.init.true)
}
ret <- cubappr(reu13.df.obs, phi.init.true, y, n,
               nIter = nIter,
               b.Init = b.Init,
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
