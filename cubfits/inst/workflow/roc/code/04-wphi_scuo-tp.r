rm(list = ls())

suppressMessages(library(pbdMPI, quietly = TRUE))
init(set.seed = FALSE)
suppressMessages(library(cubfits, quietly = TRUE))

### Set environment.
source("00-set_env.r")
set.seed(simulation$seed)
case.name <- "wphi_scuo"
case.name <- paste(model, "_", case.name, sep = "")

### Check output directory.
fn.out <- paste(prefix$output, case.name, sep = "")
if(!file.exists(fn.out)){
  comm.stop(paste(fn.out, " is not found.", sep = ""))
}

### Load data.
fn.in <- paste(prefix$data, "pre_process.rda", sep = "")
load(fn.in)
fn.in <- paste(prefix$data, "init_", model, ".rda", sep = "")
load(fn.in)

### Initial.
nIter <- run.info$nIter

### For configuration.
.CF.DP$dump <- run.info$dump
.CF.DP$prefix.dump <- run.info$prefix.dump
.CF.CT$parallel <- run.info$parallel

### Run.
if(.CF.CT$model.Phi[1] == "logmixture"){
  phi.init.SCUO <- phi.init.SCUO.emp
}
if(.CF.CONF$scale.phi.Obs){
  phi.Obs <- phi.Obs / mean(phi.Obs)
}
if(.CF.CONF$scale.phi.Obs || .CF.CONF$estimate.bias.Phi){
  phi.init.SCUO <- phi.init.SCUO / mean(phi.init.SCUO)
}
ret <- cubfits(reu13.df.obs, phi.Obs, y, n,
               nIter = nIter,
               p.nclass = p.nclass,
               phi.Init = phi.init.SCUO,
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

