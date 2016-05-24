### Generate fake data from uniform distribution.

rm(list = ls())

suppressMessages(library(cubfits, quietly = TRUE))

### Preload environment and set data.
source("00-set_env.r")
set.seed(simulation$seed)
fn.in <- paste(prefix$param, "small_bInit.rda", sep = "")
if(file.exists(fn.in)){
  load(fn.in)
} else{
  stop(paste(fn.in, " is not found.", sep = ""))
}
fn.in <- paste(prefix$param, "small_train.rda", sep = "")
if(file.exists(fn.in)){
  load(fn.in)
} else{
  stop(paste(fn.in, " is not found.", sep = ""))
}
fn.in <- paste(prefix$param, "small_length.rda", sep = "")
if(file.exists(fn.in)){
  load(fn.in)
  eval(parse(text = paste("b.Init <- b.InitList.", model, sep = "")))
} else{
  stop(paste(fn.in, " is not found.", sep = ""))
}

### Generate E[Phi_g] from log normal or set to phi.Obs.
EPhi <- phi.Obs
if(simulation$EPhi){
  EPhi <- rlnorm(length(phi.Obs))
  ### For simulation based on yassour's output if available.
  # EPhi <- rlnorm(length(phi.Obs), meanlog = -sdlog^2 / 2, sdlog = sdlog)
}
names(EPhi) <- names(phi.Obs)

### Generate E[b] or set to b.Init.
Eb <- b.Init
if(simulation$Eb){
  Eb <- lapply(Eb, function(x){
                     tmp <- list()
                     tmp$coefficients <- runif(length(x$coefficients), -1, 1)
                     tmp$coef.mat <- matrix(tmp$coefficients,
                                            nrow = nrow(x$coef.mat),
                                            byrow = TRUE)
                     names(tmp$coefficients) <- names(x$coefficients)
                     tmp
                   })
}
names(Eb) <- names(b.Init)

### For Yassour 2009 only. EPhi and simu.phi.Obs are overwrote.
# GM <- apply(yassour[, -1], 1, function(x) exp(mean(log(x[x != 0]))))
# phi.Obs.all <- yassour[, -1] / sum(GM) * 15000
# phi.Obs.all[phi.Obs.all == 0] <- NA
# X <- log(as.matrix(phi.Obs.all))
# param.init <- list(K = 2, prop = c(0.95, 0.05), mu = c(-0.59, 3.11),
#                    sigma2 = c(1.40, 0.59), sigma2.e = 0.03)
# ret <- mixnorm.optim(X, K = 2, param = param.init)
# X.simu <- simu.mixnorm(length(EPhi), ret$param)
# tmp.names <- names(EPhi)
# EPhi <- X.simu$Phi
# names(EPhi) <- tmp.names
### simu.phi.Obs <- X.simu$phi.Obs
### names(simu.phi.Obs) <- tmp.names

### Generate expression.
simu.phi.Obs <- simu.phi.Obs(EPhi, sigmaW.lim = c(1, 1))

### Generate sequences.
simu.seq <- simu.orf(length(EPhi), Eb, phi.Obs = EPhi, AA.prob = AA.prob,
                     orf.length = gene.length, orf.names = names(EPhi),
                     model = model)
simu.phi <- data.frame(ORF = names(EPhi), phi.value = simu.phi.Obs,
                       true.phi = EPhi)

### Dump files.
fn.out <- paste(prefix$data, "simu_seq_", model, ".fasta", sep = "")
write.seq(simu.seq, fn.out)
fn.out <- paste(prefix$data, "simu_phi.tsv", sep = "")
write.phi.df(simu.phi, fn.out)
fn.out <- paste(prefix$data, "simu_true_", model, ".rda", sep = "")
save(EPhi, Eb, file = fn.out)

