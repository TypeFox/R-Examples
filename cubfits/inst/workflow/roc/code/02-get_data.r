### This script reads in data and summarizes them in data structures.

rm(list = ls())

suppressMessages(library(cubfits, quietly = TRUE))

### Set environment and data.
source("00-set_env.r")

### Check output directory.
if(!file.exists(prefix$data)){
  stop(paste(prefix$data, " is not found.", sep = ""))
}

### Load sequence data.
fn.in <- file.data$fasta
if(file.exists(fn.in)){
  seq.data <- read.seq(fn.in)
} else{
  stop(paste(fn.in, " is not found.", sep = ""))
}

### Load expression data.
fn.in <- file.data$tsv
if(file.exists(fn.in)){
  phi <- read.phi.df(fn.in)
  colnames(phi)[1:2] <- c("ORF", "phi")
} else{
  ### Convert to string and get SCUO.
  seq.string <- convert.seq.data.to.string(seq.data)
  y.scuo <- gen.scuo(seq.string)
  SCUO <- calc_scuo_values(y.scuo)$SCUO

  ### A fake SCUO is used.
  phi <- data.frame(ORF = names(seq.string), phi = as.double(SCUO))
}

### Check and reorder.
seq.data <- seq.data[names(seq.data) %in% phi$ORF]
phi <- phi[phi$ORF %in% names(seq.data),]
seq.data <- seq.data[order(names(seq.data))]
phi <- phi[order(phi$ORF),]

### .CF.CONF$estimate.bias.Phi may be TRUE or not, but it should not conflict
### with .CF.CONF$scale.phi.Obs. No matter .CF.CONF$scale.phi.Obs is TRUE of FALSE.
### Here, phi is for the observed measurements.
if(.CF.CONF$scale.phi.Obs){
  phi.scale <- mean(phi[, 2])
  phi[, 2] <- phi[, 2] / phi.scale
} else{
  ### .CF.CONF$estimate.bias.Phi may be TRUE
  phi.scale <- 1
}

### Convert to string and get SCUO after subsetting and reordering.
seq.string <- convert.seq.data.to.string(seq.data)
y.scuo <- gen.scuo(seq.string)
SCUO <- calc_scuo_values(y.scuo)$SCUO

### Convert to R objects. 
reu13.df.obs <- gen.reu13.df(seq.string, phi)
y <- gen.y(seq.string)
n <- gen.n(seq.string)
phi.Obs <- gen.phi.Obs(phi)
reu13.df.obs.list <- gen.reu13.list(seq.string)
y.list <- convert.y.to.list(y)
n.list <- convert.n.to.list(n)

### Comput CAI.
CAI <- calc_cai_values(y, y.list)

### Dump files.
fn.out <- paste(prefix$data, "pre_process.rda", sep = "")
list.save <- c("reu13.df.obs", "y", "n",
               "reu13.df.obs.list", "y.list", "n.list",
               "y.scuo", "phi.Obs",
               "SCUO", "phi.scale", "CAI")
save(list = list.save, file = fn.out)

print(proc.time())
