### This compute log posterior on all amino acids.

### Get the specific function according to the options.
get.my.logdmultinomCodAllR <- function(parallel){
  if(!any(parallel[1] %in% .CF.CT$parallel)){
    stop("parallel method is not found.")
  }
  ret <- eval(parse(text = paste("my.logdmultinomCodAllR.", parallel[1],
                                 sep = "")))
  assign("my.logdmultinomCodAllR", ret, envir = .cubfitsEnv)
  ret
} # End of get.my.logdmultinomCodAllR().


### For lapply.
my.logdmultinomCodAllR.lapply <- function(b, phi, y, n, reu13.df = NULL){
  ### Returns log posterior of codon draws for all amino acids.
  ### For each element, it is a vector of length "# of genes".
  lpclist <- lapply(1:length(y),
               function(i.aa){ # i'th amino acid.
                 .cubfitsEnv$my.logdmultinomCodOne(
                   b[[i.aa]], phi, y[[i.aa]], n[[i.aa]], vec = TRUE,
                   reu13.df.aa = reu13.df[[i.aa]])
               })

  ### Return posterior which is the sum of all amino acids.
  rowSums(do.call("cbind", lpclist))
} # End of my.logdmultinomCodAllR.lapply().

### For mclapply.
my.logdmultinomCodAllR.mclapply <- function(b, phi, y, n, reu13.df = NULL){
  ### Returns log posterior of codon draws for all amino acids.
  ### For each element, it is a vector of length "# of genes".
  lpclist <- parallel::mclapply(1:length(y),
               function(i.aa){ # i'th amino acid.
                 .cubfitsEnv$my.logdmultinomCodOne(
                   b[[i.aa]], phi, y[[i.aa]], n[[i.aa]], vec = TRUE,
                   reu13.df.aa = reu13.df[[i.aa]])
               }, mc.set.seed = FALSE, mc.preschedule = FALSE)

  ### Return posterior which is the sum of all amino acids.
  rowSums(do.call("cbind", lpclist))
} # End of my.logdmultinomCodAllR.mclapply().

### For task.pull.
my.logdmultinomCodAllR.task.pull <- function(b, phi, y, n, reu13.df = NULL){
  ### pbdLapply is good enough.
  my.logdmultinomCodAllR.pbdLapply(b, phi, y, n, reu13.df = reu13.df)
} # End of my.logdmultinomCodAllR.pbdLapply().

### For pbdLapply.
my.logdmultinomCodAllR.pbdLapply <- function(b, phi, y, n, reu13.df = NULL){
  ### No need to bring everything back to master.
  lpclist <- pbdMPI::pbdLapply(1:length(y),
               function(i.aa){ # i'th amino acid.
                 .cubfitsEnv$my.logdmultinomCodOne(
                   b[[i.aa]], phi, y[[i.aa]], n[[i.aa]], vec = TRUE,
                   reu13.df.aa = reu13.df[[i.aa]])
               }, pbd.mode = "spmd", bcast = FALSE)

  ### Check non-empty list first.
  ret <- rep(0.0, length(phi))
  if(length(lpclist) > 0){
    ret <- as.double(rowSums(do.call("cbind", lpclist)))
  }
  ret <- pbdMPI::spmd.allreduce.double(ret, double(length(phi)))
  ret
} # End of my.logdmultinomCodAllR.pbdLapply().
