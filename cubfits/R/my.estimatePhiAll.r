### This file contains two type of methods: MLE and posterior mean for
### estimating appropriate initial values of expected gene expression phi
### for all genes given sequence information/counts.
### Main purpose of these functions is parallelization.

### - fitlist is a list beta estimation from fitMultinom(). amino acid -> beta.
### - reu13.list is an all gene list. gene -> amino acid -> codon -> position.
### - y is an all gene list. gene -> amino acid -> codon count.
### - n is an all gene list. gene -> amino acid -> total codon count.

### Get the specific function according to the options.
get.my.estimatePhiAll <- function(parallel){
  if(!any(parallel[1] %in% .CF.CT$parallel)){
    stop("parallel is not found.")
  }
  ret <- eval(parse(text = paste("my.estimatePhiAll.", parallel[1], sep = "")))
  assign("my.estimatePhiAll", ret, envir = .cubfitsEnv)
  ret
} # End of get.my.estimatePhiAll().


### For lapply.
my.estimatePhiAll.lapply <- function(fitlist, reu13.list, y.list, n.list,
    E.Phi = .CF.OP$E.Phi, lower.optim = .CF.OP$lower.optim,
    upper.optim = .CF.OP$upper.optim, lower.integrate = .CF.OP$lower.integrate,
    upper.integrate = .CF.OP$upper.integrate, control = list()){
  ret <- lapply(1:length(reu13.list),
           function(i.gene){ # i'th gene.
             .cubfitsEnv$my.estimatePhiOne(fitlist, reu13.list[[i.gene]],
                                          y.list[[i.gene]], n.list[[i.gene]],
                                          E.Phi = E.Phi,
                                          lower.optim = lower.optim,
                                          upper.optim = upper.optim,
                                          lower.integrate = lower.integrate,
                                          upper.integrate = upper.integrate,
                                          control = control)
           })
  names(ret) <- names(reu13.list)
  ret <- unlist(ret)
  ret
} # End of my.estimatePhiAll.lapply().

### For mclapply.
my.estimatePhiAll.mclapply <- function(fitlist, reu13.list, y.list, n.list,
    E.Phi = .CF.OP$E.Phi, lower.optim = .CF.OP$lower.optim,
    upper.optim = .CF.OP$upper.optim, lower.integrate = .CF.OP$lower.integrate,
    upper.integrate = .CF.OP$upper.integrate, control = list()){
  ret <- parallel::mclapply(1:length(reu13.list),
           function(i.gene){ # i'th gene.
             .cubfitsEnv$my.estimatePhiOne(fitlist, reu13.list[[i.gene]],
                                          y.list[[i.gene]], n.list[[i.gene]],
                                          E.Phi = E.Phi,
                                          lower.optim = lower.optim,
                                          upper.optim = upper.optim,
                                          lower.integrate = lower.integrate,
                                          upper.integrate = upper.integrate,
                                          control = control)
           }, mc.set.seed = FALSE, mc.preschedule = FALSE)
  names(ret) <- names(reu13.list)
  ret <- unlist(ret)
  ret
} # End of my.estimatePhiAll.mclapply().

### For task pull.
my.estimatePhiAll.task.pull <- function(fitlist, reu13.list, y.list, n.list,
    E.Phi = .CF.OP$E.Phi, lower.optim = .CF.OP$lower.optim,
    upper.optim = .CF.OP$upper.optim, lower.integrate = .CF.OP$lower.integrate,
    upper.integrate = .CF.OP$upper.integrate, control = list()){
  ret <- pbdMPI::task.pull(1:length(reu13.list),
           function(i.gene){ # i'th gene.
             .cubfitsEnv$my.estimatePhiOne(fitlist, reu13.list[[i.gene]],
                                          y.list[[i.gene]], n.list[[i.gene]],
                                          E.Phi = E.Phi,
                                          lower.optim = lower.optim,
                                          upper.optim = upper.optim,
                                          lower.integrate = lower.integrate,
                                          upper.integrate = upper.integrate,
                                          control = control)
           }, bcast = TRUE)
  names(ret) <- names(reu13.list)
  ret <- unlist(ret)
  ret
} # End of my.estimatePhiAll.task.pull().

### For pbdLapply.
my.estimatePhiAll.pbdLapply <- function(fitlist, reu13.list, y.list, n.list,
    E.Phi = .CF.OP$E.Phi, lower.optim = .CF.OP$lower.optim,
    upper.optim = .CF.OP$upper.optim, lower.integrate = .CF.OP$lower.integrate,
    upper.integrate = .CF.OP$upper.integrate, control = list()){
  ret <- pbdMPI::pbdLapply(1:length(reu13.list),
           function(i.gene){ # i'th gene.
             .cubfitsEnv$my.estimatePhiOne(fitlist, reu13.list[[i.gene]],
                                          y.list[[i.gene]], n.list[[i.gene]],
                                          E.Phi = E.Phi,
                                          lower.optim = lower.optim,
                                          upper.optim = upper.optim,
                                          lower.integrate = lower.integrate,
                                          upper.integrate = upper.integrate,
                                          control = control)
           }, pbd.mode = "spmd", bcast = TRUE)
  names(ret) <- names(reu13.list)
  ret <- unlist(ret)
  ret
} # End of my.estimatePhiAll.pbdLapply().

