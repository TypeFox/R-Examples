### Fit multinomial logistic regression and find Delta.t and log(mu) (selection and mutation
### effects) based on vglm() in VGAM package.
###
### These functions are for all amino acids.
###
### Main purpose of these functions is parallelization.

# Get the specific function according to the options.
get.my.fitMultinomAll <- function(parallel){
  if(!any(parallel[1] %in% .CF.CT$parallel)){
    stop("parallel method is not found.")
  }
  ret <- eval(parse(text = paste("my.fitMultinomAll.", parallel[1], sep = "")))
  assign("my.fitMultinomAll", ret, envir = .cubfitsEnv)
  ret
} # End of get.my.fitMultinomAll().


### For lapply.
my.fitMultinomAll.lapply <- function(reu13.df, phi, y, n, phi.new = NULL,
    coefstart = NULL){
  ### since vglm change seeds within it's call, backup first and restore later.
  if(exists(".Random.seed", envir = .GlobalEnv)){
    .GlobalEnv$.Random.seed.org <- .GlobalEnv$.Random.seed
  }

  ret <- lapply(1:length(reu13.df),
           function(i.aa){ # i'th amino acid.
             .cubfitsEnv$my.fitMultinomOne(reu13.df[[i.aa]], phi, y[[i.aa]],
                                           n[[i.aa]], phi.new = phi.new,
                                           coefstart = coefstart[[i.aa]])
           })

  if(exists(".Random.seed.org", envir = .GlobalEnv)){
    .GlobalEnv$.Random.seed <- .GlobalEnv$.Random.seed.org
  }

  ret
} # End of my.fitmultinomAll.lapply().

### For mclapply.
my.fitMultinomAll.mclapply <- function(reu13.df, phi, y, n, phi.new = NULL,
    coefstart = NULL){
  ### since vglm change seeds within it's call, backup first and restore later.
  if(exists(".Random.seed", envir = .GlobalEnv)){
    .GlobalEnv$.Random.seed.org <- .GlobalEnv$.Random.seed
  }

  ret <- parallel::mclapply(1:length(reu13.df),
           function(i.aa){ # i'th amino acid.
             .cubfitsEnv$my.fitMultinomOne(reu13.df[[i.aa]], phi, y[[i.aa]],
                                           n[[i.aa]], phi.new = phi.new,
                                           coefstart = coefstart[[i.aa]])
           }, mc.set.seed = FALSE, mc.preschedule = FALSE)

  if(exists(".Random.seed.org", envir = .GlobalEnv)){
    .GlobalEnv$.Random.seed <- .GlobalEnv$.Random.seed.org
  }

  ret
} # End of my.fitmultinomAll.mclapply().

### For task pull.
my.fitMultinomAll.task.pull <- function(reu13.df, phi, y, n, phi.new = NULL,
    coefstart = NULL){
  ### since vglm change seeds within it's call, backup first and restore later.
  if(exists(".Random.seed", envir = .GlobalEnv)){
    .GlobalEnv$.Random.seed.org <- .GlobalEnv$.Random.seed
  }

  ret <- pbdMPI::task.pull(1:length(reu13.df),
           function(i.aa){ # i'th amino acid.
             .cubfitsEnv$my.fitMultinomOne(reu13.df[[i.aa]], phi, y[[i.aa]],
                                           n[[i.aa]], phi.new = phi.new,
                                           coefstart = coefstart[[i.aa]])
           }, bcast = TRUE)

  if(exists(".Random.seed.org", envir = .GlobalEnv)){
    .GlobalEnv$.Random.seed <- .GlobalEnv$.Random.seed.org
  }

  ret
} # End of my.fitmultinomAll.task.pull().

### For pbdLapply.
my.fitMultinomAll.pbdLapply <- function(reu13.df, phi, y, n, phi.new = NULL,
    coefstart = NULL){
  ### since vglm change seeds within it's call, backup first and restore later.
  if(exists(".Random.seed", envir = .GlobalEnv)){
    .GlobalEnv$.Random.seed.org <- .GlobalEnv$.Random.seed
  }

  ret <- pbdMPI::pbdLapply(1:length(reu13.df),
           function(i.aa){ # i'th amino acid.
             .cubfitsEnv$my.fitMultinomOne(reu13.df[[i.aa]], phi, y[[i.aa]],
                                           n[[i.aa]], phi.new = phi.new,
                                           coefstart = coefstart[[i.aa]])
           }, pbd.mode = "spmd", bcast = TRUE)

  if(exists(".Random.seed.org", envir = .GlobalEnv)){
    .GlobalEnv$.Random.seed <- .GlobalEnv$.Random.seed.org
  }

  ret
} # End of my.fitmultinomAll.pbdLapply().
