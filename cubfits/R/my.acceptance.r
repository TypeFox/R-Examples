### Inital global storages for acceptance rate.
my.set.acceptance <- function(nSave, n.aa,
    n.p = NULL,
    n.G = NULL, n.G.pred = NULL){
  .cubfitsEnv$acceptance <- list()

  ### For acceptance rate in S/M.
  .cubfitsEnv$acceptance$b <- rep(0L, n.aa)

  ### For acceptance rate in prior.
  .cubfitsEnv$acceptance$p <- rep(0L, n.p)

  ### For acceptance rate in training.
  if(!is.null(n.G)){
    .cubfitsEnv$acceptance$phi <- rep(0L, n.G)
  }

  ### For acceptance rate in prediction.
  if(!is.null(n.G.pred)){
    .cubfitsEnv$acceptance$phi.pred <- rep(0L, n.G.pred)
  }

  invisible()
} # End of my.set.acceptance().


### Updating function to the global storages based on the variable name.
my.update.acceptance <- function(var.name, accept){
  .cubfitsEnv$acceptance[[var.name]] <-
    .cubfitsEnv$acceptance[[var.name]] + accept

  invisible()
} # End of my.update.acceptance().

my.check.acceptance <- function(var.names){
  ### Since check at the end of MCMC and curr.new was renewed.
  curr.window <- .cubfitsEnv$curr.renew - 1

  for(i.var.name in var.names){
    curr.accept <- .cubfitsEnv$adaptive[[i.var.name]][[curr.window]] /
                   .CF.AC$renew.iter
    id.accept.0 <- curr.accept == 0
    id.accept.1 <- curr.accept == 1
    accept.lower <- sum(curr.accept < .CF.AC$target.accept.lower)
    accept.upper <- sum(curr.accept > .CF.AC$target.accept.upper)

    ### Print.
    .cubfitsEnv$my.cat("- var.name: ", i.var.name, "\n", sep = "")
    .cubfitsEnv$my.cat("    ill acceptance #: none = ",
                       sum(id.accept.0), ", all = ",
                       sum(id.accept.1), "\n", sep = "")
    .cubfitsEnv$my.cat("    last renew NOT in range #: lower = ", accept.lower,
                       ", upper = ", accept.upper,
                       ", total = ", accept.lower + accept.upper,
                       "\n", sep = "")
  }

  invisible()
} # End of my.check.acceptance().

