sim_setup <- function(base, ...) UseMethod("sim_setup")

sim_setup.data.frame <- function(base, ..., simName = "") {
  
  simFunList <- sort_sim_fun(...)
  
  new("sim_setup", base = base, simName = simName, simFunList)
  
}

sort_sim_fun <- function(...) {
  dots <- list(...)
  if (length(dots) == 0) {
    list()
  } else {
    ordering <- sapply(dots, function(fun) fun@order) %>% order
    dots[ordering]
  }
  
}

sim_setup.sim_setup <- function(base, ..., simName = base@simName) {
  argList <- c(S3Part(base, strictS3=TRUE), list(...) , list(base = base@base, simName = simName))
  do.call(sim_setup, argList)
}
