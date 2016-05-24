## Print the statistics summary of the chains
print.rbugs <- function(x, ...)
{
    cat("\n The MCMC simulation has ", x$n.chains, " chain(s),\n each with ", x$n.iter,
        " iterations (first ", x$n.burnin, " discarded).\n Each chain is thinned by ", x$n.thin, 
        ".\n A total of ", x$n.sims ," iterations were saved.\n\n", sep = "")
}
