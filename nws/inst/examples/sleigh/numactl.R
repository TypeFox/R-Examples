library(nws)

numactl <- function(host, envVars, options) {
  env <- unlist(envVars)
  rsleighid <- env[grep('^RSleighID=[0-9][0-9]*$', env)]
  stopifnot(length(rsleighid) == 1)
  id <- as.integer(sub('^RSleighID=', '', rsleighid))
  paste(c('env', envVars,
    'numactl', '-C', options$core[id], '-m', options$memory[id],
    '--', file.path(options$scriptDir, options$scriptName)))
}

# set options needed by the 'numactl' function
defaultSleighOptions$core <- c(0, 4, 1, 5, 2, 6, 3, 7)
defaultSleighOptions$memory <- c(0, 1, 0, 1, 0, 1, 0, 1)

s <- sleigh(verbose=TRUE, scriptExec=numactl)
print(eachElem(s, sqrt, 1:3))
