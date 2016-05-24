# Define 'global' variables local to the package.

# module.cache is an environment that holds a list of already loaded files, and their mod times for file reloading
module.cache <- new.env()

# define a list of paths to check for lrequire'd files
paths <- c('./R_modules', './lib', '../R_modules', '../lib', gsub('\\\\', '/', path.expand('~/.R_modules')))
assign('.:module.paths',
       paths,
       envir = module.cache)

assign('.:warn.not.found', TRUE, envir = module.cache)

# exports is an empty list that can be used to deliver values back to the calling environment
exports <- list()

# module.exports is an empty variable that is delivered verbatim back to the calling environment
module.exports <- NULL

# set module.change_code to 1 to force a reload of the file next time it is lrequire-d
module.change_code <- 0
assign('.:module.change_code',
       module.change_code,
       envir = module.cache)
