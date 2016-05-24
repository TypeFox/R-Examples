# file:    unlibrary.R 
# author:  Dan Navarro
# contact: daniel.navarro@adelaide.edu.au
# changed: 29 June 2013

# unlibrary() removes a package from the search path. It's just a wrapper to detach() that
# allows the user to remove packages using the same syntax that they use to load it via 
# library()
unlibrary <- function(package) { 
  env.name <- deparse(substitute(package))   # allow input to drop the quote marks 
  env.name <- paste("package:", env.name, sep="")   # add the "package:" bit (TODO: this is hacky)
  detach(name = env.name, character.only = TRUE)   # now detach it
}