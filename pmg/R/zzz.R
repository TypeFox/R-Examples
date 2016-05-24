
## when we have a namespace, we can do this.
.onAttach <- function(...){
    packageStartupMessage("Loading pmg...")
  pmg()
}
