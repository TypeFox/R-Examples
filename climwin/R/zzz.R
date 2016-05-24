.onAttach <- function(...) {
  if (!interactive() || stats::runif(1) > 0.5) return()
  
  intro <- c("For an introduction to climwin see our new vignette. \n",
             "climwin v0.1.2 has changed some parameter names and levels. See help documentation and release notes for more details."
  )
  packageStartupMessage(intro)
}
