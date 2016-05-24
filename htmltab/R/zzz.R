.onAttach <- function(...) {
  if (!interactive() || stats::runif(1) > 0.1) return()
  tips <- c(
    "Due to new features and arguments introduced with v.0.6.0 code can break.\nUse suppressPackageStartupMessages to eliminate package startup messages.",
    "If you find a bug, please report it to 'https://github.com/crubba/htmltab/issues'\nUse suppressPackageStartupMessages to eliminate package startup messages."
  )
  tip <- sample(tips, 1)
  packageStartupMessage(tip)
}
