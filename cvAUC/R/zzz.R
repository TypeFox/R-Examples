.onAttach <- function(...) {
  packageStartupMessage(' ')
  packageStartupMessage('cvAUC version: ', utils::packageDescription('cvAUC')$Version)
  packageStartupMessage('Notice to cvAUC users: Major speed improvements in version 1.1.0')
  packageStartupMessage(' ')
}
