###
### R routines for the R package mvmeta (c) Antonio Gasparrini 2012-2014
#
.onAttach <-
function(lib, pkg) {
#
################################################################################
#
  meta <- packageDescription("mvmeta")
  attachmsg <- paste("This is mvmeta ",meta$Version,
    ". For an overview type: help('mvmeta-package').",sep="")
  packageStartupMessage(attachmsg, domain = NULL, appendLF = TRUE)
}

#
