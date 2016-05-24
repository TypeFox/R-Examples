.onAttach <- function (lib, pkg) {
	packageStartupMessage("Categorical Regression Splines (version 0.15-24)\n[vignette(\"crs_faq\") provides answers to frequently asked questions]", domain = NULL,  appendLF = TRUE)
}

.onLoad <- function (lib, pkg) {
    if(is.null(options('crs.messages')$crs.messages))
        options(crs.messages = TRUE)
    ## The package np is declared in NAMESPACE and Imports
    ## (DESCRIPTION) to support npglpreg() which eventually will
    ## migrate to the np package, at which time the following are no
    ## longer required.
    if(is.null(options('np.messages')$np.messages))
        options(np.messages = TRUE)
    if(is.null(options('np.tree')$np.tree))
        options(np.tree = FALSE)    
}

.onUnload <- function (lpath){
  library.dynam.unload("crs", libpath=lpath)
}
