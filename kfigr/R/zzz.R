.onLoad = function(libname, pkgname){ 
  anchorenv <- new.env(parent = getNamespace("kfigr"))
  assign("anchorenv", anchorenv, envir = getNamespace("kfigr"))
  assign("types", character(0), envir=anchorenv)
  assign("index", list(), envir=anchorenv)
  assign("history", list(), envir=anchorenv)
}

.onAttach = function(libname, pkgname){
  # default options
  evalq(opts_knit$set(kfigr.link = TRUE, kfigr.prefix = FALSE), 
       envir=getNamespace('knitr'))
  # anchor hook definition
  evalq(knit_hooks$set(anchor = hook_anchor), envir = getNamespace('knitr'))
  packageStartupMessage('knitr hook "anchor" is now available')
}

.onDetach = function(libname){
  # remove default options
  evalq(opts_knit$set(kfigr.link = NULL, kfigr.prefix = NULL), 
       envir=getNamespace('knitr'))
  # remove anchor hook definition
  evalq(knit_hooks$set(anchor = NULL), envir = getNamespace('knitr'))
  packageStartupMessage('knitr hook "anchor" has been removed')
}
