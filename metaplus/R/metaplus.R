.onAttach <-
  function (libname, pkgname) 
  {
    loadmsg <- "\nNote that there are changes to the names of some functions in version 0.7.5. See NEWS.\n"
    packageStartupMessage(loadmsg, domain = NULL, appendLF = TRUE)
  }

metaplus <- function(yi,sei,mods=NULL,random="normal",
      label=switch(random,"normal"="Random Normal","t-dist"="Random t-distribution", "mixture"="Random mixture"),
      plotci=FALSE,justfit=FALSE,slab=1:length(yi),
      useAGQ=FALSE,quadpoints=21,data) {
  if (!(random %in% c("normal","t-dist","mixture"))) stop("Unknown random effect type")
  if (missing(data)) 
    data <- NULL
  if (is.null(data)) {
    data <- sys.frame(sys.parent())
  }
  else {
    if (!is.data.frame(data)) {
      data <- data.frame(data)
    }
  }
  mf <- match.call()
  mf.yi <- mf[[match("yi", names(mf))]]
  mf.sei <- mf[[match("sei", names(mf))]]
  mf.slab <- mf[[match("slab", names(mf))]]
  mf.mods <- mf[[match("mods", names(mf))]]
  yi <- eval(mf.yi, data, enclos = sys.frame(sys.parent()))
  sei <- eval(mf.sei, data, enclos = sys.frame(sys.parent()))
  if (!is.null(mf.slab)) slab <- eval(mf.slab, data, enclos = sys.frame(sys.parent()))
  mods <- eval(mf.mods, data, enclos = sys.frame(sys.parent()))
  if (!is.null(mods)) {
    mods <- as.data.frame(mods)
    if (dim(mods)[2]==1) names(mods) <- deparse(mf.mods)
  }
  df <- switch(random,
                "normal"=length(yi)-1,
                "t-dist"=length(yi)-2,
                "mixture"=length(yi)-3)
  if (!is.null(mods)) df <- df-dim(mods)[2]
  if (df<=1) stop("Insufficient studies to fit model")
  if ((df<=3) & (!justfit)) warning("Very few studies. Solution may be unstable.")
  fit <- switch(random,
                "normal"=profilenorm.metaplus(yi,sei,mods=mods,justfit=justfit,plotci=plotci,slab=slab),
                "t-dist"=profilet.metaplus(yi,sei,mods=mods,justfit=justfit,plotci=plotci,slab=slab,useAGQ,quadpoints),
                "mixture"=profilemix.metaplus(yi,sei,mods=mods,justfit=justfit,plotci=plotci,slab=slab))
  fit$label <- label
  class(fit) <- "metaplus"
  return(fit)
}
  