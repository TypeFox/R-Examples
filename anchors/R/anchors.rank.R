#######################################################################
##
## Function: anchors.rank()
## Author  : Jonathan Wand <wand(at)stanford.edu>
## Created :  2008-05-01
##
## INPUT:
##   adata  : object of class anchors.data
#    fdata  : original dataframe
##   type  : B or C
##
#######################################################################

anchors.B <- function(adata, fdata, options) {
  mf0 <- match.call(expand.dots = FALSE)
  mf0[[1]] <- as.name("anchors.rank")
  mf0$type <- "B"
  rank <-  eval(mf0, parent.frame())
  return(rank)
}

anchors.C <- function(adata, fdata, options) {
  mf0 <- match.call(expand.dots = FALSE)
  mf0[[1]] <- as.name("anchors.rank")
  mf0$type <- "C"
  rank <-  eval(mf0, parent.frame())
  return(rank)
}

anchors.rank <- function(adata, fdata, type, options) {

  ra <- me <- cp <- su <- NULL
  
  if (options$debug > 0) cat("anchors.rank: Get ranks\n")
  ra <- anchors.rank.type(adata, type, debug=options$debug)

  SkipExtra <- "none" %in% options$rank.extras
  
  ## BREAKING TIES: uniform and omit
  ## and get other summary information 
  if (options$debug > 0) cat("anchors.rank: Get summary\n")
  su <- summary.anchors.rank.type( ra )

  ## BREAKING TIES: minentropy
  ## calculate allocation per minimum entropy
  if (options$debug > 0) cat("anchors.rank: Get minentropy\n")  
  if (!SkipExtra && "minentropy" %in% options$rank.extras)
    me <- minimum.entropy(ra, debug=options$debug )
  
  ## BREAKING TIES: cpolr
  ## need to do some data massaging for cpolr
  if (options$debug > 0) {
    cat("anchors.rank: Get cpolr\n")
  }
  if (!SkipExtra && "cpolr" %in% options$rank.extras && !is.null(adata$formula$cpolr)) {
    fdata <- trim.data(fdata, adata)
    fdata <- insert(fdata , ra )
    if (type == "C")
        fo <- as.formula( paste("cbind(Cs,Ce) ~",as.character( adata$formula$cpolr)[2] ))
    if (type == "B") 
        fo <- as.formula( paste("cbind(Bs,Be) ~",as.character( adata$formula$cpolr)[2] ))
    if (options$debug > 0)  print(fo)
    cp <- cpolr( fo, data=fdata, method = options$cpolr.method, debug=options$debug )
  }
  
  out <- list(rank=ra, summary = su, minentropy = me, cpolr = cp, type=type)
  
  class(out) <- "anchors.rank"
  return(out)
}
