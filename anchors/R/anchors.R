#######################################################################
##
## Function: anchors()
## Author  : Jonathan Wand <wand(at)stanford.edu>
##           http://wand.stanford.edu
## Created :  2008-04-20
##  
## DESCRIPTION: Function for analyzing data with anchoring vignettes
##
## NOTES :      Replaces separate functions chopit, anchors, vignette.order, etc
##
## OUTPUT:      object of class anchors.$method (where $method is specified at input)
##
## INPUT:
##   formula  : a list of formula, e.g., list(self=...,vign=...,tau=...,)
##              or a a single formula, e.g., self ~ vign1 + ... + vignJ 
##   data     : matrix or dataframe 
##   method   : type of analysis, choose one
##   type     : B or C for methods that use non-parametric ranks (rank, entropy)
##   subset   : (optional) logical statement based on variables in dataset
##   na.action: usual usage
## 
## MODIFIED:  
## 
#######################################################################
anchors <- function(formula, data,
                    method = c("B","C"),
                    options=anchors.options(),
                    subset,
                    combn  = FALSE,
                    na.action = na.omit) {

  ## this allows user to use method=B or =C as shortcut
  ## instead of specifying two options: e.g., method='rank',type ='B'
  method <- match.arg(method)

  if (method %in% c("B","C")) {
    type   <- method
    method <- "rank"
  }
  
  ## Process data using anchors.data() via use of eval()
  mf0 <- match.call(expand.dots = FALSE)
  m <- match(c("formula","data", "subset",
               "na.action","min.response"), names(mf0), 0)
  mf0 <- mf0[c(1, m)]
  mf0$method <- method
  
  if (!is.null(options$na.response)) mf0$na.response <- options$na.response
  if (!is.null(options$delete)) mf0$delete <- options$delete
  if (!is.null(options$debug )) mf0$debug  <- options$debug
  if (is.matrix(eval(mf0$data, parent.frame()))) mf0$data <- as.data.frame(data)
  mf0[[1]] <- as.name("anchors.data")
  adata <-  eval(mf0, parent.frame())
  if (options$debug>0) cat("anchors.data(): complete\n")

  ## just build data
  if (method=="data") return(adata)

  ## rank/non-parametric: B or C
  obj <- anchors.rank(adata=adata, fdata=data, type=type, options=options)
  
  ## do combn (optional)
  if (method == "rank" && combn == TRUE)
      obj$combn <- anchors.combn( adata=adata, fdata=data, type=type, options=options)
  
  return( obj )
  
}

