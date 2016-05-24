#######################################################################
##
## Function:  anchors.order()
## Authors:   Jonathan Wand <wand(at)stanford.edu>
##            and Dan Hopkins 
## 
## DESCRIPTION: Calculate frequency of vignette ordering
##              Treatment of ties/intransitivies is a user option
##
## INPUT:
##   data  : object of class anchors.data
##   ties  : one of,
##             "set" 
##             "nominal"
##             "random"
##             "mset"
##
## MODIFIED: 
##   2006-10-02: JW
##   - WAS anchors.plot.R
##   - changed to formula entry, avoid attach()
##   - standardize plots, make optional
##
##   2008-04-28: JW
##   - WAS vignette.order()
##   - data must be a anchors.data object, created usually by anchors.data()
##
#######################################################################
anchors.order <- function(formula, data,
                          ties = c("set","nominal","random","mset"),
                          subset, na.action = na.omit)
{
  mf0 <- match.call(expand.dots = FALSE)

  m <- match(c("formula","data", "subset","na.action"), names(mf0), 0)
  mf0 <- mf0[c(1, m)]
  mf0$method <- "order"

  if (is.matrix(eval(mf0$data, parent.frame()))) mf0$data <- as.data.frame(data)
  mf0[[1]] <- as.name("anchors.data")
  adata <-  eval(mf0, parent.frame())

  obj <-  anchors.order.calc( adata, ties = ties)
  return(obj)
  
}  


