#######################################################################
##
## Function: anchors.chopit()
## Author  : Jonathan Wand <wand@stanford.edu>
##           http://wand.stanford.edu
## Created :  2002-05-20
##
## DESCRIPTION: Estimate Pooled Ordered Probit model
##             aka Compound Hierarchical Ordered Probit (CHOPIT) 
## 
##
## INPUT:
##
##   data    = anchors.data object
##   options = list
## 
## NOTES:
##
##     The code assumes that your choice categories are
##     sequential integers from 1 to n.cat.
##
##     ZERO is assumed to be a missing response and is excluded from
##     the likelihood without deleting the who case.
## 
## EXAMPLE:
##
## fo <- list(self = cbind(qself1,qself2)        ~ age+country,       
##            vign = cbind(qvign1,qvign2,qvign3) ~            ,
##            tau  =                             ~ age+country)
## out<- chopit( fo, data=chopitsim )
##
## MODIFIED:
##   2008-04-28
##   - big internal changes:
##     parcelled out internals to modular functions
##     and creation of anchors.options() to unify options
##     but no user observable changes to output
##      
#######################################################################
chopit <- function(formula, data, subset, options=anchors.options(), na.action = na.omit) {

  mf0 <- match.call(expand.dots = FALSE)
  m <- match(c("formula","data", "subset",
               "na.action","min.response"), names(mf0), 0)
  mf0 <- mf0[c(1, m)]
  mf0$method <- "chopit"
  
  if (!is.null(options$na.response)) mf0$na.response <- options$na.response
  if (!is.null(options$delete)) mf0$delete <- options$delete
  if (!is.null(options$debug )) mf0$debug  <- options$debug
  if (is.matrix(eval(mf0$data, parent.frame()))) mf0$data <- as.data.frame(data)
  mf0[[1]] <- as.name("anchors.data")
  adata <-  eval(mf0, parent.frame())
  if (options$debug>0) cat("anchors.data(): complete\n")

#  obj <- anchors.chopit(adata, options=options)
  
  if (options$debug>0) cat("Entering anchors.chopit\n")
  
  options$normalize <- match.arg( options$normalize, c("self","vign","hilo","none"))
  options$optimizer <- match.arg( options$optimizer, c("optim","genoud","none"))

  ## COUNT THINGS...
  if (options$debug > 0) cat("count vars\n")
  count   <- anchors.data.count( adata , options )
  ## reconcile options
  if (options$debug > 0) cat("check options\n")
  options <- anchors.chopit.check( count, options )
  ## BUILD PARAMETERS list
  if (options$debug > 0) cat("build parms\n")
  parm    <- anchors.chopit.parm( adata, count, options )

  ## Test/dump data
  if (options$debug > 0) {
    cat("Contents of PARM list:\n")
    print(parm)
    cat("Contents of COUNT list:\n")
    print(count)
    cat("Contents of OPTIONS list:\n")
    print(options)
  }

  ## put it all together
  list.parm <- list()
  list.parm[[1]] <- parm
  
  out <- anchors.chopit.fit(data  = adata,
                            parm  = parm,
                            count = count,
                            options   = options )
              
  class(out) <- "anchors.chopit"
  return(out)
  
} ## END of anchors.chopit()

