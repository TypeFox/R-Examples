#######################################################################
##
## Function: anchors.options
## Author  : Jonathan Wand <wand(at)stanford.edu>
##           http://wand.stanford.edu
## Created :  2008-04-20
##
## Unify all the options so they can be used in anchors() 
#######################################################################
anchors.options <- function(...) {
  single <- FALSE
  args <- list(...)
  if (!length(args)) 
    return( .Anchors.Options )
  else {
    if ( class(args[[1]]) == "anchors.options" ) {
      tmp <- args[[1]]
      args[[1]] <- NULL
      tmp <- replace.list(tmp , args)
    } else {
      tmp <- replace.list(.Anchors.Options, args)
    }
    class(tmp) <- "anchors.options"
    return( tmp )
  }
}

.Anchors.Options <- list(## data handling
                         delete      = "minimal",
                         na.response = c(0,NA),
                         ## vignette order
                         ties        = "set",
                         ## rank
                         rank.extras = c("minentropy","cpolr"),
                         ## chopit:
                         vign.cut    = "homo",
                         vign.var    = "hetero", # replaces single.vign.var, alt: 
                         normalize   = "self",  
                         optimizer   = "optim", # replaces use.genoud
                         start       = NULL,  # was start = list(start,labels,estimated)
                         labels      = NULL,
                         estimated   = NULL,   
                         linear      = TRUE , # was use.linear
                         analytic    = TRUE , # was use.gr
                         fitted      = FALSE, # was get.prob
                         random      = FALSE, # was do.re
                         ## cpolr
                         cpolr.method = "probit",
                         ## testing / annotation
                         verbose     = FALSE,
                         silence     = FALSE,
                         debug       = FALSE, # was do.test
                         rprof       = FALSE,  # was do.profile
                         ## genoud
                         print.level = 1,
                         domain      = 5,  # was gen.domain
                         BFGS        = TRUE,
                         wait.generations=1,
                         pop.size    = 500,
                         MemoryMatrix= TRUE,
                         max.generations = 100,
                         hess        = TRUE   ,
                         print.oob   = TRUE,
                         penalty     = -1e15  ,
                         ## opt.hstep=.Machine$double.eps^(1/4),     
                         ## list of further controls for optim...
                         optim.method= "BFGS" , # was opt.meth
                         reltol      = 1e-16,
                         trace       = 5,
                         maxit       = 500,
                         REPORT      = 1,
                         ## integration
                         int.meth    ="gh",       ## local.options: gh or ai
                         int.ghorder = 16 ,       ## gh: weights and values
                         int.gh      = NULL,     
                         int.sub     = 100,          ## ai: max subdivisions
                         int.tol     = 1e-12,        ## ai: max tolerance
                         int.range   = c(-Inf,Inf),  ## ai: range for integration
                         ## other stuff
                         digits      = 16,
                         ##
                         offset      = 0,
                         ## maximun number of combn to expand
                         max.combn   =   1e5, 
                         ## initial pararameter
                         defparm.lin.gamma.constant = 0.1,
                         defparm.lin.gamma = 0.0,
                         defparm.lin.se    = 1.0,
                         defparm.non.gamma.constant = -1.0,
                         defparm.non.gamma = 0.0,
                         defparm.non.se    = 0.0,
                         defparm.beta      = 0.0,
                         defparm.theta     = 0.0
                         )

class(.Anchors.Options) <- "anchors.options"
