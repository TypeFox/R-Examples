setClass("prediction",
         representation(predictions = "list",
                        labels      = "list",
                        cutoffs     = "list",
                        fp          = "list",
                        tp          = "list",
                        tn          = "list",
                        fn          = "list",
                        n.pos       = "list",
                        n.neg       = "list",
                        n.pos.pred  = "list",
                        n.neg.pred  = "list"))

setClass("performance",
         representation(x.name       = "character",
                        y.name       = "character",
                        alpha.name   = "character",
                        x.values     = "list",
                        y.values     = "list",
                        alpha.values = "list" ))

#setMethod("plot",signature(x="performance",y="missing"),
#          function(x,y,...) {
#              .plot.performance(x,...)
#          })

setMethod("plot",signature(x="performance",y="missing"),
          function(x,y,..., avg="none", spread.estimate="none",
  spread.scale=1, show.spread.at=c(), colorize=FALSE,
  colorize.palette=rev(rainbow(256,start=0, end=4/6)),
  colorkey=colorize, colorkey.relwidth=0.25, colorkey.pos="right",
  print.cutoffs.at=c(), cutoff.label.function=function(x) { round(x,2) },
  downsampling=0, add=FALSE ) {

              .plot.performance(x,..., 
              					avg= avg, 
              					spread.estimate= spread.estimate,
								spread.scale= spread.scale, 
								show.spread.at= show.spread.at, 
								colorize= colorize,
								colorize.palette= colorize.palette,
								colorkey= colorkey, 
								colorkey.relwidth= colorkey.relwidth, 
								colorkey.pos= colorkey.pos,
								print.cutoffs.at= print.cutoffs.at, 
								cutoff.label.function= cutoff.label.function,
								downsampling= downsampling, 
								add= add)
          })


## .First.lib <- function( libname, pkgname, where) {
##     if (!require(methods)) {
##         stop("Require Methods package")
##     }
##     if (!require(gplots)) {
##         stop("Require gplots package")
##     }
    
##     where <- match(paste("package:",pkgname, sep=""), search())
## }


