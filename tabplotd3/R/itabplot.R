.e <- new.env()

# system.file <- function(..., pkg=".", package = "base", lib.loc = NULL, mustWork = FALSE){
#   pkg <- as.package(pkg)
#   if (package == pkg$package){
#     file.path(normalizePath(pkg$path, winslash="/"), "inst", ...)
#   } else {
#     base::system.file(..., package=package, lib.loc=lib.loc, mustWork=mustWork)
#   }
# }

#' Interactive tableplot
#' 
#' \code{itabplot} is an interactive version of \code{\link{tableplot}}. It starts 
#' your browser and allows for zooming, panning and sorting the tableplot. This
#' version can be used for explorative usage, while \code{tableplot} can be used for
#' publication purposes.
#' It needs the same parameters as tabplot.
#' @param x \code{data.frame} or \code{ffdf} object used for plotting.
#' @param ... parameters that will be given to \code{tableplot}
#' @seealso \code{\link{tableplot}}
#' @export
#' @importFrom httpuv runServer
#' @import Rook
#' @importFrom brew brew
#' @importFrom tabplot tableplot
itabplot <- function(x, ...){
  xlit <- deparse(substitute(x))
  tp <- tableplot(x, plot=FALSE, ...)
  args <- list(...)
  
  
  app <- Builder$new( Static$new( urls = c('/css/','/img/','/js/','/.+\\.json$') #, "/.*\\.html$")
                                , root = system.file("app", package="tabplotd3")#"."
                                )
                    , Brewery$new( url='.*\\.html$'
                                 , root= system.file("app", package="tabplotd3")
                                 , dat=x
                                 , xlit=xlit
                                 , args=args
                                 , ...
                                 )
                    , URLMap$new( '^/json' = tpjson
                                  #, ".*" = Redirect$new("/tableplot.html")
                                  , ".*" = function(env){
                                             req <- Request$new(env)
                                             body <- paste(capture.output(brew(system.file("app/index.html", package="tabplotd3"))),collapse="\n")
                                             res <- Response$new()
                                             res$write(body)
                                             res$finish()
                                           }
                                )
                    )

  browseURL("http://localhost:8100")
  id <- runServer("0.0.0.0", 8100, list(call=app$call, onWSOpen=function(ws){}, onHeaders=function(x){}))
}

#### testing
### just run this whole file after load_all()
# irisg <- iris[sample(nrow(iris), 1e4, replace=TRUE),]
#itabplot(iris)
