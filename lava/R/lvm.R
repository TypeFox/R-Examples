##' Initialize new latent variable model
##'
##' Function that constructs a new latent variable model object
##'
##' @aliases lvm print.lvm summary.lvm
##' @param x Vector of variable names. Optional but gives control of the
##' sequence of appearance of the variables. The argument can be given as a
##' character vector or formula, e.g. \code{~y1+y2} is equivalent to
##' \code{c("y1","y2")}. Alternatively the argument can be a formula specifying
##' a linear model.
##' @param \dots Additional arguments to be passed to the low level functions
##' @param silent Logical variable which indicates whether messages are turned
##' on/off.
##' @return Returns an object of class \code{lvm}.
##' @author Klaus K. Holst
##' @seealso \code{\link{regression}}, \code{\link{covariance}},
##' \code{\link{intercept}}, ...
##' @keywords models regression
##' @export
##' @examples
##'
##' m <- lvm() # Empty model
##' m1 <- lvm(y~x) # Simple linear regression
##' m2 <- lvm(~y1+y2) # Model with two independent variables (argument)
##' m3 <- lvm(list(c(y1,y2,y3)~u,u~x+z)) # SEM with three items
##'
lvm <- function(x=NULL, ..., silent=lava.options()$silent) {

  M <- C <- par <- fix <- numeric(); mu <- list()


  noderender <- list(
                  fill=c(),
                  shape=c(),
                  label=c()
                  )

  edgerender <- list(lty=c(),
                     lwd=c(),
                     col=c(),
                     textCol=c(),
                     est=c(),
                     arrowhead=c(),
                     dir=c(),
                     cex=c(),
                     futureinfo=list())
  graphrender <- list(recipEdges="distinct")

  graphdefault <- list(
                    "fill"="white",
                    "shape"="rectangle",
                    "label"=expression(NA),
                    "lty"=1,
                    "lwd"=1,
                    "col"="black",
                    "textCol"="black",
                    "est"=0,
                    "arrowhead"="open",
                    "dir"="forward",
                    "cex"=1.5,
                    "label"=expression(),
                    "futureinfo"=c())

  modelattr <- list(
                     randomslope=list(),
                     survival=list(),
                     parameter=list(),
                     categorical=list(),
                     distribution=list(),
                     functional=list(),
                     label=list())

  res <- list(M=M, par=par, cov=C, covpar=C, fix=fix, covfix=fix,latent=list(),
              mean=mu, index=NULL, exogenous=NA,
              constrain=list(), constrainY=list(),
              attributes=modelattr, noderender=noderender,
              edgerender=edgerender, graphrender=graphrender,
              graphdef=graphdefault)
  class(res) <- "lvm"

  myhooks <- gethook("init.hooks")
  for (f in myhooks) {
    res <- do.call(f, list(x=res))
  }

  myvar <- NULL
  if (!is.list(x)) x <- list(x,...)
  for (myvar in x) {
    if (inherits(myvar,"formula")) {
        ## if (length(getoutcome(myvar))>0) {
        ##   regression(res,...,silent=silent) <- myvar
        ## } else {
        ##   myvar <- all.vars(myvar)
        ## }
        ## regression(res,...,silent=silent) <- myvar
      regression(res,silent=silent) <- myvar
    }
    if (is.character(myvar)) {
      res <- addvar(res, myvar, silent=silent)  }
  }
  if (!is.null(myvar)) {
    index(res) <- reindex(res,zeroones=TRUE) }


  return(res)
}
